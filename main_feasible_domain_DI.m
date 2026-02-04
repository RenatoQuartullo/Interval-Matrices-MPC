clear
clc
close all

addpath data
addpath tools

% Double Integrator System
Ahat = [1 1; 0 1];
Bhat = [0 1]';
[n,m] = size(Bhat);

% MPC parameters
Nmax = 25;                                  % Nmax: Maximum horizon length for varaible-horizon enumeration method                                           
gamma = 1;                                  % gamma: horizon length weighting parameter

% State/input constraints
xmax = [12; 4];
umax = 2;
X = cell(Nmax+1,1);
U = cell(Nmax+1,1);
for i = 1:Nmax+1
    X{i}.A = [eye(n); -eye(n)];
    X{i}.b = [xmax; xmax];
    U{i}.A = [eye(m); -eye(m)];
    U{i}.b = [umax; umax];
end

% Gris initial conditions
x01 = linspace(-xmax(1), xmax(1), 30);
x02 = linspace(-xmax(2), xmax(2), 10);
x0 = [];
for i = 1:size(x01,2)
    for j = 1:size(x02,2)
        x0 = [x0 [x01(i);x02(j)] ];
    end
end
nx0 = size(x0,2);

% Parametric uncertainty
Delta = [0.1 0.05 0.05; 0.01 0.03 0.02];
IS = intervalMatrix([Ahat Bhat],Delta);     % IS: interval matrix set for system matrices
IS_vertices = compute_vertices_matrix(IS);
A_vertices = IS_vertices(:,1:n,:);
B_vertices = IS_vertices(:,1+n:end,:);

% Calculate K
% K = KStableVertices(VertAM, VertBM);
K = [-0.47 -1.48];
ABK_vertices = pagemtimes(IS_vertices,[eye(n);K]);
AKhat = Ahat + Bhat*K;

% Calculate RPI terminal set 
F1 = [diag(1./xmax); -diag(1./xmax)];
G1 = [diag(1./umax); -diag(1./umax)];
F = [F1; zeros(size(G1,1),n)];
G = [zeros(size(F1,1),m); G1];
h = ones(size(F,1),1);
[W,~] = eig(AKhat);
V = inv(W);
alpha_N = FindAlphaRPIParallelotope(ABK_vertices,V,W,F,G,K,h,'symmetric');
if isnan(alpha_N)
    disp('alpha_N not found')
    return
end
Xf.A = [V;-V];
Xf.b = [alpha_N(1:n); alpha_N(n+1:end)];

%% INTERVAL MATRIX MPC (IM-MPC)

% Pre compute all interval matrices needed for optimization problem
IDelta = intervalMatrix(zeros(n,m+n),Delta);
IDK = IS*[eye(n); K];
[Bc,Bsup,BKc,BKsup] = PreComputeIntervalMatrices(IDK,IDelta,K,Nmax);
for i = 1:Nmax
    [Hineq{i},gineq{i},Heq{i},geq{i},geqDyn{i},cJ{i},LB{i},UB{i}] = OPmatrices(Ahat, Bhat, i, Bc, Bsup, BKc, BKsup, X, U, Xf, K, 'LP');
end

% Calculate feasible domain 
N = nan(nx0,1);                             % N collect all optimal horizon lengths
toc_IMMPC = zeros(nx0,1);
feasibleX0 = true(nx0,1);
oPt.Display = 'off';
for i = 1:nx0

    % Solve IM-MPC problem
    tic_OP = tic;
    [sol,JN] = SolveIMMPC(x0(:,i), Nmax, Hineq, gineq, Heq, geq, geqDyn, cJ, LB, UB, gamma, 'LP', oPt);
    toc_IMMPC(i) = toc(tic_OP);
    [~,N(i)] = min(JN);

    % Check initial feasibility
    if min(JN) == inf
        feasibleX0(i) = false;
    end
end

x0feasible = x0(:,feasibleX0);
FD_IMMPC = Polyhedron('V', x0feasible');          % Calculate feasible domain as convex hull of feeasible initial conditions

%% Polytopic Tube MPC (PT-IMPC)

% Pre-compute all needed matrices for N = 1,2,...,Nmax
Q = eye(n);
R = eye(m);
F1 = [diag(1./xmax); -diag(1./xmax)];
G1 = [diag(1./umax); -diag(1./umax)];
F = [F1; zeros(size(G1,1),n)];
G = [zeros(size(F1,1),m); G1];
h = ones(size(F,1),1);
Hc = ComputeHc_PTMPC(Xf.A,F,G,K);
H = ComputeH_PTMPC(Xf.A,ABK_vertices);
for i = 1:Nmax
    E = [eye(m), zeros(m,(i-1)*m)];
    M = [zeros(m*i,m), [eye(m*(i-1)); zeros(m,m*(i-1))]];
    Psi = [AKhat, Bhat*E; zeros(m*i,n) M];
    Qhat = [Q+K'*R*K, K'*R*E; E'*R*K, E'*R*E];
    Wp = dlyap(Psi',Qhat);
    [HineqPTMPC{i},gineqPTMPC{i},galpha0PTMPC{i},HeqPTMPC{i},geqPTMPC{i},QqpPTMPC{i}] = OPmatrices_PTMPC(i, H, Hc, Xf.A, B_vertices, G, Wp, Xf.b);
end

% Calculate feasible domain 
feasibleX0PTMPC = true(nx0,Nmax);
toc_PTMPC = nan(nx0,Nmax);  
NminFeasible = Nmax*ones(nx0,1);
for i = 1:nx0

    for j = 1:Nmax

        % Solve PT-MPC problem
        gineqPTMPC_k = [galpha0PTMPC{j}*x0(:,i); gineqPTMPC{j}];
        tic_PTMPC = tic;
%         [~,~,flag] = quadprog(QqpPTMPC{j}, [], HineqPTMPC{j}, gineqPTMPC_k, HeqPTMPC{j}, geqPTMPC{j},[],[],[],oPt);
        [~,~,flag] = linprog([], HineqPTMPC{j}, gineqPTMPC_k, HeqPTMPC{j}, geqPTMPC{j},[],[],oPt);
        toc_PTMPC(i,j) = toc(tic_PTMPC);

        % Check initial feasibility
        if flag ~= 1
            feasibleX0PTMPC(i,j) = false;
        else
            NminFeasible(i) = j;
            break
        end
        
    end
end
feasICPTMPC = sum(feasibleX0PTMPC,2) > 0;
x0feasiblePTMPC = x0(:,feasICPTMPC);
FD_PTMPC = Polyhedron('V',x0feasiblePTMPC');
toc_PTMPC_Nmin = zeros(nx0,1);
for i = 1:nx0
    toc_PTMPC_Nmin(i) = toc_PTMPC(i,NminFeasible(i));
end

%% Plot and Results
% To see paper results --> load data/FD.mat
figure
hold on
load data/FD_SLSMPC.mat
plot(FD_IMMPC, 'alpha',0.2, 'color','blue','EdgeColor','b')
plot(FD_SLSMPC, 'alpha',0.2, 'color','green','LineStyle','--','EdgeColor','g')
% plot(Polyhedron(Xf.A,Xf.b),'color','green')
plot(FD_PTMPC, 'alpha',0.3, 'color','red','EdgeColor','r')
load data/FD_OTMPC.mat
plot(FD_OTMPC, 'alpha',0.5, 'color','yellow','LineStyle','--','EdgeColor','yellow')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
legend('IM-MPC','SLS-MPC','PT-MPC','OT-MPC (N = 3)','interpreter','latex')
axis([-12.5 20 -4.5 4.5])
box on

