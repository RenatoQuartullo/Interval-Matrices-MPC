clear
clc
close all

addpath tools
addpath data

% import maximal robust invariant sets computed from [1]
load 'data/RIS_coverage.mat'

% System
Ahat = [1 0.15; 0.1 1];
Bhat = [0.1; 1.1];
[n,m] = size(Bhat);

% MPC parameters
Nmax = 10;                                             % N: prediction horizon, j = 0...N
gamma = 1;

% State/input constraints
xmax = [8; 8];
umax = 4;
X = cell(Nmax+1,1);
U = cell(Nmax+1,1);
for i = 1:Nmax+1
    X{i}.A = [eye(n); -eye(n)];
    X{i}.b = [xmax; xmax];
    U{i}.A = [eye(m); -eye(m)];
    U{i}.b = [umax; umax];
end
F1 = [diag(1./xmax); -diag(1./xmax)];
G1 = [diag(1./umax); -diag(1./umax)];
F = [F1; zeros(size(G1,1),n)];
G = [zeros(size(F1,1),m); G1];

% possile values of \delta
delta_A = 0.05:0.05:0.6;

%% Begin calculating coverage
coverage_IMMPC = zeros(length(delta_A),1);
coverage_PTMPC_N3 = zeros(length(delta_A),1);
coverage_PTMPC_N10 = zeros(length(delta_A),1);
oPt.Display = 'off';

for ii = 1:length(delta_A)

    % Parametric uncertainty
    Delta = [delta_A(ii) 0 0; 0 0 0.1];
    IS = intervalMatrix([Ahat Bhat], Delta);
    IS_vertices = compute_vertices_matrix(IS);
    A_vertices = IS_vertices(:,1:n,:);
    B_vertices = IS_vertices(:,1+n:end,:);

    % Calculate K
%     K = KStableVertices(A_vertices,B_vertices);
    K = -place(Ahat,Bhat,[0.1 0.3]);
    ABK_vertices = pagemtimes(IS_vertices,[eye(n);K]);
    AKhat = Ahat + Bhat*K;

    % Terminal Set Xf
    Xf = RIS{ii};

    % Gris initial conditions
    x0 = RIS{ii}.grid(10)';
    nx0 = size(x0,2);

    %% IM-MPC
    IDelta = intervalMatrix(zeros(n,m+n),Delta);
    IDK = IS*[eye(n); K];
    [Bc,Bsup,BKc,BKsup] = PreComputeIntervalMatrices(IDK,IDelta,K,Nmax);
    for i = 1:Nmax
        [Hineq{i},gineq{i},Heq{i},geq{i},geqDyn{i},cJ{i},LB{i},UB{i}] = OPmatrices(Ahat, Bhat, i, Bc, Bsup, BKc, BKsup, X, U, Xf, K, 'LP');
    end

    feasibleX0 = true(nx0,1);
    for i = 1:nx0
        disp([num2str(ii), '    ' ,num2str(i)])

        % Solve IM-MPC problem
        [sol,JN] = SolveIMMPC(x0(:,i), Nmax, Hineq, gineq, Heq, geq, geqDyn, cJ, LB, UB, gamma, 'LP', oPt);

        % Check initial feasibility
        if min(JN) == inf
            feasibleX0(i) = false;
        end
    end
    x0feasible = x0(:,feasibleX0);
    coverage_IMMPC(ii) = size(x0feasible,2)/nx0;

    %% PT-MPC
    Hc = ComputeHc_PTMPC(Xf.A,F,G,K);
    H = ComputeH_PTMPC(Xf.A,ABK_vertices);

    % horizon length N = 10
    N_PTMPC = 10;
    [Hineq_PTMPC,gineq_PTMPC,galpha0_PTMPC,Heq_PTMPC,geq_PTMPC] = OPmatrices_PTMPC(N_PTMPC, H, Hc, Xf.A, B_vertices, G, [], Xf.b);
    feasibleX0PTMPC = true(nx0,1);
    for i = 1:nx0

        % Solve PT-MPC problem
        gineq_k = [galpha0_PTMPC*x0(:,i); gineq_PTMPC];
        solPTMPC = linprog([],Hineq_PTMPC,gineq_k,Heq_PTMPC,geq_PTMPC,[],[],oPt);
        if isempty(solPTMPC)
            feasibleX0PTMPC(i) = false;
        end
    end
    x0feasiblePTMPC = x0(:,feasibleX0PTMPC);
    coverage_PTMPC_N10(ii) = size(x0feasiblePTMPC,2)/nx0;

    % N = 3
    N_PTMPC = 3;
    [Hineq_PTMPC,gineq_PTMPC,galpha0_PTMPC,Heq_PTMPC,geq_PTMPC] = OPmatrices_PTMPC(N_PTMPC, H, Hc, Xf.A, B_vertices, G, [], Xf.b);
    feasibleX0PTMPC = true(nx0,1);
    for i = 1:nx0

        % Solve PT-MPC problem
        gineq_k = [galpha0_PTMPC*x0(:,i); gineq_PTMPC];
        solPTMPC = linprog([],Hineq_PTMPC,gineq_k,Heq_PTMPC,geq_PTMPC,[],[],oPt);
        if isempty(solPTMPC)
            feasibleX0PTMPC(i) = false;
        end
    end
    x0feasiblePTMPC = x0(:,feasibleX0PTMPC);
    coverage_PTMPC_N3(ii) = size(x0feasiblePTMPC,2)/nx0;

end
% save('data/coverage.mat','coverage_PTMPC_N3','coverage_PTMPC_N10','delta_A','coverage_IMMPC')


%% Plot and Results
% To see coverage results in the paper -> load data/coverage.mat
figure
plot(delta_A,coverage_IMMPC,'-*','LineWidth',1.5,'MarkerSize',10)
hold on
load data/coverage_SLSMPC.mat
plot(delta_A,coverage_SLSMPC_N3,'s-','Color',rgb('SeaGreen'),'LineWidth',1.5)
plot(delta_A,coverage_SLSMPC_N10,'s--','Color',rgb('SeaGreen'),'LineWidth',1.5)
load data/coverage_OTMPC.mat
plot(delta_A,coverage_OTMPC,'s-','Color',rgb('Gold'),'LineWidth',1.5)
plot(delta_A,coverage_PTMPC_N3,'s-','Color',rgb('Crimson'),'LineWidth',1.5)
plot(delta_A,coverage_PTMPC_N10,'s--','Color',rgb('Crimson'),'LineWidth',1.5)
legend('IM-MPC','SLS-MPC (N=3)','SLS-MPC (N=10)','OT-MPC (N=3)','PT-MPC (N = 3)','PT-MPC (N = 10)','interpreter','latex')
xlabel('$\delta$','Interpreter','latex')
xlim([delta_A(1) delta_A(end)])
grid on

% [1] S. Chen, V. M. Preciado, M. Morari, and N. Matni, "Robust model predictive control with polytopic model
% uncertainty through system level synthesis," Automatica, vol. 162, p. 111431, 2024.
