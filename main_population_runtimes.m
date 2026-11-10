clear 
clc
close all

addpath tools
addpath data

% System (Fertility and survivance) parameters
n = 6;
f = [0.01 0.45 0.4 0.14 0 0];
s = [0.9 0.95 0.9 0.85 0.8];
Ahat = [f; [diag(s), zeros(n-1,1)] ];
Ahat(n,n) = 0.2;
Bhat = eye(n);
Bhat = Bhat(:,[2 3 4 5]);
[~,m] = size(Bhat);
xref = [1 1 1 1 1 1]';
uref = pinv(Bhat)*(eye(n)-Ahat)*xref;
xmax = 1.5*[1 1 1 1 1 1]';
umax = [0.35 0.35 0.3 0.25]'-uref;
Nmax = 10;                                                       
gamma = 1;
Q = eye(n);
R = eye(m);
N_PTMPC = 10;

% Uncertainty
DeltaAll = zeros(n,n+m);
DeltaAll(1,1:n) = 0.08*f;
DeltaAll(2:n,1:n) = 0.03*Ahat(2:end,1:end);
DeltaAll(:,n+1:end) = 0.01*Bhat;
Delta = zeros(n,n+m);

% State/input constraints
F1 = [eye(n); -eye(n)];
G1 = [eye(m); -eye(m)];
f1 = [xmax; xmax];
g1 = [umax; umax];
F = [F1; zeros(size(G1,1),n)];
G = [zeros(size(F1,1),m); G1];
h = [f1; g1];
X = cell(Nmax+1,1);
U = cell(Nmax+1,1);
for i = 1:Nmax+1
    X{i}.A = F1;
    X{i}.b = f1;
    U{i}.A = G1;
    U{i}.b = g1;
end

%% Monte Carlo for increasing vertices number and for multiple initial conditions
nx0 = 20;
x0 = repmat((-1:0.1:1),[n,1]);
posRand_not0 = find(DeltaAll ~= 0);
n_unc = length(posRand_not0);
toc_PTMPC = nan(n_unc,nx0);
toc_IMMPC = nan(n_unc,nx0);
n_ineq_IMMPC = nan(n_unc,1);
n_eq_IMMPC = nan(n_unc,1);
n_ineq_PTMPC = nan(n_unc,1);
n_eq_PTMPC = nan(n_unc,1);
N = zeros(n_unc,nx0);
N_QP = zeros(n_unc,nx0);
oPt.Display = 'off';

for iu = 1:n_unc
    
    disp(iu)
    
    % Add uncertainty entry one by one
    posRand = posRand_not0(iu);
    Delta(posRand) = DeltaAll(posRand);  

    % Calculate vertices
    IS = intervalMatrix([Ahat Bhat],Delta);     % IS: interval matrix set for system matrices
    IS_vertices = compute_vertices_matrix(IS);
    A_vertices = IS_vertices(:,1:n,:);
    B_vertices = IS_vertices(:,1+n:end,:);
    
    % Calculate K for pole placement
    K = -place(Ahat,Bhat,[0.3 0.3 0.3 0.3 0.2 0.2]);
    
    % Compute vertices of A+BK
    ABK_vertices = pagemtimes(IS_vertices,[eye(n); K]);
    AKhat = Ahat + Bhat*K;
    % if ~CheckLMI(ABK_vertices)
    %     disp('K not satisfying LMI')
    %     break
    % end

    % Compute terminal set
    [W,~] = schur(AKhat,'real');
    V = inv(W);
    alpha_N = FindAlphaRPIParallelotope(ABK_vertices,V,W,F,G,K,h,'non-symmetric');
    if any(isnan(alpha_N))
        disp('alpha_N not found')
        break
    end
    Xf.A = [V;-V];
    Xf.b = [alpha_N(1:n); alpha_N(n+1:end)];  
    
    %% PT-MPC
    
    % Terminal weight
    E = [eye(m), zeros(m,(N_PTMPC-1)*m)];
    M = [zeros(m*N_PTMPC,m), [eye(m*(N_PTMPC-1)); zeros(m,m*(N_PTMPC-1))]];
    Psi = [AKhat, Bhat*E; zeros(m*N_PTMPC,n) M];
    Qhat = [Q+K'*R*K, K'*R*E; E'*R*K, E'*R*E];
    Wp = dlyap(Psi',Qhat);

    % Pre-compute all needed matrices
    Hc = ComputeHc_PTMPC(Xf.A,F,G,K);
    H = ComputeH_PTMPC(Xf.A,ABK_vertices);
    [Hineq_PTMPC,gineqPTMPC,galpha0PTMPC,Heq_PTMPC,geq_PTMPC,Qqp_PTMPC] = OPmatrices_PTMPC(N_PTMPC, H, Hc, Xf.A, B_vertices, G, Wp, Xf.b);

    % Calculate solution for multiple x0
    for ix = 1:nx0
        x_PTMPC = x0(:,ix);
        gineq_PTMPC_k = [galpha0PTMPC*x0(:,i); gineqPTMPC];
        tic_PTMPC = tic;
        [~,~,flag_PTMPC] = quadprog((Qqp_PTMPC+Qqp_PTMPC')/2, zeros(size(Qqp_PTMPC,2),1),Hineq_PTMPC,gineq_PTMPC_k,Heq_PTMPC,geq_PTMPC,[],[],[],oPt);
        toc_PTMPC(iu,ix) = toc(tic_PTMPC);
        if flag_PTMPC ~= 1
            disp('infeasible PT-MPC')
            break
        end       
    end

    % Calculate number of in/equalities of OCP
    n_ineq_PTMPC(iu) = size(Hineq_PTMPC,1);
    n_eq_PTMPC(iu) = size(Heq_PTMPC,1);
    
    %% IM-MPC
    ID = intervalMatrix(zeros(n,m+n),Delta);
    IDK = ID*[eye(n); K];
    [Bc,Bsup,BKc,BKsup] = PreComputeIntervalMatrices(IDK,ID,K,Nmax);
    for i = 1:Nmax
        [Hineq{i},gineq{i},Heq{i},geq{i},geqDyn{i},cJ{i},LB{i},UB{i}] = OPmatrices(Ahat,Bhat,i,Bc,Bsup,BKc,BKsup,X,U,Xf,K,'QP');
    end

    % MPC Iteration
    for ix = 1:nx0
        tic_IMMPC = tic;
        [sol,JN] = SolveIMMPC(x0(:,ix),Nmax,Hineq,gineq,Heq,geq,geqDyn,cJ,LB,UB,gamma,'QP',oPt);
        [J,N(iu,ix)] = min(JN);
        if J == inf
            disp('infeasible IM-MPC')
            break
        end
        toc_IMMPC(iu,ix) = toc(tic_IMMPC);
    end

    % Calculate number of in/equalities of OCP (calculated for the maximum N)
    n_ineq_IMMPC(iu) = size(Hineq{Nmax},1);
    n_eq_IMMPC(iu) = size(Heq{Nmax},1);
    
end
% save('data/population_runtimes.mat','toc_IMMPC','toc_PTMPC','n_ineq_IMMPC','n_eq_IMMPC','n_ineq_PTMPC','n_ineq_PTMPC','n_unc')

%% PLOT
% Too see the paper results --> load data/population_runtimes.mat
figure
hold on
plot(1:n_unc, mean(toc_IMMPC,2),'-*','LineWidth',1.5)
load data/population_runtimes_SLSMPC.mat
plot(1:n_unc, toc_SLSMPC,'-s','Color',rgb('SeaGreen'),'LineWidth',1.5)
plot(1:n_unc, mean(toc_PTMPC,2),'-s','Color',rgb('Crimson'),'LineWidth',1.5)
xlabel('Number  of uncertain entries of $A$','Interpreter','latex')
xlim([1 n_unc])
ylabel('Average CPU time [s]', 'Interpreter','latex')
grid on
legend('IM-MPC','SLS-MPC','PT-MPC', 'Interpreter','latex')
box on

figure
hold on
plot(1:n_unc, n_ineq_IMMPC+n_eq_IMMPC,'-*','LineWidth',1.5)
plot(1:n_unc, ncostr_SLSMPC,'-s','Color',rgb('SeaGreen'),'LineWidth',1.5)
plot(1:n_unc, n_ineq_PTMPC+n_eq_PTMPC,'-s','Color',rgb('Crimson'),'LineWidth',1.5)
xlabel('Number  of uncertain entries of $[A \quad B]$','Interpreter','latex')
xlim([1 n_unc])
ylabel('Numeber of constraints', 'Interpreter','latex')
grid on
legend('IM-MPC','SLS-MPC','PT-MPC', 'Interpreter','latex')
box on
set(gca, 'YScale', 'log')
