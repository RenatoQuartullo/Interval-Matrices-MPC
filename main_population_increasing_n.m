clear 
% clc
close all

addpath tools
addpath data

oPt.Display = 'off';

% System (Fertility and survivance) parameters
nstate = [6 12 18 24];

f{1} = [0.01 0.45 0.4 0.14 0 0];
f{2} = repelem(f{1}/2, 2);
f{3} = repelem(f{1}/3, 3);
f{4} = repelem(f{1}/4, 4);

s{1} = [0.9 0.95 0.9 0.85 0.8 0.2];
s{2} = [0.9 0.9 0.95 0.95 0.9 0.9 0.85 0.85 0.8 0.8 0.8 0.2];
s{3} = [0.9 0.9 0.9 0.95 0.95 0.95 0.9 0.9 0.9 0.85 0.85 0.85 0.8 0.8 0.8 0.8 0.8 0.2];
s{4} = [0.9 0.9 0.9 0.9 0.95 0.95 0.95 0.95 0.9 0.9 0.9 0.9 0.85 0.85 0.85 0.85 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.2];

BhatAll{1} = [zeros(1,4); eye(4); zeros(1,4)];
BhatAll{2} = [zeros(1,10); eye(10); zeros(1,10)];
BhatAll{3} = [zeros(1,16); eye(16); zeros(1,16)];
BhatAll{4} = [zeros(1,22); eye(22); zeros(1,22)];

polesAK{1} = [0.3*ones(1,4) 0.2*ones(1,2)];
polesAK{2} = [0.3*ones(1,8) 0.2*ones(1,4)];
polesAK{3} = [0.3*ones(1,12) 0.2*ones(1,6)];
polesAK{4} = [0.3*ones(1,16) 0.2*ones(1,8)];

nexp = length(nstate);
N_PTMPC = 10;
Nmax = 10;                                                       
gamma = 1;
nx0 = 20;

N = nan(nexp,nx0);
toc_IMMPC = nan(nexp,nx0);
toc_PTMPC = nan(nexp,nx0);

for in = 1:nexp

    disp(in)

    n = nstate(in);

    Ahat = [f{in}; [diag(s{in}(1:n-1)), zeros(n-1,1)] ];
    Ahat(n,n) = s{in}(n);
    Bhat = BhatAll{in};
    [~,m] = size(Bhat);
    xref = ones(n,1);
    uref = pinv(Bhat)*(eye(n)-Ahat)*xref;
    xmax = 1.5*ones(n,1)*10000;
    umax = 0.3*ones(size(uref))*1000-uref;
    x0 = repmat((-1:2/nx0:1),[n,1]);
    Q = eye(n);
    R = eye(m);


    % Uncertainty
    
    DeltaA = zeros(n);
    DeltaB = zeros(n,m);
    DeltaA(1,1:2) = 0.0*f{in}(1:2);
    DeltaA(end,end-1:end) = 0.0*s{in}(end-1:end);
    Delta = [DeltaA DeltaB];

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

    % Calculate vertices
    IS = intervalMatrix([Ahat Bhat],Delta);     % IS: interval matrix set for system matrices
    IS_vertices = compute_vertices_matrix(IS);
    A_vertices = IS_vertices(:,1:n,:);
    B_vertices = IS_vertices(:,1+n:end,:);

    % Calculate K for pole placement
    K = -place(Ahat,Bhat,polesAK{in});

    % Compute vertices of A+BK
    ABK_vertices = pagemtimes(IS_vertices,[eye(n); K]);
    AKhat = Ahat + Bhat*K;
    if ~CheckLMI(ABK_vertices)
        disp('K not satisfying LMI')
        break
    end

    % Compute terminal set (low-complexity) 
    [W,~] = schur(AKhat,'real');
    V = inv(W);
    alpha_N = FindAlphaRPIParallelotope(ABK_vertices,V,W,F,G,K,h,'non-symmetric');
    if any(isnan(alpha_N))
        disp('alpha_N not found')
        break
    end
    Xf.A = [V;-V];
    Xf.b = [alpha_N(1:n); alpha_N(n+1:end)];  

    %% PT-MPC (low complexity - parallelotopic shape) -> quadprog of MOSEK SUGGESTED

    % Terminal weight
    E = [eye(m), zeros(m,(N_PTMPC-1)*m)];
    M = [zeros(m*N_PTMPC,m), [eye(m*(N_PTMPC-1)); zeros(m,m*(N_PTMPC-1))]];
    Psi = [AKhat, Bhat*E; zeros(m*N_PTMPC,n) M];
    Qhat = [Q+K'*R*K, K'*R*E; E'*R*K, E'*R*E];
    % Wp = dlyap(Psi',Qhat);
    
    %%%%%%%%%%%%%%%%% Alternative to dlyap %%%%%%%%%%%%%%%%%
    Qhat_sym = (Qhat + Qhat')/2;   
    Wp = zeros(size(Qhat_sym));
    P = eye(size(Psi));   
    max_iter = 100;  
    for k = 1:max_iter
        Wp = Wp + P' * Qhat_sym * P;
        P = P * Psi;
    end
    Wp = (Wp + Wp')/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Pre-compute all needed matrices
    nv = size(IS_vertices,3);
    Phi_tilde = zeros(n,n,nv);
    Phi_tilde_m = zeros(n,n,nv);
    Phi_tilde_p = zeros(n,n,nv);
    B_tilde = zeros(n,m,nv);
    for i = 1:nv
        Phi_tilde(:,:,i) = V*ABK_vertices(:,:,i)*W;
        Phi_tilde_p(:,:,i) = max(Phi_tilde(:,:,i),0);
        Phi_tilde_m(:,:,i) = max(-Phi_tilde(:,:,i),0);
        B_tilde(:,:,i) = V*B_vertices(:,:,i);
    end
    F_tilde = F*W;
    K_tilde = K*W;
    FGK_p = max(F_tilde+G*K_tilde, 0);
    FGK_m = max(-(F_tilde+G*K_tilde), 0);
    [Hineq_PTMPC,gineq_PTMPC,Heq_PTMPC,geq_PTMPC,Qqp_PTMPC] = OPmatrices_PTMPC_lowComplexity(N_PTMPC, Phi_tilde_p, Phi_tilde_m, B_tilde, FGK_p, FGK_m, G, h, V, alpha_N, Wp);

    % Calculate solution for multiple x0
    % addpath('/Users/renatoquartullo/mosek/11.1/toolbox/r2022b', '-begin')
    for ix = 1:nx0
        x_PTMPC = x0(:,ix);
        gineq_PTMPC_k = [gineq_PTMPC; V*x0(:,ix); -V*x0(:,ix)];
        tic_PTMPC = tic;
        [~,~,flag_PTMPC] = quadprog((Qqp_PTMPC+Qqp_PTMPC')/2,zeros(size(Qqp_PTMPC,2),1),Hineq_PTMPC,gineq_PTMPC_k,Heq_PTMPC,geq_PTMPC,[],[],[],oPt);
        toc_PTMPC(in,ix) = toc(tic_PTMPC);
        if flag_PTMPC ~= 1
            disp('infeasible PT-MPC')
            break
        end       
    end



    %% IM-MPC
    ID = intervalMatrix(zeros(n,m+n),Delta);
    IDK = ID*[eye(n); K];
    [Bc,Bsup,BKc,BKsup] = PreComputeIntervalMatrices(IDK,ID,K,Nmax);
    for i = 1:Nmax
        [Hineq{i},gineq{i},Heq{i},geq{i},geqDyn{i},cJ{i},LB{i},UB{i}] = OPmatrices(Ahat,Bhat,i,Bc,Bsup,BKc,BKsup,X,U,Xf,K,'QP');
    end

    % addpath('/Applications/MATLAB_R2025b.app/toolbox/optim/optim/', '-begin')
    % MPC Iteration
    for ix = 1:nx0
        tic_IMMPC = tic;
        [sol,JN] = SolveIMMPC(x0(:,ix),Nmax,Hineq,gineq,Heq,geq,geqDyn,cJ,LB,UB,gamma,'QP',oPt);
        [J,N(in,ix)] = min(JN);
        if J == inf
            disp('infeasible IM-MPC')
            break
        end
        toc_IMMPC(in,ix) = toc(tic_IMMPC);
    end
    
end
% save('data/population_runtimes.mat','toc_IMMPC','toc_PTMPC','n_unc')

%% PLOT
% Too see the paper results --> load data/population_runtimes.mat
figure
hold on
plot(nstate, mean(toc_IMMPC,2),'-*','LineWidth',1.5)
plot(nstate, mean(toc_PTMPC,2),'-s','Color',rgb('Crimson'),'LineWidth',1.5)
load data/population_increasing_n_SLS.mat
plot(2:nexp, toc_SLSMPC(2:nexp),'-s','Color',rgb('SeaGreen'),'LineWidth',1.5)
xlabel('Number  of uncertain entries of $A$','Interpreter','latex')
xlim([1 n_unc])
ylabel('Average CPU time [s]', 'Interpreter','latex')
grid on
legend('IM-MPC','SLS-MPC','PT-MPC', 'Interpreter','latex')
box on

figure
hold on
plot(1:n_unc, n_ineq(:,end)+n_eq(:,end),'-*','LineWidth',1.5)
plot(1:n_unc, ncostr_SLSMPC,'-s','Color',rgb('SeaGreen'),'LineWidth',1.5)
plot(1:n_unc, n_ineq_PTMPC+n_eq_PTMPC,'-s','Color',rgb('Crimson'),'LineWidth',1.5)
xlabel('Number  of uncertain entries of $[A \quad B]$','Interpreter','latex')
xlim([1 n_unc])
ylabel('Numeber of constraints', 'Interpreter','latex')
grid on
legend('IM-MPC','SLS-MPC','PT-MPC', 'Interpreter','latex')
box on
set(gca, 'YScale', 'log')
