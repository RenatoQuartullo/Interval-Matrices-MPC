clear 
clc
close all

% System (Fertility and survivance) parameters
n = 6;
f = [0.01 0.45 0.4 0.14 0 0];
s = [0.9 0.95 0.9 0.85 0.8];
Ahat = [f; [diag(s), zeros(n-1,1)] ];
Ahat(n,n) = 0.2;
Bhat = eye(n);
Bhat = Bhat(:,[2 3 4 5]);
[~,m] = size(Bhat);
xref = [1 1 1 1 1 1]';                              % state and input will be relative to refernce state and input
uref = pinv(Bhat)*(eye(n)-Ahat)*xref;
summax = 6.1-sum(xref);
summin = sum(-xref);
umax = [0.35 0.35 0.3 0.25]'-uref;
umin = -uref;
Nmax = 10;  
Kmax = 2*Nmax;
gamma = 1;
nrun = 50;
x0 = zeros(n,1)-xref;

% Uncertainty
Delta = zeros(n,n+m);
Delta(1,1:n) = 0.08*f;
Delta(2:n,1:n) = 0.03*Ahat(2:end,1:end);

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
% if ~CheckLMI(VertABK)
%     disp('K not satisfying LMI')
%     break
% end

% State/input constraints
F1 = [ones(1,n); -ones(1,n)];
G1 = [eye(m); -eye(m)];
f1 = [summax; -summin];
g1 = [umax; -umin];
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

% Compute terminal set (low-complexity) 
[W,~] = schur(AKhat,'real');
V = inv(W);
alpha_N = FindAlphaRPIParallelotope(ABK_vertices,V,W,F,G,K,h,'non-symmetric');
if any(isnan(alpha_N))
    disp('alpha_N not found')
    return
end
Xf.A = [V;-V];
Xf.b = [alpha_N(1:n); alpha_N(n+1:end)];  

% True system matrices
AB = myRandomSampling(matZonotope(IS),nrun);
A = AB(:,1:n,:);
B = AB(:,n+1:end,:);

%% IM-MPC
IDelta = intervalMatrix(zeros(n,m+n),Delta);
IDK = IS*[eye(n); K];
[Bc,Bsup,BKc,BKsup] = PreComputeIntervalMatrices(IDK,IDelta,K,Nmax);
for i = 1:Nmax
    [Hineq{i},gineq{i},Heq{i},geq{i},geqDyn{i},cJ{i},LB{i},UB{i}] = OPmatrices(Ahat,Bhat,i,Bc,Bsup,BKc,BKsup,X,U,Xf,K,'LP');
end
oPt.Display = 'off';

x = zeros(n,Kmax+1,nrun);
u = zeros(m,Kmax,nrun);
J = zeros(Kmax,nrun);
N = zeros(Kmax,nrun);

for i = 1:nrun

    disp(i)

    x(:,1,i) = x0;

    for k = 1:Kmax

        [sol,JN] = SolveIMMPC(x(:,k,i),Nmax,Hineq,gineq,Heq,geq,geqDyn,cJ,LB,UB,gamma,'LP',oPt);
        [J(k,i), N(k,i)] = min(JN);

        % Check initial feasibility
        if k == 1 && J(k,i) == inf
            disp('Initially infeasible')
            break
        end

        % Compute solution
        [~,v] = getSolution(sol{N(k,i)},n,m,N(k,i));
        u(:,k,i) = v(:,1); 

        % Simulate system
        x(:,k+1,i) = A(:,:,i)*x(:,k,i) + B(:,:,i)*u(:,k,i);
    end   

end

%% PLOT
% To see paper results --> load data/population_different_AB.mat
figure
subplot(4,2,[1 2 3 4])
hold on
plot(0:Kmax,sum(xref)*ones(1,Kmax+1),'--k')
for i = 1:nrun
    plot(0:Kmax,sum(x(:,:,i))+sum(xref),'b')
end
plot(0:Kmax,(summax+sum(xref))*ones(1,Kmax+1),'--g')
xlim([0 Kmax])
ylim([0 6.5])
ylabel('$y(k) = \sum_{i=1}^n x_i(k)$','Interpreter','latex')
box on
for i = 1:m
    subplot(4,2,i+4)
    hold on
    for j = 1:nrun   
        stairs(u(i,:,j)+uref(i),'b')
    end
    plot(uref(i)*ones(1,Kmax),'--k')
    plot((umax(i)+uref(i))*ones(1,Kmax),'--g')
    plot((umin(i)+uref(i))*ones(1,Kmax),'--g')
    ylim([-0.02 1.1*(umax(i)+uref(i))])
    xlim([0 Kmax-1])
    if i>2
        xlabel('$k$','Interpreter','latex')
    end
    ylabel(['$u_',num2str(i),'(k)$'],'Interpreter','latex')
    box on
end


