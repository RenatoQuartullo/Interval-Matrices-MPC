function [Hineq,gineq,Heq,geq,geqDyn,cJ,LB,UB] = OPmatrices(A,B,N,Bc,Bsup,BKc,BKsup,X,U,Xf,K,OP)
% A and B are the (nominal) system matrices, N is the prediction horizon, 
% Bc and Bsup are the centers and supremum of matrix intervals B_k, 
% respectively (the same for BKc and BKsup for the input constrainte) 
% and X, U are the (polytopic) state/input constraints (possibly
% time-vaying). The outputs are s.t. the constraints of OP will be 
% Hineq*Y <= gineq, Heq*Y == [geqDyn*x(k); geq]. 
% The solution given by linporog will be Y = [z;v;|z|;|v|; t; |t|]. 
% t(j) = |v(j)-Kz(j)|  for the cost J = N+\gamma\sum_j|v(j)-Kz(j)|. T
% the terminal set Xf is a polytope and Xf = {x| Xf.A x <= Xf.b}

    [n,m] = size(B);
    
    % dynamic constraint quantities
    [Fm, Gm] = MPCsysmat(A,B,N);
    HeqDyn = [eye(n), zeros(n, n*N+m*N), zeros(n,m*N+n*(N+1)+2*m*N); zeros(n*N,n), eye(n*N), -Gm, zeros(n*N,m*N+n*(N+1)+2*m*N)];      
    geqDyn = [eye(n); Fm];                      % Initial and dynamic constraints will be HeqDyn*X == geqDyb*x(k)

    % constraint for absolute values
    HinAbs = [eye(n*(N+1)), zeros(n*(N+1),m*N), -eye(n*(N+1)), zeros(n*(N+1),m*N+2*m*N);...
             -eye(n*(N+1)), zeros(n*(N+1),m*N), -eye(n*(N+1)), zeros(n*(N+1),m*N+2*m*N);...
              zeros(m*N,n*(N+1)), eye(m*N), zeros(m*N,n*(N+1)), -eye(m*N), zeros(m*N,2*m*N);...
              zeros(m*N,n*(N+1)), -eye(m*N), zeros(m*N,n*(N+1)), -eye(m*N), zeros(m*N,2*m*N);...
              zeros(m*N,2*n*(N+1)+2*m*N), eye(m*N), -eye(m*N);...
              zeros(m*N,2*n*(N+1)+2*m*N), -eye(m*N), -eye(m*N)];
    ginAbs = zeros(2*(n*(N+1)+m*N)+2*m*N,1);

    % state and input inequality constraint matrix for each j
    HinState = [];
    ginState = [];
    HinInput = [];
    ginInput = [];
    for j = 1:N
        HinState = [HinState; BuildStateConstr(X{j}.A, j-1, N, Bc, Bsup)];
        ginState = [ginState; X{j}.b];
        HinInput = [HinInput; BuildInputConstr(U{j}.A, j-1, N, BKc, BKsup)];
        ginInput = [ginInput; U{j}.b];
    end
    HinState = [HinState, zeros(size(HinState,1),2*m*N)];
    HinInput = [HinInput, zeros(size(HinInput,1),2*m*N)];
    
    % terminal constraint
    HinTerminal = BuildStateConstr(Xf.A, N, N, Bc, Bsup);
    HinTerminal = [HinTerminal, zeros(size(HinTerminal,1),2*m*N)];
    ginTerminal = Xf.b;
    
    Kblk = [kron(eye(N),K), zeros(m*N,n)];
    HeqCost = [-Kblk, eye(m*N), zeros(m*N,n*(N+1)+m*N), -eye(m*N), zeros(m*N)];
    geqCost = zeros(m*N,1);
    if strcmp(OP,'LP')
        % sum-of-norms cost
        cJ = sparse([zeros(2*n*(N+1)+3*m*N,1); ones(m*N,1)]);
    else
        cJ = sparse(blkdiag(1e-7*eye(2*n*(N+1)+3*m*N),kron(eye(N),eye(m))));
    end
    
    % all quantites for quadprog
    Heq = sparse([HeqDyn; HeqCost]);
    geq = sparse(geqCost);
    Hineq = sparse([HinAbs; HinState; HinInput; HinTerminal]);
    gineq = sparse([ginAbs; ginState; ginInput; ginTerminal]);
     
    % LB and UB for Y
    LB = sparse([-inf(n*(N+1)+m*N,1); zeros(n*(N+1)+m*N,1); -inf(m*N,1); zeros(m*N,1)]);
    UB = sparse([inf(n*(N+1)+m*N,1); inf(n*(N+1)+m*N,1); inf(m*N,1); inf(m*N,1)]);
    
end

