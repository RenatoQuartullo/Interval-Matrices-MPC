function [Hineq,gineq,galpha0,Heq,geq,Qqp] = OPmatrices_PTMPC(N, H, Hc, V, B, G, Wp, alpha_N)
% Returns the quadpprog-style quantites needed for solving PT-MPC problem.
% All the quantities involved are in [1, Sec. 5.5]
% [1] B. Kouvaritakis and M. Cannon, Model predictive control. Springer, 2016.
% The solution vector is [alpha; c; s0]

    [n,m,nvert] = size(B);
    [nc,nv] = size(Hc);
    lenSol = nv*(N+1) + m*N + n;

    % state/input inequality constraints
    HStateInput = [kron(eye(N),Hc), zeros(nc*N,nv), kron(eye(N),G), zeros(nc*N,n)];
    gStateInput = ones(nc*N,1);

    % initial alpha -> NB: Halpha0 <= galpha0*xk
    Halpha0 = [-eye(nv), zeros(nv,lenSol-nv)];
    galpha0 = -V;

    % initial alpha/s0
    Hs0 = [-eye(nv), zeros(nv, lenSol-n-nv), V];
    gs0 = zeros(nv,1);

    % dynamic constraint
    Hdyn = [];
    for j = 1:nvert
        for i = 1:N
            Hdyni = zeros(nv,lenSol);
            Hdyni(:,nv*(i-1)+1:nv*i) = H(:,:,j);
            Hdyni(:,nv*i+1:nv*i+nv) = -eye(nv);
            Hdyni(:,nv*(N+1)+(i-1)*m+1:nv*(N+1)+i*m) = V*B(:,:,j);
            Hdyn = [Hdyn; Hdyni];
        end
    end
    gdyn = zeros(nv*N*nvert,1);

    % Terminal constraint 
    Heq = sparse(nv,lenSol);
    Heq(:,nv*N+1:nv*(N+1)) = eye(nv);
    geq = sparse(alpha_N);

    % Total matrices
    Hineq = sparse([Halpha0; Hdyn; HStateInput; Hs0]);
    gineq = sparse([gdyn; gStateInput; gs0]);
    
    % Cost
    if isempty(Wp)
        Qqp = [];
    else
        S = zeros(n + N*m, n + N*m);
        S(1:n, N*m+(1:n)) = eye(n);
        S(n+(1:N*m), 1:N*m) = eye(N*m);
        Qqp = S'*Wp*S;
        Qqp = sparse(blkdiag(zeros(nv*(N+1)),Qqp));
    end
    
end


