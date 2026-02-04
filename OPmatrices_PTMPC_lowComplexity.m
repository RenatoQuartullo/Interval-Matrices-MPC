function [Hineq,gineq,Heq,geq,Qqp] = OPmatrices_PTMPC_lowComplexity(N, Phi_tilde_p, Phi_tilde_m, B_tilde, FGK_p, FGK_m, G, h, V, alpha_N, Wp)

% solCannon = [alpha_l(:); alpha_u(:); c(:); s0]. NB: gineq will be
% [gineq; V*xk; -V*xk].

    [n,m,nv] = size(B_tilde);

    % inequalities
    HalphaL = zeros(n*nv*N,2*n*(N+1)+m*N+n);
    row_counter = 1;
    for j = 1:nv
        for i = 0:N-1    
            Hij = [zeros(n,n*i), -Phi_tilde_p(:,:,j), eye(n), zeros(n,(N+1)*n-(i+2)*n), zeros(n,n*i), Phi_tilde_m(:,:,j), zeros(n,(N+1)*n-(i+1)*n), zeros(n,m*i), -B_tilde(:,:,j), zeros(n,N*m-(i+1)*m), zeros(n)];
            row_idx = (row_counter-1)*n + 1 : row_counter*n;
            HalphaL(row_idx, :) = Hij;
            row_counter = row_counter + 1;
        end
    end
    HalphaU = zeros(n*nv*N,2*n*(N+1)+m*N+n);
    row_counter = 1;
    for j = 1:nv
        for i = 0:N-1
            Hij = -[zeros(n,n*i), Phi_tilde_m(:,:,j), zeros(n,(N+1)*n-(i+1)*n), zeros(n,n*i), -Phi_tilde_p(:,:,j), eye(n), zeros(n,(N+1)*n-(i+2)*n), zeros(n,m*i), -B_tilde(:,:,j), zeros(n,N*m-(i+1)*m), zeros(n)];
            row_idx = (row_counter-1)*n + 1 : row_counter*n;
            HalphaU(row_idx, :) = Hij;
            row_counter = row_counter + 1;
        end
    end
    Halpha = [HalphaL; HalphaU];
    galpha = zeros(size(Halpha,1),1);
    nY = size(Halpha,2);

    ncsi = size(FGK_m,1);
    HStateInput = [-kron(eye(N),FGK_m), zeros(ncsi*N,n), kron(eye(N),FGK_p), zeros(ncsi*N,n), kron(eye(N),G), zeros(ncsi*N,n)];
    gStateInput = repmat(h,[N,1]);

    Halpha0 = zeros(2*n,nY);
    Halpha0(1:n,1:n) = eye(n);
    Halpha0(n+1:end,(N+1)*n+1:(N+1)*n+n) = -eye(n); % galpha0 will be [V*xk; -V*xk]

    Hs0 = Halpha0;
    Hs0(1:n,2*n*(N+1)+m*N+1:2*n*(N+1)+m*N+n) = -V;
    Hs0(1+n:end,2*n*(N+1)+m*N+1:2*n*(N+1)+m*N+n) = V;
    gs0 = zeros(2*n,1);

    % equalities
    Halphaend = zeros(2*n,nY);
    Halphaend(1:n,n*N+1:n*N+n) = eye(n);
    Halphaend(1+n:end,n*(2*N+1)+1:n*(2*N+1)+n) = eye(n);
    galphaend = [-alpha_N(n+1:end); alpha_N(1:n)];

    % Total matrices
    Hineq = sparse([Halpha; HStateInput; Hs0; Halpha0]);
    gineq = sparse([galpha; gStateInput; gs0]);
    Heq = sparse(Halphaend);
    geq = sparse(galphaend);
    
    % Cost
    S = zeros(n + N*m, n + N*m);
    S(1:n, N*m+(1:n)) = eye(n);
    S(n+(1:N*m), 1:N*m) = eye(N*m);
    Qqp = S'*Wp*S;
    Qqp = sparse(blkdiag(zeros(2*n*(N+1)),Qqp));
    
end
