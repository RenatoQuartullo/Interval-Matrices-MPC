function [K,S,Y] = KStableVertices(VertAM,VertBM)
    
    [n,m,nv] = size(VertBM);
    
    if ~strcmp(cvx_solver,'Mosek')
        cvx_solver Mosek
    end
    cvx_begin sdp quiet
        variable S(n,n) symmetric
        variable Y(m,n) 
    %     variable H(n,n) symmetric
    %     variable gamma_opt(1,1)
        minimize(trace(S)+ 0.0*norm(Y,'fro') )
        subject to
            S >= 1e-6*eye(n);
    %         H == semidefinite(n);
            for j = 1:nv
    %             [[S, (VertAM(:,:,j)*S + VertBM(:,:,j)*Y)'; VertAM(:,:,j)*S + VertBM(:,:,j)*Y, S], [S*sqrt(Q), Y'*sqrt(R); zeros(n,m+n)]; zeros(n+m, 2*n) , gamma_opt*eye(n+m)] >= 1e-6*eye(2*n + n + m);
                [S, (VertAM(:,:,j)*S + VertBM(:,:,j)*Y)'; VertAM(:,:,j)*S + VertBM(:,:,j)*Y, S] == semidefinite(2*n);
            end
    cvx_end
    K = Y/S;
    
end

