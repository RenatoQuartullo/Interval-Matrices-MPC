function [c, alpha, s0] = getSolution_PTMPC(solCannon,m,n,N,nv)
% Returns the solution from quadorog of PT-MPC
    
    s0 = solCannon(end-n+1:end);
    c = reshape(solCannon(end-n-m*N+1:end-n),[m,N]);
    alpha = reshape(solCannon(1:(N+1)*nv),[nv,N+1]);

end

