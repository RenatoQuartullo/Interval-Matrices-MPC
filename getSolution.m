function [z,v,zabs,vabs,t,tabs] = getSolution(Y,n,m,N)
% Given the solution Y from linprog/quadprog, the state dimension n, the input 
% dimension m and the horizon N, returns the nominal state, input, their 
% absolute values and the quantity t = v-Kz

    z = reshape(Y(1:n*(N+1)),n,N+1);
    v = reshape(Y(1+n*(N+1):n*(N+1)+m*N),m,N);
    zabs = reshape(Y(n*(N+1)+m*N+1:n*(N+1)+m*N+n*(N+1)),n,N+1);
    vabs = reshape(Y(1+n*(N+1)+m*N+n*(N+1):n*(N+1)+m*N+n*(N+1)+m*N),m,N);
    t = Y(end-2*m*N+1:end-m*N);
    tabs = Y(end-m*N+1:end);
end

