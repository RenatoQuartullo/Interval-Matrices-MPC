function [K,P] = KStableVertices(VertAM,VertBM)
    
    [n,m,nv] = size(VertBM);

    S = sdpvar(n,n);
    Y = sdpvar(m,n);
    rho = 1;
    constr = S >= 1e-6*eye(n);
    for j = 1:nv
        mat = [rho*S, (VertAM(:,:,j)*S + VertBM(:,:,j)*Y)'; ...
        VertAM(:,:,j)*S + VertBM(:,:,j)*Y, S];
        constr = [constr, mat >= 0];
    end
    ops = sdpsettings('verbose', 0,'solver', 'mosek');
    optimize(constr, trace(S), ops);
    Svalue = value(S);
    Yvalue = value(Y);
    K = Yvalue/Svalue;
    P = inv(Svalue);

end

