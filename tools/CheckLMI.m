function flag = CheckLMI(AK_vertices)
    
    [n,~,nv]  = size(AK_vertices);
    P = sdpvar(n,n);
    constr = P >= 1e-6*eye(n);
    for i = 1:nv
        mat = P - AK_vertices(:,:,i)'*P*AK_vertices(:,:,i);
        constr = [constr, mat >= 0];
    end
    ops = sdpsettings('verbose', 0,'solver', 'mosek');
    sol = optimize(constr, [], ops);
    if sol.problem == 0
        flag = 1;
    else
        flag = 0;
    end
end
