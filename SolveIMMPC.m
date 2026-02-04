function [sol,J] = SolveIMMPC(xk,Nmax,Hineq,gineq,Heq,geq,geqDyn,cJ,LB,UB,gamma,OP,oPt)
% Solution of variable-horizon scheme IM-MPC by enumeration for N=1...Nmax
% Hineq, gineq, etc are structs for all N.

    J = nan(Nmax,1);
    sol = cell(Nmax,1);
    
    for N = 1:Nmax
        geqk = sparse([geqDyn{N}*xk; geq{N}]); 
        if strcmp(OP,'LP')
            [sol{N},fv,flag] = linprog(cJ{N},Hineq{N},gineq{N},Heq{N},geqk,LB{N},UB{N},oPt);
        elseif strcmp(OP,'QP')
            [sol{N},fv,flag] = quadprog(cJ{N},zeros(size(cJ{N},2),1),Hineq{N},gineq{N},Heq{N},geqk,LB{N},UB{N},[],oPt);
        else
            [sol{N},fv,flag] = linprog([],Hineq{N},gineq{N},Heq{N},geqk,LB{N},UB{N},oPt);
        end
        if flag ~= 1
            fv = inf;
        end
        J(N) = gamma*N + fv;
    end
end
