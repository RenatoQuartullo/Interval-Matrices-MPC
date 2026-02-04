function [flag, P] = CheckLMI(AK)
    
    [n,~,nv]  = size(AK);
    if ~strcmp(cvx_solver,'Mosek')
        cvx_solver Mosek
    end
    cvx_begin sdp quiet
    variable P(n,n) symmetric
    minimize (0)
        subject to
        P >= 1e-6*eye(n);
        for i = 1:nv
            P - AK(:,:,i)'*P*AK(:,:,i) == semidefinite(n);
        end
    cvx_end
    flag = strcmp(cvx_status,'Solved');
end

