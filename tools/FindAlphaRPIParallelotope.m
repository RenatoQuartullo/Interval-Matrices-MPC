function alpha_N = FindAlphaRPIParallelotope(VertABK,V,W,F,G,K,h,symm)
% Computes the vector alpha_N for the low-complexity (parallelotopic) RPI
% set {x| -alpha_N<=|Vx|<=alpha_N} (if symm == 'symmetric') or 
% {x| -alpha_N^L<=|Vx|<=alpha_N^U} (otheriwise) w.r.t. the parametric
% uncertainty affecting the A+BK matrix and satysfing the constraints.

    [n,~,nv] = size(VertABK);
    
    % Perron-Frobenius condition
    phibar = zeros(n,n);
    phibarik = zeros(1,nv);
    for i = 1:n
        for k = 1:n
            for j = 1:nv
                phibarik(j) = V(i,:)*VertABK(:,:,j)*W(:,k);
            end
        phibar(i,k) = max(phibarik);
        end
    end
    lambdaPF = max(abs(eig(phibar))); 
    if lambdaPF > 1
        disp('lambda PF > 1')
        alpha_N = nan(n,1);
%         return
    end
    
    % Optimization problem with symmetric parallelotope
    options = optimoptions( 'fmincon','Algorithm','interior-point','SpecifyObjectiveGradient', true, 'Display','off');
    if strcmp(symm,'symmetric')
        H = abs((F+G*K)*W);
        for j = 1:nv
            H = [H; abs(V*VertABK(:,:,j)*W)-eye(n)];
        end
        g = [h; zeros(n*nv,1)];
        [alpha_N,~,flag] = fmincon(@(x)deal(-sum(log(x)), -1./x), ones(n,1), H, g, [], [], 1e-5*ones(n,1), [], [], options);
        alpha_N = [alpha_N; alpha_N];
    else
    % Optimization problem with non-symmetric parallelotope
        H = [max((F+G*K)*W,0), max(-(F+G*K)*W,0)];
        for j = 1:nv
            VphiW = V*VertABK(:,:,j)*W;
            M = [max(VphiW,0), max(-VphiW,0); max(-VphiW,0), max(VphiW,0)];
            H = [H; M-eye(2*n)];
        end
        g = [h; zeros(2*n*nv,1)];
        [alpha_N,~,flag] = fmincon(@(x)deal(-sum(log(x)), -1./x), ones(2*n,1), H, g, [], [], 1e-5*ones(2*n,1), [], [], options);    
    end
    
    if flag~=1
        alpha_N = nan(2*n,1);
    end

end
