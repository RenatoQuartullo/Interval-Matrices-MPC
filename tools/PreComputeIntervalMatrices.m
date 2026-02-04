function [Bc,Bsup,BKc,BKsup] = PreComputeIntervalMatrices(MDK,MDelta,K,N)
% Calculates the matrices for the constraint sets tightening

    [n,nm] = size(MDelta);
    m = nm-n;
    MDeltaMZ = matZonotope(MDelta);
    MDKMdelta{1} = MDeltaMZ;                            % MDKMdelta{j+1} = MDKint^j*Mdelta = M_T(j) (report)
    KMDKMdelta{1} = K*MDeltaMZ;
    BT{1} = intervalMatrix(MDKMdelta{1});
    BKT{1} = intervalMatrix(KMDKMdelta{1});
    for j = 1:N
        MDKMdelta{j+1} = PowerMatIntTimesObj(MDK,j,MDeltaMZ);
        KMDKMdelta{j+1} = K*MDKMdelta{j+1};
        BT{j+1} = intervalMatrix(MDKMdelta{j+1}); 
        BKT{j+1} = intervalMatrix(KMDKMdelta{j+1});
    end
    Bc = zeros(n,n+m,N+1);
    Bsup = zeros(n,n+m,N+1);
    BKc = zeros(m,n+m,N+1);
    BKsup = zeros(m,n+m,N+1);
    for j=1:N+1
        [Bi,Bs] = MatIntInfSup(BT{j});
        Bc(:,:,j) = (Bi+Bs)/2;
        Bsup(:,:,j) = (-Bi+Bs)/2;
        [Bi,Bs] = MatIntInfSup(BKT{j});
        BKc(:,:,j) = (Bi+Bs)/2;
        BKsup(:,:,j) = (-Bi+Bs)/2;
    end
    
end

