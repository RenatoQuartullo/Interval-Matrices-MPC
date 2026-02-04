function MZout = MatIntTimesMatZon(MI,MZ)
% MZout = T_{MI}(MZ), where the operator T is the one in the referred paper
    
    [infMI, supMI] = MatIntInfSup(MI);
    L = (supMI + infMI)/2;          % interval center
    S = (supMI - infMI)/2;          % symmetric part
    [n,m] = size(L);
    
    C = MZ.center;
    G = MZ.G;
    [~,p,e] = size(G);
    
    MZoutCenter = L*C;
    
    MZoutGen = zeros(n,p,e+n*p);
    for i = 1:e
        MZoutGen(:,:,i) = L*G(:,:,i);
    end
    
    F = S*(abs(C) + sum(abs(G),3));
    
    for i = 1:n
        for j = 1:p
            MZoutGen(i,j,e+(i-1)*p+j) = F(i,j);
        end
    end
    
    MZout = matZonotope(MZoutCenter, MZoutGen);
    
end

