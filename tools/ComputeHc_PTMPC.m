function Hc = ComputeHc_PTMPC(V,F,G,K)

    oPt.Display = 'off';
    FGK = F + G*K;
    nc = size(FGK,1);
    nv = size(V,1);
    c = ones(nv,1);
    
    Hc = zeros(nc,nv);
    for i = 1:nc
        hiOpt = linprog(c,[],[],V',FGK(i,:)',zeros(nv,1),[],oPt);
        Hc(i,:) = hiOpt';
    end

end