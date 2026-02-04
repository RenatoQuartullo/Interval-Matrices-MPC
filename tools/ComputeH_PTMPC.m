function H = ComputeH_PTMPC(V,ABK)

    oPt.Display = 'off';
    nvert = size(ABK,3);
    nv = size(V,1);
    c = ones(nv,1);
    H = zeros(nv,nv,nvert);
    for j = 1:nvert
        for i = 1:nv
            hiOpt = linprog(c,[],[],V',(V(i,:)*ABK(:,:,j))',zeros(nv,1),[],oPt);
            H(i,:,j) = hiOpt';
        end
    end
    
end