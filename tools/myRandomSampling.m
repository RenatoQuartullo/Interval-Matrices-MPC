function  Samples = myRandomSampling(MZ,N)
    
    MZgen = MZ.G;
    g = size(MZgen,3);
    betas = 1-2*rand(g,N);
    Samples = repmat(MZ.center,[1 1 N]);
    
    for i = 1:N
        for j = 1:g
            Samples(:,:,i) = Samples(:,:,i) + betas(j,i)*MZgen(:,:,j);
        end
    end     
    
end

