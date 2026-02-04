function Ht = BuildInputConstr(H, j, N, BKc, BKsup)
    % Returns the matices Ht for the input constraint Ht*Y<=g_u, taking
    % into account the shape of the solution vector, the horizon length. BKc
    % and BKsup are for the constraint tightening
    
    [c,m] = size(H);
    [~,nm,~] = size(BKc);
    n = nm-m;
    
    Ht = zeros(c,2*(m*N+n*(N+1)));
    
    for i = 0:j-1
        Ht(:,i*nm+1:(i+1)*nm) = H*BKc(:,:,j-i-1+1);
    end
    Ht(:,j*nm+n+1:j*nm+nm) = H;
    
    for i = 0:j-1
        Ht(:,m*N+(N+1)*n+i*nm+1:m*N+(N+1)*n+(i+1)*nm) = abs(H)*BKsup(:,:,j-i-1+1);
    end
    
    % Permutate Ht
    P = build_permutation_matrix(N,n,m);
    Ht = Ht*blkdiag(P,P);
end

