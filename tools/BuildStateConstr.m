function Ht = BuildStateConstr(H, j, N, Bc, Bsup)
% Returns the matices Ht for the state constraint Ht*Y<=g_u, taking
% into account the shape of the solution vector, the horizon length. Bc
% and Bsup are for the constraint tightening

    [c,n] = size(H);
    [~,nm,~] = size(Bc);
    m = nm-n;
    
    Ht = zeros(c,2*(m*N+n*(N+1)));
    
    for i = 0:j-1
        Ht(:,i*nm+1:(i+1)*nm) = H*Bc(:,:,j-i-1+1);
    end
    Ht(:,j*nm+1:j*nm+n) = H;
    
    for i = 0:j-1
        Ht(:,m*N+(N+1)*n+i*nm+1:m*N+(N+1)*n+(i+1)*nm) = abs(H)*Bsup(:,:,j-i-1+1);
    end
    
    % Permutate Ht
    P = build_permutation_matrix(N,n,m);
    Ht = Ht*blkdiag(P,P);
end

