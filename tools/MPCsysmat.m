function [F,G] = MPCsysmat(A,B,N)
% Calculation of concatenated system matrices F and G, given A, B and an
% horizon length N

    [n,m] = size(B);

    F = [];
    Gt = [];
    G = zeros(n*N,m*N);

    for j = 1:N
        F = [F;A^j];
        Gt = [A^(j-1)*B Gt];
        G(1+n*(j-1):j*n,1:m*j) = Gt;
    end
end