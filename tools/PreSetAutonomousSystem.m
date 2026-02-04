function PreS = PreSetAutonomousSystem(S,vertAB,K)
    [n,m,nv] = size(vertAB);
    F = S.A; f = S.b;  
    Fbar = []; fbar = [];
    for i = 1:nv
        A = vertAB(:,1:n,i); 
        B = vertAB(:,n+1:end,i);
        constr = F*(A + B*K);
        Fbar = [Fbar; constr]; 
        fbar = [fbar; f];
    end
    PreS = Polyhedron(Fbar, fbar);
end