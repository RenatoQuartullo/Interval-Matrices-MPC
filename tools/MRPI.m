function [Xf,flag,volList] = MRPI(X,U,vertAB,K,Niter,minVolume)
% Calculates the MRPI set for the autonomous system x+ = (A+BK)x as in [1]
% [1] F. Borrelli, A. Bemporad, and M. Morari, Predictive control for 
% linear and hybrid systems. Cambridge University Press, 2017

    F = X.A;
    f = X.b;
    G = U.A;
    g = U.b;
    Xinit = Polyhedron([F; G*K], [f; g]);
    volList(1) = Xinit.volume();

    for i = 1:Niter
    
        PreS = PreSetAutonomousSystem(Xinit,vertAB,K);
        X_new = and(Xinit, PreS);
        X_new.minHRep();
    
        if X_new == Xinit
            flag = 1;
            Xf = X_new;
            return
        end

        Xinit = X_new;
        volList(i+1) = Xinit.volume();
    
        if Xinit.volume() < minVolume
            flag = 2;
            return
        end
    
        if abs(volList(i+1) - volList(i)) < 1e-2
            flag = 3;
            return
        end
    end
end