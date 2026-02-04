function ObjOut = PowerMatIntTimesObj(IM,j,Obj)
% j Successive applicatio of the T-operator, referred as T^j_IM(Obj)
    
    ObjOut = Obj;
    
    if isa(Obj,'matZonotope')
        for i = 1:j
            ObjOut = MatIntTimesMatZon(IM,ObjOut);
        end
    else
        for i = 1:j
            ObjOut = IM*ObjOut;
        end
    end
    
end

