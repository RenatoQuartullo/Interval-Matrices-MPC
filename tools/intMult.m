function IMult = intMult(I1,I2)

    [C1,D1] = MatIntSymmetric(I1);
    [C2,D2] = MatIntSymmetric(I2);

    IMult = intervalMatrix(C1 * C2,  abs(C1)*D2 + D1*abs(C2) + D1*D2);
    
end