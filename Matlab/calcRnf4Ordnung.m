function [GRnf, ARnf, bRnf, cTRnf, TRnf] = calcRnf4Ordnung(A,bv,cT,d)
    % Ausgabe: [GRnf, ARnf, bRnf, cTRnf, TRnf] 
    
    qT = [0 0 0 1] * inv([bv A*bv A^2*bv A^3*bv]);
    TRnf = inv([
            qT;
            qT * A;
            qT * A^2;
            qT * A^3
            ]);
    bRnf=inv(TRnf)*bv;
    cTRnf=cT*TRnf;
    ARnf=inv(TRnf)*A*TRnf;
    [nenn, zeahl] = ss2tf(A, bv, cT, d);
    GRnf = tf(nenn, zeahl);

end