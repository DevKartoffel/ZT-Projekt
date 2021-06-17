function [kT, kTRnf, aRnf, aDachRnf, ARnf, bRnf, cTRnf , TRnf] = ackerSelfmade(A, b, cT, d, eigenWunsch)
    % Ausgabe: [kT, kTRnf, aRnf, aDachRnf, ARnf, bRnf, cTRnf , TRnf]
    [~, ARnf, bRnf, cTRnf , TRnf] = calcRnf4Ordnung(A,b,cT, d);

    % Die Koeffizienten des charakteristisches Polyom berechnen
    carPoly = poly(eigenWunsch);
    aDach = flip(carPoly(1,2:end));

    % Letze Zeile der RNF Matrix herausfilter und somit die Koeffizienten des
    % charakteristischen Polynoms der Systemmatrix erhalten
    a = - ARnf(size(ARnf,2),:);

    % kTrnf berechnen
    kTRnf = aDach - a;
    kT = kTRnf * inv(TRnf);
    aRnf = a;
    aDachRnf = aDach;
end