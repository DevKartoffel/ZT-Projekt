function [kT, kTRnf, aRnf, aDachRnf, ARnf, bRnf, cTRnf , TRnf] = ackerSelfmade(A, b, cT, d, eigenWunsch)
    % Ausgabe: [kT, kTRnf, aRnf, aDachRnf, ARnf, bRnf, cTRnf , TRnf]
    [~, ARnf, bRnf, cTRnf , TRnf] = calcRnf4Ordnung(A,b,cT, d);

    % Eigenwerte auf die Diagonale einer Matrix bringen
    lambdaDach = zeros(size(eigenWunsch,1),size(eigenWunsch,1));
    for n=1:size(eigenWunsch,1)
        lambdaDach(n,n) = eigenWunsch(n);
    end
    % Die Koeffizienten des charakteristisches Polyom berechnen
    carPoly = charpoly(lambdaDach);
    % Die wunschkoeffizienten in einen transformierten Vektor schreiben und von
    % a0 bis an-1 laufen lassen
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