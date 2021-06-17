clear all;
clc;
%% ############  Praktikum 1  ############
%% Variablen
D=0.017;
Rm=0.02;    %Metall Radius
Rg=0.01875; %Gummi Radius
mM=0.27;    %Masse Metall
mG=0.025;   %Masse Gummi

% Kugelauswahl
m=mM;
R=Rm;
%m=mG;
%R=Rg;

% Variablen der Ruhelage
x10 = 0;
x20 = 0;
x30 = 0;
x40 = 0;
x0 = [x10 x20 x30 x40];

% Eingangsgrößen
u0 = 0.2;
z0 = 0;
wsoll = 1;

ZRMMetall = calcZRM(mM,Rm, x0, u0);
% Zustandraummodel erstellen
ZRM = ZRMMetall;

%% 3.1
x0 = [5 0 0 0];
figure(2);
initial(ZRM, x0);
xlim([0 10]);

%% 3.2
figure(3);
step(ZRM);
%% 4 Analyse
clc;
n = length(ZRM.A);
display(['n = ', num2str(n)]);

%% 4.1 Vollständige Steuerbarkeit
Ss = [ZRM.B ZRM.A*ZRM.B ZRM.A^2*ZRM.B ZRM.A^3*ZRM.B];
%detSs = det(Ss);
rankSs = rank(Ss);
display(['Rang(S_s) = ', num2str(rankSs)]);
% Das Syste ist vollständig steuerbar, da det(Ss) ~= 0
% Das Syste ist vollständig steuerbar, da Rang(Ss) = n

%% 4.2 Vollständige Beobachtbarkeit
c0 = [ZRM.B ZRM.A*ZRM.B ZRM.A^2*ZRM.B ZRM.A^3*ZRM.B];
rangc0 = rank(c0);
ob = [
    ZRM.C;
    ZRM.C*ZRM.A;
    ZRM.C*ZRM.A^2;
    ZRM.C*ZRM.A^3
    ];
rankSb = rank(ob);
display(['Rang(S_b) = ', num2str(rankSb)]);
% Das Syste ist vollständig beobachtbar, da Rang(Sb) = n

%% 4.3
eigenA = eig(ZRM.A);
maxEigen = max(real(eigenA));
disp(['max{\lambda _{i}} = ', num2str(maxEigen)]);
% Das System ist nicht zustandsstabil, da der Realteil des größten Pols in 
% der rechten Halbebene liegt. Die Eigenwerte sind noch kein Beweis auf
% E/A-Stabilität.

%% 4.4
[nenn, zeahl] = ss2tf(ZRM.A, ZRM.B, ZRM.C, ZRM.D);
G = tf(nenn, zeahl)
% Erstes Hurwitz-Kriterium ist nicht erfüllt, da nicht alle a_i > 0. Das
% Sytem hat kein integriendes Verhalten, da S^k = 1 bzw k = 0.

%% 4.5
poleG = pole(G)
maxPol = max(real(poleG));
disp(['max{\lambda _{i}} = ', num2str(maxPol)]);
% Das System ist nicht E/A-stabil, weil der Realteil des größten Pols in 
% der rechten Halbebene liegt

%% 5 Normalformen
clc;

%% 5.1 Kanonische Normalform
% Spezielles Vorgehen
KnfSpez = canon(ZRM)
% In A32 und A23 sind die imaginären Teile der konugiert komplexen Pole zu
% sehen. Demensprechend handelt es sich hier nicht um eine Diagonalform.

% Normale Vorgehensweise
[Veig, eigenA] = eig(ZRM.A);
AKnfNormal = inv(Veig) * ZRM.A * Veig
bvKnfNormal = inv(Veig) * ZRM.B
cTKnfNormal = ZRM.C * Veig
% Alle Eigenwerte liegen auf der Diagonele. Die A-Matrix ist in die
% Diagonalform transformiert.


%% 5.2 RNF
qT = [0 0 0 1] * inv([ZRM.B ZRM.A*ZRM.B ZRM.A^2*ZRM.B ZRM.A^3*ZRM.B]);
TRnf = inv([
        qT;
        qT * ZRM.A;
        qT * ZRM.A^2;
        qT * ZRM.A^3
        ]);
bvRnf=inv(TRnf)*ZRM.B
cTRnf=ZRM.C*TRnf
ARnf=inv(TRnf)*ZRM.A*TRnf
[nenn, zeahl] = ss2tf(ZRM.A, ZRM.B, ZRM.C, ZRM.D);
GRNF = tf(nenn, zeahl)
% Es ist zu sehen, dass der Koeffizient a1 in der Regelungsnormalform auf 0
% gerundet wurde.

%% 5.3 BNF
% Über Tansformationsmatrix
r = inv([
    ZRM.C;
    ZRM.C*ZRM.A;
    ZRM.C*ZRM.A^2;
    ZRM.C*ZRM.A^3
    ])*[0;0;0;1];
TBnf = [r ZRM.A*r ZRM.A^2*r ZRM.A^3*r];
bvBnf = inv(TBnf)*ZRM.B
cTBnf = ZRM.C*TBnf
ABnf = inv(TBnf)*ZRM.A*TBnf

% Über Dualität der Regelungsnormalform
ABnfDual = ARnf'
bvBnfDual = cTBnf
cTBnfDual = bvBnf


%% ############  Praktikum 3  ############

% Vorbereitung: 
% Es gibt 4 unterschiedliche Güteanforderungen
% 1) Stabilitätsanforderung:
% Regelkreis muss stabil sein => y geht gegen null, wenn kein Signal
% anliegt
% 2)Forderung nach Störkompensation und Sollwertfolge:
% Regeldifferenz soll gleich null sein, sprich y ist w im eingeschwungenen
% Zustand
% 3)Dynamikanforderungen:
% Die dynamik soll eine gewünschte Form annehmen
% => Anstiegszeit T_r, Überschwingzeit T_m, Einschwungzeit T_e,
% Überschwingweite deltah, Toleranzband müssen erfüllt bzw. eingehalten werden
% 4)Robustheitsanforderungen:
% Unempfindlichkeit des Regelkreises gegenüber Änderungen der Strecke

% Lage der Pole:
% Je weiter sich der Realteil eines Pols in der linken Halbene befindet,
% desto geringer ist der Einfluss auf die Dynamik des Systems. Somit kann
% man die Dynamik des Systems über die Eigenwerte verändern. Pole in der
% der linken Halbene, deren Realteil nahe am Nullpunkt liegen werden
% dominierende Pole genannt.

%% ## 3. Entwurf Zustandrückführung

%% dominierende Pol berechnen eines PT2-Gliedes
% Anfoderungen
deltah = 0.1;
t5Proz = 1;

% Berechnung Pole der diminierenden Dynamik
temp = log(deltah)^2/(log(deltah)^2+pi^2);
ZRM.D = sqrt(temp);
omega0 = 3 / (ZRM.D * t5Proz);
% Übertragungsfunktion der PT2-Gliedes
gwsoll = tf(omega0^2, [1 2*ZRM.D*omega0 omega0^2]);

% Wunscheigenwerte
eigenWunsch = eig(gwsoll);
eigenWunsch(3,1) = -20;
eigenWunsch(4,1) = -20;


%% k^t über RNF
% Kaptiel mit RNF muss vorher durchgeführt werden

% Die Koeffizienten des charakteristisches Polyom berechnen
WunschPoly = poly(eigenWunsch);
% Die wunschkoeffizienten in einen transformierten Vektor schreiben und von
% a0 bis an-1 laufen lassen
aDach = flip(WunschPoly(1,2:end));

% Letze Zeile der RNF Matrix herausfilter und somit die Koeffizienten des
% charakteristischen Polynoms der Systemmatrix erhalten
a = - ARnf(size(ARnf,2),:);

% kTrnf berechnen
kTRnf = aDach - a;
kT = kTRnf * inv(TRnf);

%% k^t über acker
kTacker = acker(ZRM.A, ZRM.B, [eigenWunsch]);

%% ## 5 Entwurf eines Vorfilters

V = -(ZRM.C*(ZRM.A - ZRM.B * kT)^-1 * ZRM.B)^-1;
VRnf = aDach(1,1) / cTRnf(1,1);

%% ## 4. Dynamikanforderungen Testen

 out = sim('Simulink_ZT', 10);

[peakOhneFilter, IOhneFilter] = max(out.SimuYtRNFFilter.Data);
[peakMitFilter, IMitFilter] = max(out.SimuYtRNFFilterULimit.Data);

figure()
hold on;
title(['x_0 = ', num2str(x0(1))])
plot(out.SimuYtRNFFilter.Time, out.SimuYtRNFFilter.Data)
plot( out.SimuYtRNFFilterULimit.Time, out.SimuYtRNFFilterULimit.Data);
plot(out.SimuYtRNFFilter.Time(IOhneFilter), peakOhneFilter, 'x');
plot(out.SimuYtRNFFilterULimit.Time(IMitFilter), peakMitFilter, 'x');
legend('Ohne Limit','Mit Limit');
grid minor;
hold off;



%% Aufgabe 6

% Anfangsauslenkungen
% x0 = 0.05: deltah (mit Beschränkung) wird kleiner im Vgl. mit dem deltah
% ohne Anfangsauslenkung. Bleibt jedoch größer als das deltah bei 0.05
% ohne Limit. (TeO= 0.8953, TeM=0.953; ymaxO=0.2147, ymaxM=0.2223)
% x0 = 0.5: die Einschwingzeit wird weiterhin nicht erfüllt (TeO=0.839,
% TeM=0.88)
% Aber das deltah bleibt unterhalb des vorgegeben Wertes von 0.1 und erfüllt somit
% ebenso die Güteanforderung (yminO =0.1794, yminM = 0.1702 ).
% x0 = -0.5: Mit einer Limitierung kann der Regler das Ziel nicht mehr
% erreichen

%% ## 7 Störgröße

figure()
hold on;
title(['Störfunktion bei x_0 = ', num2str(x0(1))])
plot(out.SimuYtRNFStoerung.Time, out.SimuYtRNFStoerung.Data)
grid minor;
hold off;

% Die Störung kann nicht geregelt werden, da sie hinter dem Regler
% angreift. Es bleibt eine Regelabweichung um die Höhe der Strörung übrig.


%% ############  Praktikum 4  ############
% Vorbereitung
%
% 8.2 Regelentwurf durch Polvorgabe
% Vorraussetzung: vollständig steuerbare Regelstrecke
% Bei einer Ausgangsrückführung kann die Systemmatrix durch Ky verändert
% werden. Dabei können die Eigenwerte so beeinflusst werden, dass die
% Führungsübertragungsfunktion eine gewünschte dynamik erreicht und
% bestimmten Güteanforderungen erfüllt sind. Dabei ist es wichtig die Lage
% des dominierenden Pols zu bestimmen. Dies wurde schon im letzen Praktikum
% brechnet.
%
% 8.3 Sicherung der stattionären Genaugkeit
% Bei der ArF und ZRF sind meistens Vorfilter nötig, um eine station.
% Genaugkeit zu erlangen.
% ZRF: V = -(cT * (A-b*K^T)^-1 * b)^-1
% ARF: V = -(cT * (A-b*Ky*cT)^-1 * b)^-1
%
% 8.4 Störgrößenkompensation, Parameterrobustheit
% Um eine größere Robustheit zu erhalten und auch um Störgrößen komensieren
% zu können, wird der Regelkreis zuruückgeführt. Die Rückführung geschieht
% vor dem Vorfilter (e(t) = w(t) - y(t)). Dieser Regler wird dann PZ-Regler
% genannt. Weist die Regelstrecke kein integrierendes Verhalten auf, wird
% bei einem PZ-Regler jedoch eine bleibende Regelabweichung bei Störungen
% auftreten.
% Wird der Vorfilter durch einen PI-Regler ersetzt, so handelt es
% sich um einen PIZ-Regler. Der PIZ-Regler wird dann benötigt, wenn der
% Regelkreis kein integrierendes Verhalten aufweist. Dann werden auch
% Störungen ohne bleibende Regelabweichung geregelt.

%% ## 3 Güteanforderung
% 1) PIZ-Regler: Weil die Regelstrecke kein integrales Verhalten aufweißt
% und eine stationäre Genauigkeit erreicht werden soll. 
% 2)

% Anfoderungen
deltah = 0.1;
t5Proz = 1.5; % epsilon = 0.05;


% Berechnung Pole der diminierenden Dynamik
temp = log(deltah)^2/(log(deltah)^2+pi^2);
p4d = sqrt(temp);
p4omega0 = 3 / (p4d * t5Proz);

% Übertragungsfunktion der PT2-Gliedes
P4gwsoll = tf(p4omega0^2, [1 2*p4d*p4omega0 p4omega0^2]);

% Wunscheigenwerte
P4eigenWunsch = eig(P4gwsoll);
realteil = max(real(P4eigenWunsch));
P4eigenWunsch(3,1) = realteil * 5; % 5 Weil PIZ
P4eigenWunsch(4,1) = realteil * 5.1;
% +1 Pol für PI-Regler
P4eigenWunsch(5,1) = realteil * 5.2;

% Die Koeffizienten des charakteristisches Polyom berechnen
P4Poly = poly(P4eigenWunsch);
P4aDach = flip(P4Poly(1,2:end));



%% ## 4 Simulation des geschlossenes Regelkreises

%% 4.1 und 4.2 Positionierung

ZRM = ZRMMetall;
% Bestimmung der RNF Parameter
[~, P4ARnf, P4bRnf, P4cTRnf, P4TRnf] = calcRnf4Ordnung(ZRM.A, ZRM.B, ZRM.C,ZRM.D);
P4aRnf = - P4ARnf(size(P4ARnf,2),:);


% 4.1
P4VKeinSprung = 0;
% 4.2
P4VSprung = P4aRnf(1) / P4cTRnf(1);

% Zuordnung des Vorfilters
%P4V = P4VSprung;
P4V = P4VKeinSprung;

% Berechnung der unbekannten Parameter V und kT
P4Temp = [[0; 1; 0; 0; 0] [0; 0; 1; 0; 0] [0; 0; 0; 1; 0] [0; 0; 0; 0; 1] ];
P4Vk = inv([[P4cTRnf'; 0] P4Temp ])...
    * (P4aDach' - [ 0; (P4aRnf' + P4V*P4cTRnf')]);

% Zordnung V und kTRnf
P4VI = P4Vk(1);
P4kTschRnf = P4Vk(2:end)';

% Berechnung kT
P4kTsch = P4kTschRnf * P4TRnf;


%% 5.3 Gummiball

% Kalkuliere neue Strecke mit dem Gummiball
ZRMGummi = calcZRM(mG,Rg);
ZRM = ZRMGummi;