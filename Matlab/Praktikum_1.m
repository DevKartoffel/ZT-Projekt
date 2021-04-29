clear all;
clc;
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

r=sqrt(R^2-(D/2)^2);
Ib=4.32*10^-5;
Iw=0.14025;
l=0.49;
b=1;
K=0.001;
g=9.81;

a1 = (m + (Ib/r^2));
a2 = (m*r^2+Ib)*1/r;
a3 = m * g;

b1 = Ib + Iw;
b2 = 2 * m;
b3 = b * l^2;
b4 = K * l^2;
b5 = a2;
b6 = a3;

%% Variablen der Ruhelage
x10=0;
x20=0;
x30=0;
x40=0;
x0=[x10 x20 x30 x40];

% Eingangsgr��e
u0= 0;

%% A Matrix linearisiertes Modell
A11=0;
A12=1;
A13=0;
A14=0;

A21= (-a2*b6* (a1* (m*((x10^2))+b1) -b5*a2) + 2*a1*a2*m*x10* (b6*x10+u0*l ))...
    / ((a1*(m*((x10^2))+b1) -b5*a2)^2);
A22= 0;
A23= (a3*(m*((x10^2))+b1)+a2*b4) / (a1*(m*(x10^2)+b1) -b5*a2);
A24= (a2*b3) / (a1*(m*(x10^2)+b1) -b5*a2);

A31=0;
A32=0;
A33=0;
A34=1;

A41= b6* (-m* (x10^2) +b1)/(m* (x10^2) +b1)^2 ...
    - (1+a2*b5* ( 2*a1*m* (x10^2) +2*a1*b1-a2*b5) / ( a1* (m* (x10^2) +b1) -a2*b5)^2 )...
    *2*m*l*x10*u0 / (m*(x10^2)+b1)^2 - a2*b5*b6*( m*(x10^2)*(3*a1*m*(x10^2) + 2*a1*b1 - a2*b5) + b1*(-a1*b1+a2*b5) )/((m*(x10^2)+b1) * (a1*(m*(x10^2)+b1) - a2*b5))^2;
A42=0;
A43= (-b4) / (m*(x10^2)+b1) - (b5*a2*b4) / ((m*(x10^2)+b1) * (a1* (m*(x10^2)+b1) -b5*a2)) - (b5*a3) / (a1* (m*(x10^2)+b1) -b5*a2);
A44= (-b3) / (m*(x10^2)+b1) - (b5*a2*b3) / ((m*(x10^2)+b1) * (a1* (m*(x10^2)+b1) -b5*a2));

%b-Vektor linearisiertes Modell
B1=0;
B2=(-a2*l*cos(x30))/(a1*(m*((x10^2))+b1)-b5*a2);
B3=0;
B4=(1+(b5*a2)/(a1*(m*((x10^2))+b1)-b5*a2)) * (l*cos(x30))/(m*((x10^2))+b1);


% Zustandsraum Elemente A,b,cT,d
A=[ A11 A12 A13 A14;
    A21 A22 A23 A24;
    A31 A32 A33 A34;
    A41 A42 A43 A44];
bv=[B1;
    B2;
    B3;
    B4];
cT=[1 0 0 0];
d=0;

%% Zustandraummodel erstellen
ZRM = ss(A,bv,cT,d);

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
n = length(A);
display(['n = ', num2str(n)]);

%% 4.1 Vollst�ndige Steuerbarkeit
Ss = [bv A*bv A^2*bv A^3*bv];
%detSs = det(Ss);
rankSs = rank(Ss);
display(['Rang(S_s) = ', num2str(rankSs)]);
% Das Syste ist vollst�ndig steuerbar, da det(Ss) ~= 0
% Das Syste ist vollst�ndig steuerbar, da Rang(Ss) = n

%% 4.2 Vollst�ndige Beobachtbarkeit
c0 = [bv A*bv A^2*bv A^3*bv];
rangc0 = rank(c0);
ob = [
    cT;
    cT*A;
    cT*A^2;
    cT*A^3
    ];
rankSb = rank(ob);
display(['Rang(S_b) = ', num2str(rankSb)]);
% Das Syste ist vollst�ndig beobachtbar, da Rang(Sb) = n

%% 4.3
eigenA = eig(A);
maxEigen = max(real(eigenA));
disp(['max{\lambda _{i}} = ', num2str(maxEigen)]);
% Das System ist nicht zustandsstabil, da der Realteil des gr��ten Pols in 
% der rechten Halbebene liegt. Die Eigenwerte sind noch kein Beweis auf
% E/A-Stabilit�t.

%% 4.4
[nenn, zeahl] = ss2tf(A, bv, cT, d);
G = tf(nenn, zeahl)
% Erstes Hurwitz-Kriterium ist nicht erf�llt, da nicht alle a_i > 0. Das
% Sytem hat kein integriendes Verhalten, da S^k = 1 bzw k = 0.

%% 4.5
poleG = pole(G)
maxPol = max(real(poleG));
disp(['max{\lambda _{i}} = ', num2str(maxPol)]);
% Das System ist nicht E/A-stabil, weil der Realteil des gr��ten Pols in 
% der rechten Halbebene liegt

%% 5 Normalformen
clc;

%% 5.1 Kanonische Normalform
% Spezielles Vorgehen
KnfSpez = canon(ZRM)
% In A32 und A23 sind die imagin�ren Teile der konugiert komplexen Pole zu
% sehen. Demensprechend handelt es sich hier nicht um eine Diagonalform.

% Normale Vorgehensweise
[V, eigenA] = eig(A);
AKnfNormal = inv(V) * A * V
bvKnfNormal = inv(V) * bv
cTKnfNormal = cT * V
% Alle Eigenwerte liegen auf der Diagonele. Die A-Matrix ist in die
% Diagonalform transformiert.


%% 5.2 RNF
qT=cT * inv([bv A*bv A^2*bv A^3*bv]);
T=inv([
    qT;
    qT * A;
    qT * A^2;
    qT * A^3
    ]);
bvRnf=inv(T)*bv
cTRnf=cT*T
ARnf=inv(T)*A*T
[nenn, zeahl] = ss2tf(A, bv, cT, d);
G = tf(nenn, zeahl)
% Es ist zu sehen, dass der Koeffizient a1 in der Regelungsnormalform auf 0
% gerundet wurde.


%% 5.3 BNF
% �ber Tansformationsmatrix
r=inv([
    cT;
    cT*A;
    cT*A^2;
    cT*A^3
    ])*[0;0;0;1];
T=[r A*r A^2*r A^3*r];
bvBnf=inv(T)*bv
cTBnf=cT*T
ABnf=inv(T)*A*T

% �ber Dualit�t der Regelungsnormalform
ABnfDual = ARnf'
bvBnfDual = cTBnf
cTBnfDual = bvBnf
