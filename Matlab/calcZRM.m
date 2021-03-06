function [ZRM] = calcZRM(m, R, x0, u0)

    x10 = x0(1);
    x30 = x0(3);

    D=0.017;
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

    % A Matrix linearisiertes Modell
    A11 = 0;
    A12 = 1;
    A13 = 0;
    A14 = 0;

    A21 = (-a2*b6* (a1* (m*((x10^2))+b1) -b5*a2) + 2*a1*a2*m*x10* (b6*x10+u0*l ))...
        / ((a1*(m*((x10^2))+b1) -b5*a2)^2);
    A22 = 0;
    A23 = (a3*(m*((x10^2))+b1)+a2*b4) / (a1*(m*(x10^2)+b1) -b5*a2);
    A24 = (a2*b3) / (a1*(m*(x10^2)+b1) -b5*a2);

    A31 = 0;
    A32 = 0;
    A33 = 0;
    A34 = 1;

    A41 = b6* (-m* (x10^2) +b1)/(m* (x10^2) +b1)^2 ...
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

    % Zustandraummodel erstellen
    ZRM = ss(A,bv,cT,d);

end