clear all;
clc;

D=0.017;
R=0.02;
m=0.27;
r=sqrt(R^2-(D/2)^2);
Ib=4.32*10^-5;
Iw=0.14025;
l=49;
b=1;
K=0.001;


a1=m+(Ib/r^2);
a2=(m*r^2+Ib)/r;
a3=m*g;
b1=Ib+Iw;
b2=2*m;
b3=b*l^2;
b4=K*l^2;
b5=(m*r^2+Ib)/r;
b6=m*g;



ut=0;
xt=0;

X0t=0;
X1t=0;

A0t=0;
A1t=0;


% Zustandsraum Elemente A,b,c,d

A=

b=

c=

d=



