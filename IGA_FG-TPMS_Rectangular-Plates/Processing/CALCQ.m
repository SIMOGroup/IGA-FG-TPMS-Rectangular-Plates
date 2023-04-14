function [Q,Qs]=CALCQ(k,theta,E,nu,G)
E1=E(1,:);
E2=E(2,:);
nu12=nu(1,:);
nu21=nu(2,:);
G12=G(1,:);
G23=G(2,:);
G13=G(3,:);
Q11(k)=E1(k)/(1-nu12(k)*nu21(k));               %k: number of layers
Q12(k)=E2(k)*nu12(k)/(1-nu12(k)*nu21(k));
Q22(k)=E2(k)/(1-nu12(k)*nu21(k));
Q66(k)=G12(k);
Q44(k)=G23(k);Q55(k)=G13(k);
Q=[Q11(k) Q12(k) 0;
   Q12(k) Q22(k) 0;
   0      0     Q66(k)];
Qs=[Q44(k) Q55(k)];