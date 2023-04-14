function [Qs_Bars]=CALCQsBAR(k,theta,Qs)
Q44_bar(k)=Qs(1,1)*cos(theta(k))^2+Qs(1,2)*sin(theta(k))^2;
Q45_bar(k)=(Qs(1,2)-Qs(1,1))*cos(theta(k))*sin(theta(k));
Q55_bar(k)=Qs(1,1)*sin(theta(k))^2+Qs(1,2)*cos(theta(k))^2;
Qs_Bars=[Q55_bar(k) Q45_bar(k);
        Q45_bar(k) Q44_bar(k)];
