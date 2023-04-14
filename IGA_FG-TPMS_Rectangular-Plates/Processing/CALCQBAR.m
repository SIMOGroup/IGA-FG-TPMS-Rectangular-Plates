function [Q_Bar]=CALCQBAR(k,theta,Q)
Q11_bar(k)=Q(1,1)*cos(theta(k))^4+2*(Q(1,2)+2*Q(3,3))*sin(theta(k))^2*cos(theta(k))^2+Q(2,2)*sin(theta(k))^4;
Q12_bar(k)=(Q(1,1)+Q(2,2)-4*Q(3,3))*sin(theta(k))^2*cos(theta(k))^2+Q(1,2)*(sin(theta(k))^4+cos(theta(k))^4);
Q22_bar(k)=Q(1,1)*sin(theta(k))^4+2*(Q(1,2)+2*Q(3,3))*sin(theta(k))^2*cos(theta(k))^2+Q(2,2)*cos(theta(k))^4;
Q16_bar(k)=(Q(1,1)-Q(1,2)-2*Q(3,3))*sin(theta(k))*cos(theta(k))^3+(Q(1,2)-Q(2,2)+2*Q(3,3))*sin(theta(k))^3*cos(theta(k));
Q26_bar(k)=(Q(1,1)-Q(1,2)-2*Q(3,3))*sin(theta(k))^3*cos(theta(k))+(Q(1,2)-Q(2,2)+2*Q(3,3))*sin(theta(k))*cos(theta(k))^3;
Q66_bar(k)=(Q(1,1)+Q(2,2)-2*Q(1,2)-2*Q(3,3))*sin(theta(k))^2*cos(theta(k))^2+Q(3,3)*(sin(theta(k))^4+cos(theta(k))^4);
Q_Bar=[Q11_bar(k) Q12_bar(k) Q16_bar(k);
       Q12_bar(k) Q22_bar(k) Q26_bar(k);
       Q16_bar(k) Q26_bar(k) Q66_bar(k)];

   