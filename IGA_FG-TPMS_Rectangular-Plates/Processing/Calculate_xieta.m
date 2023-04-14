function [coord]=Calculate_xieta (c_ele,c_point,ien,gcoord,u_knot,v_knot,b_net)

coord = zeros(1,2);
ksi   = 0;
eta   = 0;
iter  = 30;

sctr = ien(c_ele,:);
nodes = gcoord(sctr,:);
% xm=mean(xx);
% ym=mean(yy);
inc = 1;
EPS = 1e-10;
err = 1;
while (err > EPS)
    [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(c_ele,coord(1,1),coord(1,2),u_knot,v_knot,b_net);
    x = N'*nodes(:,1);
    y = N'*nodes(:,2);
    dx_dxi = dNdxi(:,1)' * nodes(:,1);
    dy_deta = dNdxi(:,2)' * nodes(:,1);
    dy_dxi = dNdxi(:,1)' * nodes(:,2);
    dx_deta = dNdxi(:,2)' * nodes(:,2);
 
    deltaX = x - c_point(1);
    deltaY = y - c_point(2);
    delta=[deltaX;deltaY];
    F=[dx_dxi dy_deta; dy_dxi dx_deta];
    invF=inv(F);
    
    ksi = ksi - invF(1,:)*delta;
    eta = eta - invF(2,:)*delta;
    
    coord(1) = ksi;
    coord(2) = eta;
    inc = inc + 1;
    err = sqrt((deltaX)^2+(deltaY)^2);
   
end
