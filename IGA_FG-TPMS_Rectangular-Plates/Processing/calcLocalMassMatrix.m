function M=calcLocalMassMatrix(ngauss,nel,inn,ien,b_net,Inetia)
global uKnot vKnot sdof

M=sparse(sdof,sdof);
% get gaussian points and weights;
[gp,gw]=genGP_GW(ngauss);
%**************************************************************************
I1=Inetia(1);I2=Inetia(2);I3=Inetia(3);
I4=Inetia(4);I5=Inetia(5);I6=Inetia(6);
tol=1e-8;% in order to check u_knot(ni) matches u_knot(ni+1) or not ...
% loop over elements;
for iel = 1:nel
    sctr=ien(iel,:);           % element scatter vector
    nn=length(sctr);
    sctrBb(1:5:5*nn-4) = 5.*sctr-4 ;
    sctrBb(2:5:5*nn-3) = 5.*sctr-3 ;
    sctrBb(3:5:5*nn-2) = 5.*sctr-2 ;
    sctrBb(4:5:5*nn-1) = 5.*sctr-1 ;
    sctrBb(5:5:5*nn)   = 5.*sctr ;
    
    %  check to see if mlv current element has nonzero area;
    ni = inn(ien(iel,1),1);% get NURBS coordinates
    nj = inn(ien(iel,1),2);
    % element has positive area in the parametric domain
    if(abs(uKnot(ni)-uKnot(ni+1))>tol)&&(abs(vKnot(nj)-vKnot(nj+1))>tol)
        da =(uKnot(ni+1) - uKnot(ni))*(vKnot(nj+1) - vKnot(nj))/4;   
        %--------------------------------------------------------
        % loop over integration points(ngauss in each direction);      
        for igauss = 1: ngauss
            for jgauss = 1: ngauss
                [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(iel,gp(igauss),gp(jgauss),uKnot,vKnot,b_net);
                % calculate given element stiffness matrix and force vector;
                gwt = gw(igauss)*gw(jgauss)*da;
                
                N0_shape=zeros(3,5*nn);
                N0_shape(1,1:5:5*nn-4)= N';
                N0_shape(2,2:5:5*nn-3)= N';
                N0_shape(3,3:5:5*nn-2)= N';
                
                N1_shape=zeros(3,5*nn);
                N1_shape(1,3:5:5*nn-2)= -dNdxy(:,1)';
                N1_shape(2,3:5:5*nn-2)= -dNdxy(:,2)';
                
                N2_shape=zeros(3,5*nn);
                N2_shape(1,4:5:5*nn-1)= N';
                N2_shape(2,5:5:5*nn)  = N';
                
                M(sctrBb,sctrBb)=M(sctrBb,sctrBb) +(I1*N0_shape'*N0_shape+ I2*N0_shape'*N1_shape+ I4*N0_shape'*N2_shape...
                                                  + I2*N1_shape'*N0_shape+ I3*N1_shape'*N1_shape+ I5*N1_shape'*N2_shape...
                                                  + I4*N2_shape'*N0_shape+ I5*N2_shape'*N1_shape+ I6*N2_shape'*N2_shape)*gwt*detj;
            end
        end
    end
end
