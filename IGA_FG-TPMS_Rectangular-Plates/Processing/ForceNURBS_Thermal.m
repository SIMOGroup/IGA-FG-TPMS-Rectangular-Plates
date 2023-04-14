function [RHSF]=ForceNURBS_Thermal(ngauss,nel,inn,ien,b_net,N_th,M_th,P_th)
global u_knot v_knot nsd nshl sdof ndof

RHSF=zeros(sdof,1);
% get gaussian points and weights;
[gp,gw]=genGP_GW(ngauss);

area = 0;% Volume of the solid (for debugging)
nel_nza = 0;% Elements of non-zero area

%**************************************************************************
tol=1e-8;% in order to check u_knot(ni) matches u_knot(ni+1) or not ...
% loop over elements;
for iel = 1:nel
    sctr=ien(iel,:);           % element scatter vector
    nn=length(sctr);
    
    sctrBb(1:5:5*nn-4) = ndof.*sctr-4 ;
    sctrBb(2:5:5*nn-3) = ndof.*sctr-3 ;
    sctrBb(3:5:5*nn-2) = ndof.*sctr-2 ;
    sctrBb(4:5:5*nn-1) = ndof.*sctr-1 ;
    sctrBb(5:5:5*nn)   = ndof.*sctr ;
    
    %  check to see if mlv current element has nonzero area;
    ni = inn(ien(iel,1),1);% get NURBS coordinates
    nj = inn(ien(iel,1),2);
    % element has positive area in the parametric domain
    if(abs(u_knot(ni)-u_knot(ni+1))>tol)&&(abs(v_knot(nj)-v_knot(nj+1))>tol)
        nel_nza = nel_nza + 1;
        % used in calculating quadrature points. The factor of 4 comes from mapping from the [-1,1]
        % line onto a real segment...(in jacobian det)
        da =(u_knot(ni+1) - u_knot(ni))*(v_knot(nj+1) - v_knot(nj))/4;
        f_ther=zeros((nsd+3)*nshl,1);
        %--------------------------------------------------------
        % loop over integration points(ngauss in each direction);
        for igauss = 1: ngauss
            for jgauss = 1: ngauss
                [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
                % calculate given element stiffness matrix and force vector;
                gwt = gw(igauss)*gw(jgauss)*da;
                area = area + detj*gwt;
                
                B1=zeros(3,5*nn);
                B1(1,1:5:5*nn-4)   =  dNdxy(:,1)';
                B1(2,2:5:5*nn-3)   =  dNdxy(:,2)';
                B1(3,1:5:5*nn-4)   =  dNdxy(:,2)';
                B1(3,2:5:5*nn-3)   =  dNdxy(:,1)';
                
                B2=zeros(3,5*nn);
                B2(1,3:5:5*nn-2)   = -dN2dxy(:,1)';
                B2(2,3:5:5*nn-2)   = -dN2dxy(:,2)';
                B2(3,3:5:5*nn-2)   = -2*dN2dxy(:,3)';
                
                B3=zeros(3,5*nn);
                B3(1,4:5:5*nn-1)   = dNdxy(:,1)';
                B3(2,5:5:5*nn)     = dNdxy(:,2)';
                B3(3,4:5:5*nn-1)   = dNdxy(:,2)';
                B3(3,5:5:5*nn)     = dNdxy(:,1)';
                
                f_ther= f_ther +(B1'*N_th+B2'*M_th+B3'*P_th)*gwt*detj;
            end
        end
        RHSF(sctrBb)=RHSF(sctrBb) + f_ther;
    end
end
