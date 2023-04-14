function [RHSF,nel_nza]=FmatNURBS_Sin(ngauss,nel,inn,ien,b_net,P,L,nodeCoordinates)
global u_knot v_knot sdof

RHSF=zeros(sdof,1);
% get gaussian points and weights;
[gp,gw]=genGP_GW(ngauss);

area = 0;% Volume of the solid (for debugging)
nel_nza = 0;% Elements of non-zero area
tol=1e-8;% in order to check u_knot(ni) matches u_knot(ni+1) or not ...
% loop over elements;
for iel = 1:nel
    sctr=ien(iel,:);           % element scatter vector
    nn=length(sctr);
    
    sctrF = 5.*sctr-2 ;
    %  check to see if mlv current element has nonzero area;
    ni = inn(ien(iel,1),1);% get NURBS coordinates
    nj = inn(ien(iel,1),2);
    % element has positive area in the parametric domain
    if(abs(u_knot(ni)-u_knot(ni+1))>tol)&&(abs(v_knot(nj)-v_knot(nj+1))>tol)
        nel_nza = nel_nza + 1;
        % used in calculating quadrature points. The factor of 4 comes from mapping from the [-1,1]
        % line onto a real segment...(in jacobian det)
        da =(u_knot(ni+1) - u_knot(ni))*(v_knot(nj+1) - v_knot(nj))/4;
        
        %--------------------------------------------------------
        % loop over integration points(ngauss in each direction);
        
        for igauss = 1: ngauss
            for jgauss = 1: ngauss
                [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(iel,gp(igauss),gp(jgauss),u_knot,v_knot,b_net);
                % calculate given element stiffness matrix and force vector;
                gwt = gw(igauss)*gw(jgauss)*da;
                x=nodeCoordinates(sctr,1)'*N;
                y=nodeCoordinates(sctr,2)'*N;
                S1=sin(pi*x/L);
                S2=sin(pi*y/L);
                RHSF(sctrF)=RHSF(sctrF)+N*P*S2*S1*gwt*detj; %stiffness matrix
                clear N dNdxi dNdxy detj
            end
        end
    end
end