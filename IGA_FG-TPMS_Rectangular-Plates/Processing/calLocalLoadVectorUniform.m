function RHSF=calLocalLoadVectorUniform(ngauss,nel,inn,ien,b_net,q_uniform)
global uKnot vKnot sdof

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
  %nn=length(sctr);
  sctrF = 5.*sctr-2 ;
  %  check to see if mlv current element has nonzero area;
  ni = inn(ien(iel,1),1);% get NURBS coordinates
  nj = inn(ien(iel,1),2);
% element has positive area in the parametric domain
  if(abs(uKnot(ni)-uKnot(ni+1))>tol)&&(abs(vKnot(nj)-vKnot(nj+1))>tol)
     nel_nza = nel_nza + 1;
     % used in calculating quadrature points. The factor of 4 comes from mapping from the [-1,1]
     % line onto a real segment...(in jacobian det)
     da =(uKnot(ni+1) - uKnot(ni))*(vKnot(nj+1) - vKnot(nj))/4;  
     %--------------------------------------------------------
     % loop over integration points(ngauss in each direction);
     
     for igauss = 1: ngauss
        for jgauss = 1: ngauss
            [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(iel,gp(igauss),gp(jgauss),uKnot,vKnot,b_net);
            % calculate given element stiffness matrix and force vector;
              gwt = gw(igauss)*gw(jgauss)*da;
              area = area + detj*gwt;
              RHSF(sctrF)=RHSF(sctrF)+N*q_uniform*gwt*detj; %stiffness matrix
              clear N dNdxi dNdxy detj 
         end
     end % end integration loop for bending and mass matrix    
  end %end if
end % end loop over elements
return
