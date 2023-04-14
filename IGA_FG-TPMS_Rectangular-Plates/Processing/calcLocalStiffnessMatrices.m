function [K,Kg]=calcLocalStiffnessMatrices(ngauss,nel,inn,ien,b_net,A,B,D,E,F,H,Ds,BucklingLoadCase)
global uKnot vKnot ndof sdof BucklingFactor

K=sparse(sdof,sdof);
Kg=sparse(sdof,sdof);
% ==== Get gaussian points and weights;
[gp,gw]=genGP_GW(ngauss);
nel_nza = 0; % Elements of non-zero area
tol=1e-8;
for iel = 1:nel
    
    sctr=ien(iel,:) ;         % element scatter vector
    nn=length(sctr);
    sctrBb(1:ndof:ndof*nn-4) = ndof.*sctr-4 ;
    sctrBb(2:ndof:ndof*nn-3) = ndof.*sctr-3 ;
    sctrBb(3:ndof:ndof*nn-2) = ndof.*sctr-2 ;
    sctrBb(4:ndof:ndof*nn-1) = ndof.*sctr-1 ;
    sctrBb(ndof:ndof:ndof*nn-0) = ndof.*sctr ;
    % Check to see if mlv current element has nonzero area;
    ni = inn(ien(iel,1),1); % Get NURBS coordinates
    nj = inn(ien(iel,1),2);
    % Element has positive area in the parametric domain
    if(abs(uKnot(ni)-uKnot(ni+1))>tol)&&(abs(vKnot(nj)-vKnot(nj+1))>tol)
        nel_nza = nel_nza + 1;
        % Used in calculating quadrature points. The factor of 4 comes from mapping from the [-1,1]
        % line onto a real segment...(in jacobian det)
        da =(uKnot(ni+1) - uKnot(ni))*(vKnot(nj+1) - vKnot(nj))/4;
        %--------------------------------------------------------
        % loop over integration points(ngauss in each direction);
        for igauss = 1: ngauss
            for jgauss = 1: ngauss
                [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(iel,gp(igauss),gp(jgauss),uKnot,vKnot,b_net);
                % calculate given element stiffness matrix and force vector;
                gwt = gw(igauss)*gw(jgauss)*da;
                % ==== Membrane matrix
                B1=zeros(3,ndof*nn);
                B1(1,1:ndof:ndof*nn-4)   =  dNdxy(:,1)';
                B1(2,2:ndof:ndof*nn-3)   =  dNdxy(:,2)';
                B1(3,1:ndof:ndof*nn-4)   =  dNdxy(:,2)';
                B1(3,2:ndof:ndof*nn-3)   =  dNdxy(:,1)';
                % ==== Bending matrix
                B2=zeros(3,ndof*nn);
                B2(1,3:ndof:ndof*nn-2)   = -dN2dxy(:,1)';
                B2(2,3:ndof:ndof*nn-2)   = -dN2dxy(:,2)';
                B2(3,3:ndof:ndof*nn-2)   = -2*dN2dxy(:,3)';
                
                B3=zeros(3,ndof*nn);
                B3(1,4:ndof:ndof*nn-1)   = dNdxy(:,1)';
                B3(2,ndof:ndof:ndof*nn)  = dNdxy(:,2)';
                B3(3,4:ndof:ndof*nn-1)   = dNdxy(:,2)';
                B3(3,ndof:ndof:ndof*nn)  = dNdxy(:,1)';
                % ==== Shear matrix 
                Bs=zeros(2,ndof*nn);
                Bs(1,4:ndof:ndof*nn-1)   = N';
                Bs(2,ndof:ndof:ndof*nn)  = N';
                % ==== Geometric matrix
                Bg=zeros(2,ndof*nn);
                Bg(1,3:ndof:ndof*nn-2)   = dNdxy(:,1)';
                Bg(2,3:ndof:ndof*nn-2)   = dNdxy(:,2)';
                % ==== Buckling load
                if BucklingLoadCase==0
                N0=BucklingFactor*[0 0;0 0]; % Initial stress matrix
                elseif BucklingLoadCase==1
                N0=BucklingFactor*[1 0;0 0]; % Initial stress matrix
                elseif BucklingLoadCase==2
                N0=BucklingFactor*[0 0;0 1]; % Initial stress matrix
                elseif BucklingLoadCase==3
                N0=BucklingFactor*[1 0;0 1]; % Initial stress matrix
                end
                %======
                k1=(B1'*A*B1+ B1'*B*B2+ B1'*E*B3)*gwt*detj;
                k2=(B2'*B*B1+ B2'*D*B2+ B2'*F*B3)*gwt*detj;
                k3=(B3'*E*B1+ B3'*F*B2+ B3'*H*B3)*gwt*detj;
                k4=Bs'*Ds*Bs*gwt*detj; 
                
                K(sctrBb,sctrBb)=K(sctrBb,sctrBb)+k1+k2+k3+k4;
                % ====
                Kg(sctrBb,sctrBb)=Kg(sctrBb,sctrBb)+ Bg'*N0*Bg*gwt*detj;
            end
        end
    end
end
