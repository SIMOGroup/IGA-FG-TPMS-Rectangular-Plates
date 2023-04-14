% Isogeometric analysis for FG-TPMS Plates
% Written by Nguyen-Xuan, Huan Nguyen, Quoc K Tran

addpath Nurbs
addpath Processing
addpath Boundary

clear all; %close all;

global p q mcp ncp nsd nshl numnode nel;
global uKnot vKnot 
global gcoord L W sdof ndof h BucklingFactor
global p_min p_max type_TPMS Density_pattern

%===== Plate geometric properties ======
L = 1.5; W = 1;
% h = L/(1/0.01);
h = W/10;
p_min = 1.00; %relative density at bottom
p_max = 1.00; %relative density at top
%%========================== Density Pattern ================================ 
% 1: A 2: B
% Density_pattern = 2; n = 3.943;
Density_pattern = 1; n = 0;
%===========================Static loading ===============================
q_uniform = -1;
%=========================== Buckling Load type ===========================
BucklingLoadCase=1;
BucklingFactor=1;
%%========================== Types of TPMS ================================ 
% 0: FG Porous   1: Primitive    2: Gyroid    3: IWP
type_TPMS = 1;
%%============================ Material ====================================
imat = 1;
switch imat % [Ec Em nuC nuM rhoC rhoM]
    case 1  % Steel
        material=[200e9 0.3 8000];
    case 2  % Alumium
        material=[70e9  0.3 2702];
    case 3  % Titanium
        material=[106e9 0.3 4510];
    case 4  % Copper
        material=[117e9 0.3 8960];
    case 5  % Brass
        material=[90e9  0.3 8730];
end
%========================= Boundary condition=============================
% Option of bounary condition (ibc)
% 1: SSSS; 2: CCCCC; 3:
ibc = 1;
% Option of shear function (imix)
%[1]: gg=0  [2]: Reddy  [3]: Shimpi  [4]: H. Nguyen-Xuan  [5]: Hoang X. Nguyen
%[6]: Tuan N. Nguyen  [7]: Chien H. Thai   [8]:Huu-Tai Thai
imix = 5;

% ==========================Problem types==================================
% 1: Static    2: Vibration    3: Buckling
problemType = 3;

% ===== Extension, Bending and Shear stiffness matrix =====
%[Ab,Bb,Db,Eb,Fb,Hb,Ds,I,N_Th,M_Th,P_Th]=getMaterialThermalMatrices(imdis,h);
%[Ab,Bb,Db,Eb,Fb,Hb,Ds,I,N_Th,M_Th,P_Th]=getMaterial_Matrices(imdis,h);
[Ab,Bb,Db,Eb,Fb,Hb,Ds,I]=getMaterialMatrices(material,h,n,imix);

%=========================== NURBS functions ==============================
Deg = 3;    % Choose degree 
Ref = 11;    % Refinement
[CP,U,V,p,q]=Square_Coarse_Mesh(L,W,Deg); 
R1 = Refinement_vec(U,Ref); R2 = Refinement_vec(V,Ref);
[CP,uKnot,vKnot] = Knot_Refine_Surf(p,q,U,V,CP,R1,R2);

% ============================ Plot Mesh ===================================
% set(gcf,'color','white')
% Plot_NURBS_Surf(uKnot,vKnot,CP); hold on; axis equal
% plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); hold on
% plotNURBS_surf_El_CP_Note(p,q,uKnot,vKnot,CP); hold on
% plot_ctrlnet(CP,'bo');
% return
% Show node index
% count =0;
% for i=1: size(CP,1)
%    for j=1:size(CP,2) 
%        count=count +1;
%        X=CP(i,j,1);Y=CP(i,j,2);
%        text(X,Y,num2str(count));
%    end
% end


%===== Input Control Net =====
B_net(:,:,1)= CP(:,:,1); B_net(:,:,2)= CP(:,:,2); B_net(:,:,3)= CP(:,:,4);
mcp = length(B_net(:,1,1)); ncp = length(B_net(1,:,1));
Numx = mcp-p; Numy = ncp-q; 

% ===== Global Coordinate of Control Points ======
gcoord(:,1)=reshape(B_net(:,:,1),mcp*ncp,1);
gcoord(:,2)=reshape(B_net(:,:,2),mcp*ncp,1);
gcoord(:,3)=reshape(B_net(:,:,3),mcp*ncp,1);

% ===== Generate Connectivities =====    
[Ien,Inn]=Gen_IEN_INN_2D(p,q,mcp,ncp);

% ===== Number of elements, degrees of freedom
numnode=mcp*ncp;          % Number of control point
nshl = (p+1)*(q+1);     % Number of local shape functions
nel=(mcp-p)*(ncp-q);    % Number of element
nsd=2;                  % Number of spatial dimension
ndof=5; 
sdof=numnode*ndof;
ngauss= p+1;            % Number of gauss point in integration

% ======================= Building matrices ===============================    
[K,KG]=calcLocalStiffnessMatrices(ngauss,nel,Inn,Ien,B_net,Ab,Bb,Db,Eb,Fb,Hb,Ds,BucklingLoadCase);
M = calcLocalMassMatrix(ngauss,nel,Inn,Ien,B_net,I);
F=calLocalLoadVectorUniform(ngauss,nel,Inn,Ien,B_net,q_uniform); % extenal loading

%========================Imposing boundary conditions =====================
switch ibc
    case 1
        [bcdof, bcval]=Bcdof_SSSS(ndof,mcp,ncp); % Full simply supported
    case 2
        [bcdof, bcval]=Bcdof_CCCC(ndof,mcp,ncp); % CCCC
end

Em = material(1); rhom = material(3); num = material(2);

switch problemType % [Static  Vibration  Buckling]
    case 1  % Static
        freedof = setdiff((1:ndof * (numnode))', bcdof');
        U = zeros(ndof * numnode, 1); 
        U(bcdof') = bcval';
        F(freedof) = F(freedof) - K(freedof, bcdof') * bcval';
        U(freedof) = K(freedof, freedof) \ F(freedof);
        
        c_ele    = floor((nel+1)/2); % central element
        c_point  = [L/2,W/2]; % coord of center 
        %Calculate xi, eta value (in parametric space) of central point 
        coord=Calculate_xieta(c_ele,c_point,Ien,gcoord,uKnot,vKnot,B_net);
        
        %Central deflection
        sctr = Ien(c_ele,:);
        sctrW = ndof.*sctr-2;
        [N,dNdxi,dNdxy,dN2dxy,detj]=Kine_Shape_2nd(c_ele,coord(1,1),coord(1,2),...
                                                                   uKnot,vKnot,B_net);
        
       def_cen=N'*U(sctrW);
       %Normalised
%        h = h/0.35;
       Dm = Em*h^3/(12*(1-num^2));
       % Normalisation
       cen_def_norm = def_cen*100*Dm/(q_uniform*L^4)
%        size_factor =3;
%        count =0;
%         for i=1: size(CP,1)
%           for j=1:size(CP,2) 
%             count=count +1;
%             CP(i,j,3)=CP(i,j,3)+ size_factor*sign(q_uniform)*U(ndof*count-2)*100*Dm/(q_uniform*L^4);
%           end
%         end

%        plotNURBS_surf_El_CP(p,q,uKnot,vKnot,CP); 
       
     case 2  % Vibration
%         h = h/0.35;
        KK=full(K);MM=full(M);
        norm1 = L^2/h*sqrt(rhom/Em);
        norm2 = h*sqrt(rhom/Em);
        D = Em*h^3/12/(1-num^2);
        Lamda0 = Eigen(KK,MM,bcdof);
        Lamda0 = sort(Lamda0);
%         Lamda0(1)^0.25
        Lamda0 = (Lamda0(1:2)*rhom*L^4*h/D).^0.25
%         Lamda0 = sqrt(Lamda0(1:1,1))*norm2
    case 3  % Buckling
        % Buckling analysis
        KK=full(K); KKG=full(KG);
        [Lamda,ModeShape]=Eigen(KK,KKG,bcdof);
        Pcr_Temp=Lamda(1:5,1);
        for ii=1:length(Pcr_Temp)
            if Pcr_Temp(ii)>0
              Pcr=Pcr_Temp(ii);
             break
            end
        end
%         h = h/0.35;
        D = Em*h^3/12/(1-num^2);
        norm = W^2/(pi^2*D);
%         norm = L^2/(Em*h^3);
        PcrX = Pcr*norm
 end


