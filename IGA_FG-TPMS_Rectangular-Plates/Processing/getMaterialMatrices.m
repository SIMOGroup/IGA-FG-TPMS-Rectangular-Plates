function [Ab,Bb,Db,Eb,Fb,Hb,Ds,I]=getMaterialMatrices(matrial,h,n,imix)

global p_min p_max type_TPMS Density_pattern

% ===== Material matrix =====
Ab=zeros(3,3); Bb=Ab; Db=Ab; Eb=Ab; Fb=Ab; Hb=Ab;
Ds=zeros(2,2);I=zeros(1,6);

% % ============ material property============================
% get integration points
%[Wg,Qg] = quadrature(8,'GAUSS',1)

x1 = 1;
x2 = -1 ;
[Qg,Wg]=gauleg(x1,x2,20) ;
zu = h/2;              % upper
zl = -h/2;             % lower

% volume fraction of TPMS
%Vc_up = 1; Vc_lo = 0;

% loop over the integration points
for igp = 1:size(Wg,1)
    pt = Qg(igp,:) ;
    
    % map the point to global
    gpt = (zu+zl)/2 + (zu-zl)/2*pt ;                 % value of z-direction
    wt = (zu-zl)/2*Wg(igp) ;
    
    % Young's modulus, poisson's ratio
    E = matrial(1);
    nu = matrial(2);     % nu: constant (0.3)
    G = E/(2*(1+nu));
    rho = matrial(3);
    
    switch Density_pattern
        case 1
            p_z = p_min + (p_max - p_min)*(gpt/h + 1/2)^n; % Eq (1) in Viet and Zaki Composite Structures 274 (2021) 114342
        case 2
            p_z = p_min + (p_max - p_min)*(1 - cos(pi*gpt/h))^n;
    end
        
    switch type_TPMS
        case 100 %Cellular solid
            E_z = E*p_z^2; G_z = G*p_z^2; nu_z = nu;
        case 0 %FG Porous %HS Upper bound
            K = E/(3*(1 - 2*nu));
            K_z = (4*G*K*p_z) / (4*G + 3*K*(1-p_z));
            G_z = ((9*K + 8*G)*G*p_z) / (20*G + 15*K - 6*(K+2*G)*p_z);
        case 1 %Primitive
            k_e = 0.25; C_1e = 0.317; n_1e = 1.264; n_2e = 2.006;
            k_g = 0.25; C_1g = 0.705; n_1g = 1.189; n_2g = 1.715;
            k_nu = 0.55; a_1 = 0.314; b_1 = -1.004; a_2 = 0.152;
        case 2 %Gyroid
            k_e = 0.45; C_1e = 0.596; n_1e = 1.467; n_2e = 2.351;
            k_g = 0.45; C_1g = 0.777; n_1g = 1.544; n_2g = 1.982;
            k_nu = 0.50; a_1 = 0.192; b_1 = -1.349; a_2 = 0.402;
        case 3 %IWP
            k_e = 0.35; C_1e = 0.597; n_1e = 1.225; n_2e = 1.782;
            k_g = 0.35; C_1g = 0.529; n_1g = 1.287; n_2g = 2.188;
            k_nu = 0.13; a_1 = 2.597; b_1 = -0.157; a_2 = 0.201;
    end
    
    switch type_TPMS
        case 0
            E_z = (9*G_z*K_z) / (3*K_z + G_z);
            nu_z = (3*K_z - 2*G_z) / (2*(3*K_z + G_z));
        case {1,2,3}
            C_2e = (C_1e*k_e^(n_1e) - 1)/(k_e^(n_2e) - 1); C_3e = 1 - C_2e;
            C_2g = (C_1g*k_g^(n_1g) - 1)/(k_g^(n_2g) - 1); C_3g = 1 - C_2g;
            d_1 = nu - a_1*exp(b_1*k_nu); b_2 = - a_2*(k_nu + 1); d_2 = nu - a_2*(1)^2 - b_2(1);
            
            e_z = (p_z <= k_e)*(C_1e*p_z^(n_1e)) + ...
                  (p_z >  k_e)*(C_2e*p_z.^(n_2e) + C_3e);
            g_z = (p_z <= k_g)*(C_1g*p_z^(n_1g)) + ...
                  (p_z >  k_g)*(C_2g*p_z.^(n_2g) + C_3g);
            nu_z =(p_z <= k_nu)*(a_1.*exp(b_1*p_z) + d_1) + ...
                  (p_z >  k_nu)*(a_2*p_z^(2) + b_2*p_z + d_2); 
            
            E_z = E*e_z;
            G_z = G*g_z;
    end

    rho_z = rho*p_z;

    % % ============= Bending-extension stiffness matrix ===================
    Qb = [E_z/(1-nu_z^2) E_z*nu_z/(1-nu_z^2) 0;
          E_z*nu_z/(1-nu_z^2) E_z/(1-nu_z^2) 0;
          0 0 G_z];
    % shear stiffness
    Qs =G_z*[1 0;
             0 1];
    % choos model for displacement field
    switch imix
        case 1
            %E. Reissner, The effect of transverse shear deformation on the bending of elastic plates
            gg = 0;% 3DOF with no higher-order term %
            dgg= 0;
        case 2
            %J. N. Reddy, A Simple Higher-Order Theory for Laminated Composite Plates,
            gg = -4*gpt^3/(3*h^2);% Reddy
            dgg= -4*gpt^2/h^2;
        case 3
            %R. P. Shimpi, Refined Plate Theory and Its Variants
            gg = -5*gpt^3/(3*h^2)+gpt/4;%Shimpi
            dgg=-5*gpt^2/h^2+1/4;
        case 4
            %H. Nguyen-Xuan, Isogeometric finite element analysis of composite sandwich plates using a higher order shear deformation theory
            gg = -1/8*gpt-2/h^2*gpt^3+2/h^4*gpt^5;% H. Nguyen-Xuan
            dgg= -1/8-6/h^2*gpt^2+10/h^4*gpt^4;
        case 5
            % H. X. Nguyen, T. N. Nguyen, A refined quasi-3D isogeometric analysis for functionally graded microplates based on the modified couple stress theory
            gg = -9*gpt+(10/h^2)*gpt^3+6/(5*h^4)*gpt^5+8/(7*h^6)*gpt^7;%Hoang X. Nguyen
            dgg=-9+(10/h^2)*3*gpt^2+6/(5*h^4)*5*gpt^4+8/(7*h^6)*7*gpt^6;
        case 6
            %T. N. Nguyen, On the general framework of high order shear deformation theories for laminated composite plate structures: A novel unified approach
            gg = - (17*gpt^3)/(10*h^2) + (22*gpt^5)/(25*h^4);%Tuan N. Nguyen
            dgg= - (17*3*gpt^2)/(10*h^2) + (22*5*gpt^4)/(25*h^4);
        case 7
            %C. H. Thai, Generalized shear deformation theory for functionally graded isotropic and sandwich plates based on isogeometric approach
            gg = atan(sin(pi/h*gpt))-gpt;%Chien H. Thai (Computers & Structures)
            dgg= pi/h*cos(pi/h*gpt)/(1+(sin(pi/h*gpt))^2)-1;
        otherwise
            disp('**************************************************');
            disp('Do not appropriate models ');
            disp('**************************************************');
            pause
    end
    
    Ab = Ab + Qb*wt ;        Bb = Bb + Qb*gpt*wt ;     Db = Db + Qb*gpt^2*wt ;
    Eb = Eb + Qb*gg*wt ;     Fb = Fb + Qb*gg*gpt*wt ;  Hb = Hb + Qb*gg^2*wt ;
    Ds = Ds + Qs*(dgg)^2*wt;         
    % ================= Inertia terms matrix =======================
    I = I + (rho_z*wt).* [1 gpt gpt^2 gg gg*gpt gg^2];
end