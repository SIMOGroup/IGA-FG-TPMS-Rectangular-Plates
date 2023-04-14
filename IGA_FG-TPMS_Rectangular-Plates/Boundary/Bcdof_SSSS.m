function [bcdof, bcval]=Bcdof_SSSS(ndof,mcp, ncp)

bcdof=[];
bcval=[];
% u, v, w, betax, betay
right1=[mcp:mcp:mcp*ncp];               % on perimeter
right2=[mcp-1:mcp:mcp*ncp-1];           % inside

left1=[1:mcp:mcp*ncp-mcp+1];            % on perimeter
left2=[2:mcp:mcp*ncp-mcp+2];            % inside

lower1=[1:mcp];                         % on perimeter
lower2=[mcp+1:2*mcp];                   % inside

upper1=[mcp*ncp-mcp+1:mcp*ncp];         % on perimeter
upper2=[mcp*ncp-2*mcp+1:mcp*ncp-mcp];   % inside

% right 
% bcdof = [bcdof ndof*right1-4];     % u0 (Immovable edge (SSSS2))
bcdof = [bcdof ndof*right1-3];     % v0
bcdof = [bcdof ndof*right1-2];     % W
bcdof = [bcdof ndof*right1];       % betay

% left 
% bcdof = [bcdof ndof*left1-4];     % u0 (Immovable edge (SSSS2))
bcdof = [bcdof ndof*left1-3];     % v0
bcdof = [bcdof ndof*left1-2];     % W
bcdof = [bcdof ndof*left1];       % betay

% upper 
bcdof = [bcdof ndof*upper1-4];     % u0
% bcdof = [bcdof ndof*upper1-3];      % v0 (Immovable edge (SSSS2))
bcdof = [bcdof ndof*upper1-2];     % W
bcdof = [bcdof ndof*upper1-1];     % betax

% lower 
bcdof = [bcdof ndof*lower1-4];     % u0
% bcdof = [bcdof ndof*lower1-3];      % v0 (Immovable edge (SSSS2))
bcdof = [bcdof ndof*lower1-2];     % W
bcdof = [bcdof ndof*lower1-1];     % betax

bcval = [bcval, zeros(1,length(bcdof))];
return