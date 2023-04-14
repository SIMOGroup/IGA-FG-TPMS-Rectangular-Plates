function [bcdof, bcval]=Bcdof_FullySupportCircularPlates(ndof,mcp, ncp)

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

bcdof = unique([bcdof ndof*right1-2 ndof*left1-2 ndof*upper1-2 ndof*lower1-2]);     % w - soft supported boundary

% bcdof = unique([bcdof ndof*right1-4 *right1-3 ndof*right1-2 ... 
%     ndof*left1-4 ndof*left1-3 ndof*left1-2 ...
%     ndof*upper1-4 ndof*upper1-3 ndof*upper1-2 ...
%     ndof*lower1-4 ndof*upper1-3 ndof*lower1-2]);     % w - soft supported boundary

% bcdof = unique([bcdof ndof*right1-4 *right1-3 ndof*right1-2 ... 
%     ndof*left1-4 ndof*left1-3 ndof*left1-2 ...
%     ndof*upper1-4 ndof*upper1-3 ndof*upper1-2 ...
%     ndof*lower1-4 ndof*upper1-3 ndof*lower1-2]);     % w - soft supported boundary
bcval = [bcval, zeros(1,length(bcdof))];
return