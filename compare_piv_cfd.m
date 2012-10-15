function compare_piv_cfd(pred_spd)
% Relates cfd data to piv data for a particular predator speed
%
% pred_spd - integer with value of 2, 11 or 20
% 
% Note that the cfd data are massive.  It therefore makes sense to crop the
% data and/or downsample prior to many analyses or visualizations.

% Approximate number of nodes in each dimension
num_nodes = 50;

% Default speed
if nargin < 1
    pred_spd = 2;
end

% Downsample interval 
down_int = 1;


y_range = 0.1;


%% Define paths 

% % Root paths for piv & cfd data
cfd_root = '/Volumes/Docs/Dropbox/Projects/Robopredator/cfd';
piv_root = '/Volumes/Docs/Dropbox/Projects/Robopredator/piv';
%cfd_root = '/Users/mmchenry/Dropbox/Projects/Robopredator/cfd';
%piv_root = '/Users/mmchenry/Dropbox/Projects/Robopredator/piv';

% Define subdirectories 
if pred_spd == 2
    
    cfd_path = [cfd_root filesep 'flow_02cmps_around_zebrafish.ascii'];
    piv_path = [piv_root filesep 'compiled data' filesep 'S02_piv.mat'];
    
elseif pred_spd == 11
    
    cfd_path = [cfd_root filesep 'flow_11cmps_around_zebrafish.ascii'];
    piv_path = [piv_root filesep 'compiled data' filesep 'S11_piv.mat'];
    
elseif pred_spd == 20
    
    cfd_path = [cfd_root filesep 'flow_20cmps_around_zebrafish.ascii'];
    piv_path = [piv_root filesep 'compiled data' filesep 'S20_piv.mat'];
    
else
    
    error('predator speed must be defined as 2, 11 or 20')
    
end


%% Load data

% Load cfd velocity data
c = cfd_import_vel(cfd_path);

% Load piv data ('p')
load(piv_path)

% Index that defines domain covered by piv data
i_piv = (c.x>min(p.x)) & (c.x<max(p.x)) & ...
        (c.y>min(p.y)) & (c.y<max(p.y)) & ...
        (c.z>min(p.z)) & (c.z<max(p.z));
    
% Trim the cfd data
x_cfd = c.x(i_piv);
y_cfd = c.y(i_piv);
z_cfd = c.z(i_piv);
u_cfd = c.u(i_piv);
v_cfd = c.v(i_piv);
w_cfd = c.w(i_piv);

clear c i_piv


%% Compare on the saggital plane

% Index for saggital plane values
i_saggital = isnan(p.v) & (p.y==0);

% PIV coordinates and vectors for the saggital plane
x = p.x(i_saggital);
y = p.y(i_saggital);
z = p.z(i_saggital);

u_piv = p.u(i_saggital);
w_piv = p.w(i_saggital);

% Downsample the piv data
x = x(1:down_int:end);
y = y(1:down_int:end);
z = z(1:down_int:end);
u_piv = u_piv(1:down_int:end);
w_piv = w_piv(1:down_int:end);

% % Trim the range of y-values for the CFD data
% i_range = (y_cfd<(max(y_cfd)*y_range)) & (y_cfd>(min(y_cfd)*y_range));
% 
% x_cfd = x_cfd(i_range);
% y_cfd = y_cfd(i_range);
% z_cfd = z_cfd(i_range);
% u_cfd = u_cfd(i_range);
% v_cfd = v_cfd(i_range);
% w_cfd = w_cfd(i_range);
% 
% clear i_range

% Define nodes at which to evaluate
xNode = linspace(min(x),max(x),num_nodes);

% Interval between nodes along x-axis
dNode = mean(diff(xNode));

% Evenly divide by same interval along y-axis
%yNode = min(y):dNode:max(y);

% Ditto for the z-axis
zNode = min(z):dNode:max(z);

[X_piv,Z_piv] = meshgrid(xNode,zNode);

Up_func = TriScatteredInterp(x,z,u_piv);
Wp_func = TriScatteredInterp(x,z,w_piv);

U_piv = Up_func(X_piv,Z_piv);
W_piv = Wp_func(X_piv,Z_piv);

%[X_cfd,Y_cfd,Z_cfd] = meshgrid(xNode,xNode.*0,zNode);

warning off
Uc_func = TriScatteredInterp(x_cfd,y_cfd,z_cfd,u_cfd);
Vc_func = TriScatteredInterp(x_cfd,y_cfd,z_cfd,v_cfd);
Wc_func = TriScatteredInterp(x_cfd,y_cfd,z_cfd,w_cfd);
warning on

U_cfd = Uc_func(X_piv,X_piv.*0,Z_piv);
%V_cfd = Vc_func(X_piv,X_piv.*0,Z_piv);
W_cfd = Wc_func(X_piv,X_piv.*0,Z_piv);


figure
subplot(1,2,1)
quiver(X_piv,Z_piv,U_piv,W_piv);
axis equal

subplot(1,2,2)
quiver(X_piv,Z_piv,U_cfd,W_cfd);
axis equal

return

[x,z,u_piv] = meshgrid(x,z,u_piv);


U = interp2(x,z,u_piv,X,Z);

%TODO: use coordinates from piv to find interpolated values from CFD
%dataset





%u_cfd = interp3(d.x,d.y,d.z,d.u,xNode,yNode,zNode);


subplot(2,2,1)
quiver(x,z,u_piv,w_piv)





[X,Y,Z] = meshgrid(p.x,p.y,p.z);

[U,V,W] = meshgrid(p.u,p.v,p.w);

%% Define nodes to evaluate cfd and piv data

% Evenly divide along x-axis
xNode = linspace(min(p.x):max(p.x),num_nodes);





clear dNode

[Xnode,Ynode,Znode] = meshgrid(xNode,yNode,zNode);


% Evaluate 



ttt=3;


