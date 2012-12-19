function ana_cfd(root_path)
% Refines CFD data for visualization and analysis of behavioral results
% This includes the creation of meshgridded flow data


%% Parameters

% Define boundaries of CFD volume for analysis (cm)
rangeX = [-1 2];
rangeY = [-1.5 1.5];
rangeZ = [-1.5 1.5];

% Number of values along x-axis
numX = 100;


%% Paths

% Query for root path, if not given
if nargin < 1
    root_path = uigetdir(pwd,'Select root directory (holds "cfd" & "behavior")');
end

% Paths to CFD data
cfd_path{1}  = [root_path filesep 'cfd' filesep 'flow_02cmps_around_zebrafish.mat'];
cfd_path{2}  = [root_path filesep 'cfd' filesep 'flow_11cmps_around_zebrafish.mat'];
cfd_path{3}  = [root_path filesep 'cfd' filesep 'flow_20cmps_around_zebrafish.mat'];

% Filenames used when saving data
fname{1} = 'flow_02cmps_reggrid';
fname{2} = 'flow_11cmps_reggrid';
fname{3} = 'flow_20cmps_reggrid';


%% Variables 

% Spacing between nodes 
dx = range(rangeX)/numX;

% Values along each dimension
xs = rangeX(1):dx:(rangeX(2)-dx);
ys = rangeY(1):dx:(rangeY(2)-dx);
zs = rangeZ(1):dx:(rangeZ(2)-dx);


%% Interpolate for each speed 

for i = 1:3;
    
    % Update user
    disp(' ');disp(['Working on ' fname{i} ' ...'])
    
    % Create x, y & z matrices
    [cR.x,cR.y,cR.z] = meshgrid(xs,ys,zs);
    
    % load cfd data ('c')
    load(cfd_path{1})
    
    %Index of CFD node values within volumetric range
    idx = ((c.x >= min(rangeX)) & (c.x <= max(rangeX))) & ...
        ((c.y >= min(rangeY)) & (c.y <= max(rangeY))) & ...
        ((c.z >= min(rangeZ)) & (c.z <= max(rangeZ)));
    
    
    % Interpolate velocity components
    warning off
    cR.u = griddata(c.x(idx),c.y(idx),c.z(idx),c.u(idx),cR.x,cR.y,cR.z);
    cR.v = griddata(c.x(idx),c.y(idx),c.z(idx),c.v(idx),cR.x,cR.y,cR.z);
    cR.w = griddata(c.x(idx),c.y(idx),c.z(idx),c.w(idx),cR.x,cR.y,cR.z);
    
    % Interpolate velocity gradients
    du_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dx(idx),cR.x,cR.y,cR.z);
    du_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dy(idx),cR.x,cR.y,cR.z);
    du_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dz(idx),cR.x,cR.y,cR.z);
    dv_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dx(idx),cR.x,cR.y,cR.z);
    dv_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dy(idx),cR.x,cR.y,cR.z);
    dv_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dz(idx),cR.x,cR.y,cR.z);
    dw_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dx(idx),cR.x,cR.y,cR.z);
    dw_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dy(idx),cR.x,cR.y,cR.z);
    dw_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dz(idx),cR.x,cR.y,cR.z);
    
    warning on
    
    % Flow speed
    cR.spd = sqrt(cR.u.^2 + cR.v.^2 + cR.w.^2);
    
    % Shear deformation
    cR.sh_def =  (2*(du_dx).^2 + 2*(dv_dy).^2 ...
        + 2*(dw_dz).^2 + (du_dy+dv_dx).^2 ...
        + (du_dz+dw_dx).^2 + (dv_dz+dw_dy).^2).^0.5;
    
    % Save data
    save([root_path filesep 'cfd' filesep fname{i}],'cR')
    
    clear cR idx
    
    disp('                                  ... Done!')
    
end

  

