function ana_body_flow(root_path)
% Analyzes the flow along the body of larvae from body position data in 
% behavioral experiments and flow data from CFD simulations.
%
% root_path - the base directory that contains 'cfd' and 'behavior' directories 


%% Parameters

% Period between detection and fast start initiation (s)
latency = 5e-3;

% Number of points along the body
numPoints = 100;


%% Paths

if nargin < 1
    root_path = uigetdir(pwd,'Select root directory (holds "cfd" & "behavior")');
end

% Load latest behavior coordinates ('b' structure)
load([root_path filesep 'behavior' filesep 'Transformed_Prey_Coords.mat'])

% Paths to CFD data
cfd2_path   = [root_path filesep 'cfd' filesep 'flow_02cmps_around_zebrafish.mat'];
cfd11_path  = [root_path filesep 'cfd' filesep 'flow_11cmps_around_zebrafish.mat'];
cfd20_path  = [root_path filesep 'cfd' filesep 'flow_20cmps_around_zebrafish.mat'];


%% Set up variables

% Store latency assumed in the analysis
f.latency = latency;

% Flow speed
f.spd      = nan(size(b.preyx,1),numPoints);

% Velocity gradient along body
f.velgrad = nan(size(b.preyx,1),numPoints);

% Shear deformation
f.shrdef    = nan(size(b.preyx,1),numPoints);

% Body position (arclength)
f.s         = nan(size(b.preyx,1),numPoints);

% Body position (inertial FOR)
f.xbod      = nan(size(b.preyx,1),numPoints);
f.ybod      = nan(size(b.preyx,1),numPoints);
f.zbod      = nan(size(b.preyx,1),numPoints);

% Number of sequences
num_seq = length(b.preyx(:,1));


%% Interpolate for 2 cm/s sequences

% Import 2 cm/s cfd data ('c' structure)
load(cfd2_path);

% Index for 2 cm/s
idx = find((b.speed(1:num_seq)==2) & (b.LL(1:num_seq)==1) ...
           & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0));

for i = 1:length(idx)
    
    % Sequence number
    seqnum = idx(i);
    
    % Offset in x, due to latency
    lat_offset = latency * 2;
 
    % Find flow conditions
    flow = interpflow(c,b,seqnum,numPoints,lat_offset);
    
    % Store results
    f.s(seqnum,:)        = flow.s;
    f.spd(seqnum,:)      = flow.spd;
    f.shrdef(seqnum,:)   = flow.sh_def;
    f.velgrad(seqnum,:)  = flow.vel_grad;
    f.xbod(seqnum,:)     = flow.bod_x;
    f.ybod(seqnum,:)     = flow.bod_y;
    f.zbod(seqnum,:)     = flow.bod_z;
    
    clear flow x_vals y_vals z_vals seqnum
end

clear c idx

disp(' ');disp('Done 2 cm/s')


%% Interpolate for 11 cm/s sequences

% Import 11 cm/s cfd data ('c' structure)
load(cfd11_path);

% Index for 11 cm/s
idx = find((b.speed(1:num_seq)==11) & (b.LL(1:num_seq)==1) ...
           & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0));

for i = 1:length(idx)
    
    % Sequence number
    seqnum = idx(i);
    
    % Offset in x, due to latency
    lat_offset = latency * 11;
 
    % Find flow conditions
    flow = interpflow(c,b,seqnum,numPoints,lat_offset);
    
    f.s(seqnum,:)        = flow.s;
    f.spd(seqnum,:)      = flow.spd;
    f.shrdef(seqnum,:)   = flow.sh_def;
    f.velgrad(seqnum,:)  = flow.vel_grad;
    f.xbod(seqnum,:)     = flow.bod_x;
    f.ybod(seqnum,:)     = flow.bod_y;
    f.zbod(seqnum,:)     = flow.bod_z;
    
    clear flow x_vals y_vals z_vals seqnum
end

clear c idx

disp(' ');disp('Done 11 cm/s')


%% Interpolate for 20 cm/s sequences

% Import 2 cm/s cfd data ('c' structure)
load(cfd20_path);

% Index for 2 cm/s
idx = find((b.speed(1:num_seq)==20) & (b.LL(1:num_seq)==1) ...
           & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0));

for i = 1:length(idx)
    
    % Sequence number
    seqnum = idx(i);
    
    % Offset in x, due to latency
    lat_offset = latency * 20;
 
    % Find flow conditions
    flow = interpflow(c,b,seqnum,numPoints,lat_offset);
    
    f.s(seqnum,:)        = flow.s;
    f.spd(seqnum,:)      = flow.spd;
    f.shrdef(seqnum,:)   = flow.sh_def;
    f.velgrad(seqnum,:)  = flow.vel_grad;
    f.xbod(seqnum,:)     = flow.bod_x;
    f.ybod(seqnum,:)     = flow.bod_y;
    f.zbod(seqnum,:)     = flow.bod_z;
    
    clear flow x_vals y_vals z_vals seqnum
end

clear c idx

disp(' ');disp('Done 20 cm/s')


%% Save data

save([root_path filesep 'behavior' filesep 'flow_along_body.mat'],'f')



function out = interpflow(c,b,seqnum,numPoints,lat_offset)
% Interpolates CFD flow field to find flow conditions at input coordinates
% c - structure of CFD data
% b - structure of 

% Scaling of body length for CFD volume to interpolate over
sclfactr = 1.5;

% Body length of larva
blength = norm([b.preyx(seqnum,3)-b.preyx(seqnum,1) ...
                b.preyy(seqnum,3)-b.preyy(seqnum,1) ...
                b.preyz(seqnum,3)-b.preyz(seqnum,1)]);

% Origin of local FOR
b_origin(1,1) = b.preyx(seqnum,1);
b_origin(1,2) = b.preyy(seqnum,1);
b_origin(1,3) = b.preyz(seqnum,1);

% Retrieve local x axis to determine coordinate system
b_xaxis(1,1) = b.preyx(seqnum,3) - b_origin(1);
b_xaxis(1,2) = b.preyy(seqnum,3) - b_origin(2);
b_xaxis(1,3) = b.preyz(seqnum,3) - b_origin(3);

% Normalize to create a unit vector
b_xaxis = b_xaxis./norm(b_xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
b_yaxis = [-b_xaxis(2) b_xaxis(1) 0];

% Normalize to create a unit vector
b_yaxis = b_yaxis./norm(b_yaxis);

%Determine local z axis
b_zaxis = cross(b_xaxis,b_yaxis);

% Normalize to create a unit vector
b_zaxis = b_zaxis./norm(b_zaxis);

%Create rotation matrix (from inertial axes to local axes)
% [x y z] * R       - inertial to local
% [x y z] * inv(R)  - local to inertial
R = [b_xaxis' b_yaxis' b_zaxis'];

% Body position values
s = [blength.*linspace(0,1,numPoints)]';

% Body position rotated in inertial FOR
bod = [s s.*0 s.*0] * inv(R);

% Translate in inertial FOR
bod_x = bod(:,1) + b_origin(1) + lat_offset;
bod_y = bod(:,2) + b_origin(2);
bod_z = bod(:,3) + b_origin(3);

clear bod b_xaxis b_yaxis b_zaxis b_origin 

% Mean position of body
meanX = mean([b.preyx(seqnum,1) b.preyx(seqnum,3)]) + lat_offset;
meanY = mean([b.preyy(seqnum,1) b.preyy(seqnum,3)]);
meanZ = mean([b.preyz(seqnum,1) b.preyz(seqnum,3)]);

% Range of volume to interrogate
rangeX = [meanX-blength*sclfactr meanX+blength*sclfactr];
rangeY = [meanY-blength*sclfactr meanY+blength*sclfactr];
rangeZ = [meanZ-blength*sclfactr meanZ+blength*sclfactr];

%Index of CFD node values within volumetric range
idx = ((c.x >= min(rangeX)) & (c.x <= max(rangeX))) & ...
      ((c.y >= min(rangeY)) & (c.y <= max(rangeY))) & ...
      ((c.z >= min(rangeZ)) & (c.z <= max(rangeZ)));

% griddata raises pointless alarms  
warning off

% Interpolate velocity components
u = griddata(c.x(idx),c.y(idx),c.z(idx),c.u(idx),bod_x,bod_y,bod_z);
v = griddata(c.x(idx),c.y(idx),c.z(idx),c.v(idx),bod_x,bod_y,bod_z);
w = griddata(c.x(idx),c.y(idx),c.z(idx),c.w(idx),bod_x,bod_y,bod_z);

% Interpolate velocity gradients
du_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dx(idx),bod_x,bod_y,bod_z);
du_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dy(idx),bod_x,bod_y,bod_z);
du_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dz(idx),bod_x,bod_y,bod_z);
dv_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dx(idx),bod_x,bod_y,bod_z);
dv_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dy(idx),bod_x,bod_y,bod_z);
dv_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dz(idx),bod_x,bod_y,bod_z);
dw_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dx(idx),bod_x,bod_y,bod_z);
dw_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dy(idx),bod_x,bod_y,bod_z);
dw_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dz(idx),bod_x,bod_y,bod_z);

% Turn warnings back on
warning on

% Calculate flow in local FOR
for j = 1:numPoints
    %Convert velocity gradient into local coordinate system
    new_dV = [du_dx(j) du_dy(j) du_dz(j);...
              dv_dx(j) dv_dy(j) dv_dz(j);...
              dw_dx(j) dw_dy(j) dw_dz(j)] * R;
    new_du_dx = new_dV(1,1);
    new_du_dy = new_dV(1,2);
    new_du_dz = new_dV(1,3);
    new_dv_dx = new_dV(2,1);
    new_dv_dy = new_dV(2,2);
    new_dv_dz = new_dV(2,3);
    new_dw_dx = new_dV(3,1);
    new_dw_dy = new_dV(3,2);
    new_dw_dz = new_dV(3,3);
    
    %save velocity gradient info across body
    out.vel_grad(j) = new_du_dx;
    out.sh_def(j) = (2*(new_du_dx).^2 + 2*(new_dv_dy).^2 ...
        + 2*(new_dw_dz).^2 + (new_du_dy+new_dv_dx).^2 ...
        + (new_du_dz+new_dw_dx).^2 + (new_dv_dz+new_dw_dy).^2).^0.5;
end

% Flow speed
out.spd = sqrt(u.^2 + v.^2 + w.^2)';

% Body position values (arclength)
out.s = s;

% Inertial FOR body coorindates
out.bod_x = bod_x;
out.bod_y = bod_y;
out.bod_z = bod_z;

% 
% % Shear deformation (in inertial FOR)
% out.sh_def =  (2*(out(1).du_dx).^2 + 2*(out(1).dv_dy).^2 ...
%               + 2*(out(1).dw_dz).^2 + (out(1).du_dy+out(1).dv_dx).^2 ...
%               + (out(1).du_dz+out(1).dw_dx).^2 + (out(1).dv_dz+out(1).dw_dy).^2).^0.5;
          

