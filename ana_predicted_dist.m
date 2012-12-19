function ana_predicted_dist(root_path,mode)
% Analyzes the flow along the body of larvae from body position data in 
% behavioral experiments and compares the response distance to what is
% predicted for the measured orientation of the larva


%% Parameters

% Period between detection and fast start initiation (s)
latency = 5e-3;

% Number of points along the body
numPoints = 100;

% Hypothetical threshold speed (cm/s)
thresh_spd = .1;

% Hypothetical threshold shear 
thresh_shear = .3;

% Hypothetical threshold speed gradient (1/s)
thresh_velgrad = .2;

% Step size in body position (cm)
b_pos_step = 0.01;


%% Paths

if nargin < 1
    root_path = uigetdir(pwd,'Select root directory (holds "cfd" & "behavior")');
end

% Load latest behavior coordinates ('b' structure)
load([root_path filesep 'behavior' filesep 'Transformed_Prey_Coords.mat'])

% Paths to CFD data
cfd_path{1}   = [root_path filesep 'cfd' filesep 'flow_02cmps_around_zebrafish.mat'];
cfd_path{2}  = [root_path filesep 'cfd' filesep 'flow_11cmps_around_zebrafish.mat'];
cfd_path{3}  = [root_path filesep 'cfd' filesep 'flow_20cmps_around_zebrafish.mat'];

% Mode can bel 'vel', 'shear' or 'velgrad'
if nargin < 2
mode = 'velgrad';
end

%% Set up variables

% Store latency assumed in the analysis
r.latency = latency;

% Measured and predicted response distances
r.resp_pos_meas = nan(size(b.preyx,1),1);
r.resp_pos_pred = nan(size(b.preyx,1),1);

% Number of sequences
num_seq = length(b.preyx(:,1));

% Indices for 2 cm/s
idx{1} = find((b.speed(1:num_seq)==2) & (b.LL(1:num_seq)==1) ...
           & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0));

% Indices for 11 cm/s       
idx{2} = find((b.speed(1:num_seq)==11) & (b.LL(1:num_seq)==1) ...
           & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0));

% Indices for 20 cm/s       
idx{3} = find((b.speed(1:num_seq)==20) & (b.LL(1:num_seq)==1) ...
           & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0));

% Speed of predator's approach
pred_spd{1} = 2;
pred_spd{2} = 11;
pred_spd{3} = 20;


%% Interpolate for all sequences

if strcmp(mode,'vel') && ...
        ~isempty(dir([root_path filesep 'behavior' filesep 'vel_thresh_test.mat']))
    
    runn = questdlg('Rerun vel threshold test?','','No','Yes',...
                      'Cancel','No');
                  
elseif strcmp(mode,'shear') && ...
        ~isempty(dir([root_path filesep 'behavior' filesep 'shear_thresh_test.mat']))
    
    runn = questdlg('Rerun shear threshold test?','','No','Yes',...
                      'Cancel','No');
                  
elseif strcmp(mode,'velgrad') && ...
        ~isempty(dir([root_path filesep 'behavior' filesep 'velgrad_thresh_test.mat']))
    
    runn = questdlg('Rerun velgrad threshold test?','','No','Yes',...
                      'Cancel','No');    
else
    runn = 'Yes';
end


if strcmp(runn,'Cancel')
    return

elseif strcmp(runn,'Yes');
    
% Loop through speeds
for i = 1:3;
    
    % Import 2 cm/s cfd data ('c' structure)
    load(cfd_path{i});
    
    % Offset in x, due to latency & pred speed
    lat_offset = latency * pred_spd{i};
    
    % Loop through sequences for present speed
    for j = 1:length(idx{i})
        
        % Sequence number
        seqnum = idx{i}(j);
        
        % Get properties of the prey body
        [blength,R,body] = prey_props(b,seqnum,numPoints,lat_offset);
        
        % Distance of each body point to predator
        dist = sqrt( body(:,1).^2 + body(:,2).^2 + body(:,3).^2 );
        
        % Index for body position closest to the predator
        i_bod_sense = find(dist == min(dist),1,'first');
        
        % Response distance
        resp_dist      = dist(i_bod_sense);
        r.resp_pos_meas(seqnum)  = body(i_bod_sense,1);
        
        % Starting position of that part of the body
        b_pos = min([0 min(body(i_bod_sense,1))]);
        
        while true
            
            % Translate body position to bod_pos
            body(:,1) = body(:,1) - body(i_bod_sense,1) + b_pos;
            
            % Find flow conditions
            flow = interpflow(c,body,blength,R);
            
            % Check that position isn't excessively far
            if b_pos > blength*10;
                warning(['Predicted response distance not found for ' ...
                    ' seq # ' num2str(seqnum)])
                r.resp_pos_pred(seqnum) = nan;
                break
            end
            
            % Check if above threshold
            if strcmp(mode,'vel') && (max(flow.spd) < thresh_spd)
                r.resp_pos_pred(seqnum) = b_pos;
                break
                
            elseif strcmp(mode,'shear') && (max(flow.sh_def) < thresh_shear)
                r.resp_pos_pred(seqnum) = b_pos;
                break
             
            elseif strcmp(mode,'velgrad') && (max(flow.sh_def) < thresh_velgrad)
                r.resp_pos_pred(seqnum) = b_pos;
                break    
            end
               
            b_pos = b_pos + b_pos_step;
        end
        
        %r.resp_pos_meas(seqnum),r.resp_pos_pred(seqnum)
        
        clear flow x_vals y_vals z_vals seqnum
    end
    
    clear c
    
    disp(' ');disp(['Done ' num2str(pred_spd{i}) ' cm/s'])
    
end

if strcmp(mode,'vel')
    save([root_path filesep 'behavior' filesep 'vel_thresh_test'],'r');
    
elseif strcmp(mode,'shear')
    save([root_path filesep 'behavior' filesep 'shear_thresh_test'],'r');
    
elseif strcmp(mode,'velgrad')
    save([root_path filesep 'behavior' filesep 'velgrad_thresh_test'],'r');
end    

else
    % Load 'r' structure
    if strcmp(mode,'vel')
        load([root_path filesep 'behavior' filesep 'vel_thresh_test'])
    elseif strcmp(mode,'velgrad')
        load([root_path filesep 'behavior' filesep 'velgrad_thresh_test'])
    elseif strcmp(mode,'shear')
        load([root_path filesep 'behavior' filesep 'shear_thresh_test'])
    end

end



%% Plot results


figure;
plot(r.resp_pos_meas(idx{1}),r.resp_pos_pred(idx{1}),'ro',...
     r.resp_pos_meas(idx{2}),r.resp_pos_pred(idx{2}),'bo',...
     r.resp_pos_meas(idx{3}),r.resp_pos_pred(idx{3}),'go',...
     [min(r.resp_pos_meas) max(r.resp_pos_meas)],...
     [min(r.resp_pos_meas) max(r.resp_pos_meas)],'--k')
xlabel('x-position measured (cm)'); 
ylabel('x-position predicted')
axis equal
title(mode)





function [blength,R,body] = prey_props(b,seqnum,numPoints,lat_offset)
% Returns properties of the prey's body
% bod - position of body coordinate when larva responds

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
body(:,1) = bod(:,1) + b_origin(1) + lat_offset;
body(:,2) = bod(:,2) + b_origin(2);
body(:,3) = bod(:,3) + b_origin(3);



function out = interpflow(c,body,blength,R)
% Interpolates CFD flow field to find flow conditions at input coordinates
% c - structure of CFD data
% b - structure of 

% Scaling of body length for CFD volume to interpolate over
sclfactr = 1.5;

% Mean position of body
meanX = mean(body(:,1));
meanY = mean(body(:,2));
meanZ = mean(body(:,3));

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
u = griddata(c.x(idx),c.y(idx),c.z(idx),c.u(idx),...
             body(:,1),body(:,2),body(:,3));
v = griddata(c.x(idx),c.y(idx),c.z(idx),c.v(idx),...
             body(:,1),body(:,2),body(:,3));
w = griddata(c.x(idx),c.y(idx),c.z(idx),c.w(idx),...
             body(:,1),body(:,2),body(:,3));

% Interpolate velocity gradients
du_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dx(idx),...
                 body(:,1),body(:,2),body(:,3));
du_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dy(idx),...
                 body(:,1),body(:,2),body(:,3));
du_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.du_dz(idx),...
                 body(:,1),body(:,2),body(:,3));
dv_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dx(idx),...
                 body(:,1),body(:,2),body(:,3));
dv_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dy(idx),...
                 body(:,1),body(:,2),body(:,3));
dv_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.dv_dz(idx),...
                 body(:,1),body(:,2),body(:,3));
dw_dx = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dx(idx),...
                 body(:,1),body(:,2),body(:,3));
dw_dy = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dy(idx),...
                 body(:,1),body(:,2),body(:,3));
dw_dz = griddata(c.x(idx),c.y(idx),c.z(idx),c.dw_dz(idx),...
                 body(:,1),body(:,2),body(:,3));

% Turn warnings back on
warning on

% Calculate flow in local FOR
for j = 1:size(body,1)
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
    out.vel_grad(j) = abs(new_du_dx);
    out.sh_def(j) = (2*(new_du_dx).^2 + 2*(new_dv_dy).^2 ...
        + 2*(new_dw_dz).^2 + (new_du_dy+new_dv_dx).^2 ...
        + (new_du_dz+new_dw_dx).^2 + (new_dv_dz+new_dw_dy).^2).^0.5;
end

% Flow speed
out.spd = sqrt(u.^2 + v.^2 + w.^2)';



% 
% % Shear deformation (in inertial FOR)
% out.sh_def =  (2*(out(1).du_dx).^2 + 2*(out(1).dv_dy).^2 ...
%               + 2*(out(1).dw_dz).^2 + (out(1).du_dy+out(1).dv_dx).^2 ...
%               + (out(1).du_dz+out(1).dw_dx).^2 + (out(1).dv_dz+out(1).dw_dy).^2).^0.5;
          

