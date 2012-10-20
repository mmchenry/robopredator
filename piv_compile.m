function piv_compile(d_path)
% Imports piv data for each speed and dumps it into individual data files

%% Parameters

% Predator speeds at which experiments were conducted
%spds = [2 4 6 8 11 12 14 16 18 20];

spds = [2 11 20];

% Listing of perpectives of camera
perspective   = {'D','V','L','R'};

% Pixel dimensions of video frame (assumed square)
fr_size = 1024;

% Size in pixels of header at top of video frames
yOffset = 48;

% Bounds and number of points for uniform grid
x_num = 60;
x_min = -0.6;
x_max = 1.6;
y_min = -.4;
y_max = 1.8;


%% Check directories

% Browse to directory, if not given
if nargin < 1
    d_path = uigetdir(pwd,'Select root directory');
end

% Check that subdirectories are present
c_path = [d_path filesep 'derived data'];
if isempty(dir(c_path))
    error([c_path ' is missing']);
end

c_path = [d_path filesep 'calibration data'];
if isempty(dir(c_path))
    error([c_path ' is missing']);
end

c_path = [d_path filesep 'compiled data'];
if isempty(dir(c_path))
    error([c_path ' is missing']);
end

clear c_path


%% Compile calibration data into 'cal' structure

a = dir([d_path filesep 'calibration data' filesep 'calibration_data*.mat']);

iV = 1; iD = 1; iL = 1; iR = 1;
for i = 1:length(a)
    
    % Define current perspective and plane
    c_per = a(i).name(end-4);
    c_pos = str2num((a(i).name((end-7):(end-6))));
    
    % Load data (conv_ratio, xprednosePix, yprednosePix)
    load([d_path filesep 'calibration data' filesep a(i).name])
    
    % Store data into appropriate field
    if strcmp(c_per,'V')
        cal.V.position(iV)      = c_pos;
        cal.V.conv_ratio(iV)    = conv_ratio;
        cal.V.xprednosePix(iV)  = xprednosePix;
        cal.V.yprednosePix(iV)  = yprednosePix;
        iV = iV + 1;
        
    elseif strcmp(c_per,'D')
        cal.D.position(iD)      = c_pos;
        cal.D.conv_ratio(iD)    = conv_ratio;
        cal.D.xprednosePix(iD)  = xprednosePix;
        cal.D.yprednosePix(iD)  = yprednosePix;
        iD = iD + 1;
        
    elseif strcmp(c_per,'L')
        cal.L.position(iL)      = c_pos;
        cal.L.conv_ratio(iL)    = conv_ratio;
        cal.L.xprednosePix(iL)  = xprednosePix;
        cal.L.yprednosePix(iL)  = yprednosePix;
        iL = iL + 1;   
        
    elseif strcmp(c_per,'R')
        cal.R.position(iR)      = c_pos;
        cal.R.conv_ratio(iR)    = conv_ratio;
        cal.R.xprednosePix(iR)  = xprednosePix;
        cal.R.yprednosePix(iR)  = yprednosePix;
        iR = iR + 1;
        
    else
        error('perspective code not recognized')
    end
        
    clear c_pos c_per conv_ratio xprednosePix yprednosePix
    
end

% Save data
save([d_path filesep 'compiled data' filesep 'calibrations'],'cal')

% Clear variables
clear iV iD iL iR



%% Define coordinate systems

% DORSAL perspective ------------------------------------------------------
xAxis   = [-1 0 0];
yAxis   = [0 -1 0];
zAxis   = [0 0 1];
 
% Transformation matrix
S_D = [xAxis' yAxis' zAxis'];

clear xAxis yAxis zAxis

% VENTRAL perspective -----------------------------------------------------
% Note: the field viewed from this perspective should
% reside mostly in negative y-values.  However, we are fliping around the
% y-axis so that the ventral and dorsal fields are assumed to record flow
% wrt the left side of the predator's body.

S_V = S_D;

% RIGHT perspective -------------------------------------------------------
xAxis   = [-1 0 0];
yAxis   = [0 0 -1];
zAxis   = [0 -1 0];
 
% Transformation matrix
S_R = [xAxis' yAxis' zAxis'];

clear xAxis yAxis zAxis

% LEFT perspective -------------------------------------------------------
xAxis   = [-1 0 0];
yAxis   = [0 0 1];
zAxis   = [0 1 0];
 
% Transformation matrix
S_L = [xAxis' yAxis' zAxis'];

clear xAxis yAxis zAxis


%% Organize dorso-ventral (DV) data


% Define common coordinate values wrt predator with cm & cm/s units
xVals = linspace(x_min,x_max,x_num)';
dVals = mean(diff(xVals));
yVals = [y_min:dVals:y_max]';

clear dVals x_min y_min x_max y_max

% Loop through speeds
for i = 1:length(spds)
    
    % String that specifies speed
    str_spd = ['0' num2str(spds(i))];
    str_spd = str_spd(end-1:end);
    
    % Initialize index for z values
    z_idx = 1;
 
    % Load VENTRAL data at all positions ----------------------------------
    for j = length(cal.V.position):-1:1
        
        % Calibration constant
        conv_ratio = cal.V.conv_ratio(j);
        
        % x position of origin (flipped in l/r direction, in cm)
        predx = ( cal.V.xprednosePix(j));
        
        % y_vid position of origin (flipped in u_vid/d direction, header removed, in cm)
        predy = (fr_size -(cal.V.yprednosePix(j) - yOffset));
        
        % z position of origin (plane position: cm -> mm)
        slice_pos = cal.V.position(j)/10;
        
        % String that specifies depth
        str_depth = ['0' num2str(cal.V.position(j))];
        str_depth = str_depth(end-1:end);
        
        % Current data filename
        f_name = ['timeAvg_H' str_depth '_V_S' str_spd '.mat'];
        
        % Load 'C' matrix (derived piv data)
        load([d_path filesep 'derived data' filesep f_name])
        
        % Define coordinates and velocity wrt orgin in video FOR
        coord_vid  = [conv_ratio.*(C(:,1) - predx) ...
                      conv_ratio.*(C(:,2) - predy) ...
                      C(:,1).*0];
        vel_vid    = [C(:,3) C(:,4) C(:,3).*0];   
            
        % Rotate coordinates into global FOR
        coord_gbl = [inv(S_V)'*coord_vid']';
        vel_gbl   = [inv(S_V)'*vel_vid']'; 
        
        % Define current transformed coordinates
        x_c = coord_gbl(:,1);
        y_c = coord_gbl(:,2);
        u_c = vel_gbl(:,1);
        v_c = vel_gbl(:,2);
        
        % Meshgrid and interpolate the transformed data
        [x_c,y_c,u_c,v_c] = mesh_interp(x_c,y_c,u_c,v_c,xVals,yVals);

        % Define coordinates in the fish FOR
        x(:,:,z_idx)   = x_c;
        y(:,:,z_idx)   = y_c;
        z(z_idx)       = slice_pos;
        
        % Define velocity components in fish FOR
        u(:,:,z_idx)   = u_c;
        v(:,:,z_idx)   = v_c; 
        
        z_idx = z_idx + 1;
 
        clear conv_ratio predx predy predz str_depth f_name coord_gbl 
        clear coord_gbl vel_gbl coord_vid vel_vid C slice_pos
        clear x_c y_c z_c u_c v_c w_c
    end
    
    clear j
    
    
    % Load DORSAL data at all positions -----------------------------------
    for j = 1:length(cal.D.position)
        
        % Calibration constant
        conv_ratio = cal.D.conv_ratio(j);
        
        % x position of origin (flipped in l/r direction, in cm)
        predx = ( cal.D.xprednosePix(j));
        
        % y position of origin (flipped in u/d direction, header removed, in cm)
        predy =  (fr_size -(cal.D.yprednosePix(j) - yOffset));
        
        % z position of origin (plane position: cm -> mm)
        slice_pos = cal.D.position(j)/10;
        
        % String that specifies depth
        str_depth = ['0' num2str(cal.D.position(j))];
        str_depth = str_depth(end-1:end);
        
        % Current data filename
        f_name = ['timeAvg_H' str_depth '_D_S' str_spd '.mat'];
        
        % Load 'C' matrix (derived piv data)
        load([d_path filesep 'derived data' filesep f_name])
 
        % Define coordinates and velocity wrt orgin in video FOR
        coord_vid  = [conv_ratio.*(C(:,1) - predx) ...
                      conv_ratio.*(C(:,2) - predy) ...
                      C(:,1).*0];
        vel_vid    = [C(:,3) C(:,4) C(:,3).*0];   
            
        % Rotate coordinates into global FOR
        coord_gbl = [inv(S_D)'*coord_vid']';
        vel_gbl   = [inv(S_D)'*vel_vid']';
        
        % Translate coordinates along z-axis wrt global origin
        coord_gbl(:,3) = coord_gbl(:,3) + slice_pos;
        
        % Define current transformed coordinates
        x_c = coord_gbl(:,1);
        y_c = coord_gbl(:,2);
        u_c = vel_gbl(:,1);
        v_c = vel_gbl(:,2);
        
        % Meshgrid and interpolate the transformed data
        [x_c,y_c,u_c,v_c] = mesh_interp(x_c,y_c,u_c,v_c,xVals,yVals);

        % Define coordinates in the fish FOR
        x(:,:,z_idx)   = x_c;
        y(:,:,z_idx)   = y_c;
        z(z_idx)   = slice_pos;
        
        % Define velocity components in fish FOR
        u(:,:,z_idx)   = u_c;
        v(:,:,z_idx)   = v_c;
        
        z_idx = z_idx + 1; 
        
        clear conv_ratio predx predy predz str_depth f_name coord_gbl 
        clear coord_gbl vel_gbl coord_vid vel_vid C slice_pos
        clear x_c y_c z_c u_c v_c w_c
    end
    
    % Store values. The x and y coords should be the same for all z.
    p.spd = spds(i);
    p.x = x(:,:,1);
    p.y = y(:,:,1);
    p.z = z;
    p.u = u;
    p.v = v;
    p.w = nan;
    
    % Save data
    save([d_path filesep 'compiled data' filesep 'DV view S' str_spd '_piv'],'p')
    
    clear X_tmp Y_tmp x y z u v z_idx str_spd j p
    
end

clear i

%% Organize lateral (LR) data

% Loop through speeds
for i = 1:length(spds)
    
    % String that specifies speed
    str_spd = ['0' num2str(spds(i))];
    str_spd = str_spd(end-1:end);
        
    % Initialize index for y values
    y_idx = 1;
    
    % Load RIGHT data at all positions ------------------------------------
    for j = 1:length(cal.R.position)
        
        % Calibration constant
        conv_ratio = cal.R.conv_ratio(j);
        
        % x position of origin (flipped in l/r direction, in cm)
        predx = ( cal.R.xprednosePix(j));
        
        % y position of origin (in cm)
        slice_pos = -cal.R.position(j)/10;
        
        % z position of origin (plane position: cm -> mm)
        predy = (fr_size -(cal.R.yprednosePix(j) - yOffset));
        
        % String that specifies depth
        str_depth = ['0' num2str(cal.R.position(j))];
        str_depth = str_depth(end-1:end);
        
        % Current data filename
        f_name = ['timeAvg_H' str_depth '_R_S' str_spd '.mat'];
        
        % Load 'C' matrix (derived piv data)
        load([d_path filesep 'derived data' filesep f_name])
        
        % Define coordinates and velocity wrt orgin in video FOR
        coord_vid  = [conv_ratio.*(C(:,1) - predx) ...
                      conv_ratio.*(C(:,2) - predy) ...
                      C(:,1).*0];
        vel_vid    = [C(:,3) C(:,4) C(:,3).*0];   
            
        % Rotate coordinates into global FOR
        coord_gbl = [inv(S_R)'*coord_vid']';
        vel_gbl   = [inv(S_R)'*vel_vid']';
        
        % Define current transformed coordinates
        x_c = coord_gbl(:,1);
        z_c = coord_gbl(:,3);
        u_c = vel_gbl(:,1);
        w_c = vel_gbl(:,3);
        
        % Meshgrid and interpolate the transformed data
        [x_c,z_c,u_c,w_c] = mesh_interp(x_c,z_c,u_c,w_c,xVals,yVals);

        % Define coordinates in the fish FOR
        x(:,y_idx,:)   = x_c;
        y(y_idx)       = slice_pos;
        z(:,y_idx,:)   = z_c;
        
        % Define velocity components in fish FOR
        u(:,y_idx,:)   = u_c;
        w(:,y_idx,:)   = w_c;
        
        y_idx = y_idx + 1; 
        
        clear conv_ratio predx predy predz str_depth f_name coord_gbl 
        clear coord_gbl vel_gbl coord_vid vel_vid C slice_pos
        clear x_c y_c z_c u_c v_c w_c
    end
    
    % Load LEFT data at all positions -------------------------------------
    for j = 1:length(cal.L.position)
        
        % Calibration constant
        conv_ratio = cal.L.conv_ratio(j);
        
        % x position of origin (flipped in l/r direction, in cm)
        predx = (cal.L.xprednosePix(j));
        
        % y position of origin (flipped in u/d direction, header removed, in cm)
        slice_pos = cal.L.position(j)/10;
        
        % z position of origin (plane position: cm -> mm)
        predy = (fr_size -(cal.L.yprednosePix(j) - yOffset));
        
        % String that specifies depth
        str_depth = ['0' num2str(cal.L.position(j))];
        str_depth = str_depth(end-1:end);
        
        % Current data filename
        f_name = ['timeAvg_H' str_depth '_L_S' str_spd '.mat'];
        
        % Load 'C' matrix (derived piv data)
        load([d_path filesep 'derived data' filesep f_name])
        
        % Define coordinates and velocity wrt orgin in video FOR
        coord_vid  = [conv_ratio.*(C(:,1) - predx) ...
                      conv_ratio.*(C(:,2) - predy) ...
                      C(:,1).*0];
        vel_vid    = [C(:,3) C(:,4) C(:,3).*0];   
            
        % Rotate coordinates into global FOR
        coord_gbl = [inv(S_L)'*coord_vid']';
        vel_gbl   = [inv(S_L)'*vel_vid']';
        
        % Define current transformed coordinates
        x_c = coord_gbl(:,1);
        z_c = coord_gbl(:,3);
        u_c = vel_gbl(:,1);
        w_c = vel_gbl(:,3);
        
        % Meshgrid and interpolate the transformed data
        [x_c,z_c,u_c,w_c] = mesh_interp(x_c,z_c,u_c,w_c,xVals,yVals);

        % Define coordinates in the fish FOR
        x(:,y_idx,:)   = x_c;
        y(y_idx)       = slice_pos;
        z(:,y_idx,:)   = z_c;
        
        % Define velocity components in fish FOR
        u(:,y_idx,:)   = u_c;
        w(:,y_idx,:)   = w_c;
        
        y_idx = y_idx + 1; 
        
        clear conv_ratio predx predy predz str_depth f_name coord_gbl 
        clear coord_gbl vel_gbl coord_vid vel_vid C slice_pos
        clear x_c y_c z_c u_c v_c w_c
    end
     
    % Store values. The x and z coords should be the same for all y.
    p.spd = spds(i);
    p.x = x(:,1,:);
    p.z = z(:,1,:);
    p.y = y;
    p.u = u;
    p.v = nan;
    p.w = w;
    
    % Save data
    save([d_path filesep 'compiled data' filesep 'LR view S' str_spd '_piv'],'p')
    
    clear X_tmp Y_tmp x y z u v z_idx str_spd j p
    
end


function [X,Y,U,V] = mesh_interp(x,y,u,v,X_vals,Y_vals)
% Interpolates and meshgrids the piv data

% Indicies for the start of each row of x values, This assumes that no mask
% occurs in on the left-hand side of the video frame
i_rowstart = find(x == x(1));

% Check that a row isn't getting added to another row, due to a left-hand
% side mask
if max(diff(i_rowstart)) > 1.3*mean(diff(i_rowstart))
    warning(['This code assumes that no mask occurs on the left-hand ' ...
             'side of the video frame.  Check that assumption'])
end

% Number of unique x and y values that repeat in a grid
num_xvals = length(unique(x));
num_yvals = length(i_rowstart);

% New unique x & y values
%X_vals = linspace(min(x),max(x),num_xvals*.75);
%Y_vals = linspace(min(y),max(y),num_yvals*.75)';

% Step through y-v
for i = 1:num_yvals
    
    % Indicies for current row
    if i==num_yvals
        i_c = i_rowstart(i):length(x);
    else
        i_c = i_rowstart(i):(i_rowstart(i+1)-1);
    end
    
    % Current values
    x_c = x(i_c);
    y_c = y(i_c);
    u_c = interp1(x(i_c),u(i_c),X_vals);
    v_c = interp1(x(i_c),v(i_c),X_vals);
    
    % Ditch nans (for interp1, below)
    u_c(isnan(u_c)) = 0;
    v_c(isnan(v_c)) = 0;
    
    % Interpolate y, u & v and store in current row
    X1(i,1:length(X_vals)) = X_vals;
    Y1(i,1:length(X_vals)) = X_vals.*0 + mean(y_c(~isnan(y_c)));
    U1(i,1:length(X_vals)) = u_c;
    V1(i,1:length(X_vals)) = v_c;
    
    clear i_c x_c y_c u_c v_c
end

clear i


% Next, create even spacing along the y-values by interpolating along the 
% y dimension at each x-coordinate
for i = 1:size(X1,2)
    
    % Define new x and y coordiates
    Y(:,i) = Y_vals;
    X(:,i) = Y_vals.*0 + X1(1,i);
    
    % Define current u and v values
    u_c = interp1(Y1(:,i),U1(:,i),Y_vals);
    v_c = interp1(Y1(:,i),V1(:,i),Y_vals);
    
    % Replace nans with zeros
    u_c(u_c==0) = nan;
    v_c(v_c==0) = nan;
    
    % Store results
    U(:,i) = u_c;
    V(:,i) = v_c;
    
    clear u_c v_c
    
end


