function d = piv_compile(d_path)
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


%% Load & organize data


% Loop through speeds
for i = 1:length(spds)
    
    % String that specifies speed
    str_spd = ['0' num2str(spds(i))];
    str_spd = str_spd(end-1:end);
    
    % Initialize vectors to store data
    x = []; y = []; z = [];
    u = []; v = []; w = [];
    
    % Load VENTRAL data at all positions ----------------------------------
    for j = 1:length(cal.V.position)
        
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
        
        % Translate coordinates along z-axis wrt global origin
        coord_gbl(:,3) = coord_gbl(:,3) + slice_pos;
        
        % Define coordinates in the fish FOR
        x   = [x; coord_gbl(:,1)];
        y   = [y; coord_gbl(:,2)];
        z   = [z; coord_gbl(:,3)];
        
        % Define velocity components in fish FOR
        u   = [u; vel_gbl(:,1)];
        v   = [v; vel_gbl(:,2)];
        w   = [w; nan(length(vel_gbl),1)];   
        
        clear conv_ratio predx predy predz str_depth f_name coord_gbl 
        clear coord_gbl vel_gbl coord_vid vel_vid C slice_pos
    end
    
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
        
        % Define coordinates in the fish FOR
        x   = [x; coord_gbl(:,1)];
        y   = [y; coord_gbl(:,2)];
        z   = [z; coord_gbl(:,3)];
        
        % Define velocity components in fish FOR
        u   = [u; vel_gbl(:,1)];
        v   = [v; vel_gbl(:,2)];
        w   = [w; nan(length(vel_gbl),1)];   
        
        clear conv_ratio predx predy predz str_depth f_name coord_gbl 
        clear coord_gbl vel_gbl coord_vid vel_vid C slice_pos
    end
        
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
        
        % Translate coordinates along z-axis wrt global origin
        coord_gbl(:,2) = coord_gbl(:,2) + slice_pos;
        
        % Define coordinates in the fish FOR
        x   = [x; coord_gbl(:,1)];
        y   = [y; coord_gbl(:,2)];
        z   = [z; coord_gbl(:,3)];
        
        % Define velocity components in fish FOR
        u   = [u; vel_gbl(:,1)];
        v   = [v; nan(length(vel_gbl),1)];
        w   = [w; vel_gbl(:,3)];   
        
        clear conv_ratio predx predy predz str_depth f_name coord_gbl 
        clear coord_gbl vel_gbl coord_vid vel_vid C slice_pos
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
        
        % Translate coordinates along z-axis wrt global origin
        coord_gbl(:,2) = coord_gbl(:,2) + slice_pos;
        
        % Define coordinates in the fish FOR
        x   = [x; coord_gbl(:,1)];
        y   = [y; coord_gbl(:,2)];
        z   = [z; coord_gbl(:,3)];
        
        % Define velocity components in fish FOR
        u   = [u; vel_gbl(:,1)];
        v   = [v; nan(length(vel_gbl),1)];
        w   = [w; vel_gbl(:,3)];   
        
        clear conv_ratio predx predy predz str_depth f_name coord_gbl 
        clear coord_gbl vel_gbl coord_vid vel_vid C slice_pos
    end
     
    % Store in 'p' structure
    p.x = x;
    p.y = y;
    p.z = z;
    p.u = u;
    p.v = v;
    p.w = w;
    
    % Save data
    save([d_path filesep 'compiled data' filesep 'S' str_spd '_piv'],'p')
    
    % Clear for next iteration
    clear x_vid y_vid z u_vid v w str_spd
    
end


