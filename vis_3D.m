function vis_3D(root_path)
% Makes 3D plots and animations of CFD and behavioral data


%% Code execution

% Interpolates CFD data with a regular mesh
do_cfd_field = 0;

% Creates 3D renderings 
do_3Dani = 1;

% Static 3D rendering of predator and all prey
do_3D = 0;


%% Parameters

% Period between detection and fast start initiation (s)
latency = 5e-3;

% Modeled threshold of shearing for a response
sh_thresh = .5;

% Modeled threshold of speed
spd_thresh = .1;

% Relative position of COM along length of the body (this is Matt's guess)
COM_pos = .25;

% Alpha transparency for the isosurface plots
bow_alpha = .3;

% Inerval for skipping plotting of normal vectors
norm_skip = 100;

% Scaling of body length for CFD volume examine for isosurface
sclfactr = .25;

% Scaling of body length for CFD volume around whole body
sclfactr_bod = 1.1;

% Number of bins used in rose plots
num_bin = 20;

% Define boundaries of CFD volume for analysis (cm)
rangeX = [-1 1];
rangeY = [-1.5 1.5];
rangeZ = [-1.5 1.5];

% Number of values along x-axis
numX = 100;

% Spacing between nodes 
dx = range(rangeX)/numX;

% Values along each dimension
xs = rangeX(1):dx:(rangeX(2)-dx);
ys = rangeY(1):dx:(rangeY(2)-dx);
zs = rangeZ(1):dx:(rangeZ(2)-dx);

% Number of points around the periphery of the prey body
numPts_circ = 50;

preyColor = [.5 .5 .5];

%% Paths

if nargin < 1
    root_path = uigetdir(pwd,'Select root directory (holds "cfd" & "behavior")');
end

% Paths to CFD data arranged in regular grid
cfd_path{1}  = [root_path filesep 'cfd' filesep 'flow_02cmps_biggrid.mat'];
cfd_path{2}  = [root_path filesep 'cfd' filesep 'flow_11cmps_biggrid.mat'];
cfd_path{3}  = [root_path filesep 'cfd' filesep 'flow_20cmps_biggrid.mat'];


% Paths to CFD data
raw_cfd_path{1}  = [root_path filesep 'cfd' filesep 'flow_02cmps_around_zebrafish.mat'];
raw_cfd_path{2}  = [root_path filesep 'cfd' filesep 'flow_11cmps_around_zebrafish.mat'];
raw_cfd_path{3}  = [root_path filesep 'cfd' filesep 'flow_20cmps_around_zebrafish.mat'];


%% Load data & define variables

% Load behavior data ('b')
load([root_path filesep 'behavior' filesep 'Transformed_Prey_Coords.mat'])

% Number of sequences
num_seq = length(b.preyx(:,1));

% Values of speeds
spds = [2 11 20];

% Indices for each speed of sequences in the dark, with lateral line intact
index{1} = (b.speed(1:num_seq)==2) & (b.LL(1:num_seq)==1) ...
     & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0);

index{2} = (b.speed(1:num_seq)==11) & (b.LL(1:num_seq)==1) ...
     & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0);
 
index{3} = (b.speed(1:num_seq)==20) & (b.LL(1:num_seq)==1) ...
     & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0);

% cBeate groups by speed
groups = [b.speed(index{1}); b.speed(index{2}); b.speed(index{3})];


%% Interpolate CFD data for each speed 

if do_cfd_field
    
for i = 1:3;
    
    % Update user
    disp(' ');disp(['Working on ' num2str(spds(i)) ' cm/s ...'])
    
    % cBeate x, y & z matrices
    [cB.x,cB.y,cB.z] = meshgrid(xs,ys,zs);
    
    % load cfd data ('c')
    load(raw_cfd_path{1})
    
    %Index of CFD node values within volumetric range
    idx = ((c.x >= min(rangeX)) & (c.x <= max(rangeX))) & ...
        ((c.y >= min(rangeY)) & (c.y <= max(rangeY))) & ...
        ((c.z >= min(rangeZ)) & (c.z <= max(rangeZ)));
    
    
    % Interpolate velocity components
    warning off
    u = griddata(c.x(idx),c.y(idx),c.z(idx),c.u(idx),cB.x,cB.y,cB.z);
    v = griddata(c.x(idx),c.y(idx),c.z(idx),c.v(idx),cB.x,cB.y,cB.z);
    w = griddata(c.x(idx),c.y(idx),c.z(idx),c.w(idx),cB.x,cB.y,cB.z);
   
    warning on
    
    % Flow speed
    cB.spd = sqrt(u.^2 + v.^2 + w.^2);

    % Save data
    save(cfd_path{i},'cB')
    
    clear cB idx
    
    disp('                                  ... Done!')
    
end

end


%% Plot 3D spatial distibution of responders

if do_3D
    
    do_pred = 1;
    
    do_bow = 1;
    
    
    p_clr = .5.*[1 1 1];
    
    % Designates speed index (e.g. 2 is 11 cm/s)
    i = 2;
    
    % Index of some interesting individuals for 11 cm/s sequences
    %prey = [1 3 8 12 15 24 27];
    
    % Load predator data (pred3Dshape... matricies)
    if do_pred
        load([root_path filesep 'morphology' filesep 'Pred3Dbodyshape.mat'])
    end
    
    % Load data, define dimensions ('m')
    load([root_path filesep 'morphology' filesep '6_02_L19_metrics.mat']);
    
    % Load CFD data in 'cB' strcuture
    if do_bow, load(cfd_path{i}); end
    
    % Define 3d data for prey in local FOR
    [pX,pY,pZ] = prey_surf(m.s,m.h,m.w,m.c,numPts_circ); clear m
    
    % Choose index for sequences that have stage 2 data
    idx = index{i} & ~isnan(b.preyx2(:,2)) & b.preyx(:,1)>0  ...
        & b.preyx(:,2)>0  & b.preyx(:,3)>0;
    
    % Get body points for all seqnences for current speed
    [rost0,com0,tail0,rost1,com1,tail1,rost2,com2,tail2,i_vent,...
        i_wrongLR,i_wrongDV,prey_dir] = give_points(b,idx,latency,spds(i));
    
    % New Figure window
    hF = figure;
    set(hF,'DoubleBuffer','on')
    
    % Render the predator -------------------------------------------------
    h1 = patch(real(pred3DshapeX),real(pred3DshapeY),...
        real(pred3DshapeZ),real(pred3DshapeX)*0);
    
    set(h1,'FaceLighting','gouraud',...
        'LineStyle','none',...
        'BackFaceLighting','reverselit',...
        'FaceColor',p_clr,...
        'AmbientStrength',.5);
    
    
    % Lighting, & view control -------------------------------------------
    axis tight
    
    lighting gouraud
    set(gca,'XColor','w','YColor','w','ZColor','w')
    
    %view([165 35])
    %view([180 0])

    view([140 40])

    % Render the prey -----------------------------------
    
    for j = 1:size(rost0,1)
        
        % Body length of prey
        blength = norm([tail0(j,1)-rost0(j,1) ...
            tail0(j,2)-rost0(j,2) ...
            tail0(j,3)-rost0(j,3)]);
        
        % Returns coordinates for prey body in global FOR
        [pXg,pYg,pZg] = prey_global(rost1(j,:),com1(j,:),tail1(j,:),...
            pX,pY,pZ);
        
        % Render the prey
        h1 = patch(real(pXg),real(pYg),real(pZg),real(pZg)*0);
        
        % Set properties
        set(h1,'FaceLighting','gouraud',...
            'LineStyle','none',...
            'BackFaceLighting','reverselit',...
            'FaceColor',preyColor,...
            'AmbientStrength',.5);
        hold on
        
        % Returns coordinates for prey body in global FOR
        [pXg2,pYg2,pZg2] = prey_global(rost2(j,:),com2(j,:),tail2(j,:),...
            pX,pY,pZ);
        
        % Render the prey at end of stage 2
        h2 = patch(real(pXg2),real(pYg2),real(pZg2),real(pZg2)*0);
        
        % Set properties
        set(h2,'FaceLighting','gouraud',...
            'LineStyle','none',...
            'BackFaceLighting','reverselit',...
            'FaceColor','r',...
            'AmbientStrength',.5);
        hold on
        
        
        axis equal
        
    end
    %         plot_prey_pos(f,index{i})
    hold on
    
    % Isosurface of bow wave
    if do_bow
        [faces,verts] = isosurface(cB.x,cB.y,cB.z,smooth3(cB.spd),spd_thresh);
        p = patch('Vertices', verts, 'Faces', faces, ...
            'FaceColor','interp', ...
            'edgecolor', 'interp');
        set(p,'FaceColor','b','EdgeColor','none');
        
        alpha(p,bow_alpha)
        
        % This is supposed to smooth the surface plot
        %isonormals(cB.x,cB.y,cB.z,cB.spd,p);
        
    end
    
    hold off
    
end


%% Animate experiment 3D spatial distibution of responders
% This is currently configured to focus on a few individual prey for 
% a 11 cm/s approach

if do_3Dani
    
    % Parameters ----------------------------------------------------------
    
    % Render the predator
    do_pred = 1;
    
    % Render the bow wave
    do_bow = 1;
    
    % frames per second playback speed
    play_rate = 30;
    
    % Color of prey body
    p_clr = .5.*[1 1 1];
    
    % Color of predator body
    pd_clr = .5.*[1 1 1];
    
    % Shadow color
    sh_clr = .9.*[1 1 1];
    
    % Designates speed index (e.g. 2 is 11 cm/s)
    i = 2;
    
    % Index of some interesting individuals for 11 cm/s sequences
    prey = [1 3 8 12 15 24 27];
    
    % Alpha tranparency of the bow wave
    bow_alpha = .4;
    
    % Starting position of predator (cm)
    start_pred = -2;
    
    % Spacing between prey (cm)
    prey_space = 1;
    
    % z-position of wall (for casting shadows)
    z_wall = -1.5;
    
    
    % Load data -----------------------------------------------------------
    
    % Load predator data (pred3Dshape... matricies)
    if do_pred
        load([root_path filesep 'morphology' filesep 'Pred3Dbodyshape.mat'])
    end
    
    % Load data, define dimensions ('m')
    load([root_path filesep 'morphology' filesep '6_02_L19_metrics.mat']);
    
    % Load CFD data in 'cB' strcuture
    if do_bow, load(cfd_path{i}); end
    
    % Define 3d data for prey in local FOR
    [pX,pY,pZ] = prey_surf(m.s,m.h,m.w,m.c,numPts_circ); clear m
    
    % Choose index for sequences that have stage 2 data
    idx = index{i} & ~isnan(b.preyx2(:,2)) & b.preyx(:,1)>0  ...
        & b.preyx(:,2)>0  & b.preyx(:,3)>0;
    

    % Get prey body points for all seqnences for current speed
    [rost0,com0,tail0,rost1,com1,tail1,rost2,com2,tail2,i_vent,...
        i_wrongLR,i_wrongDV,prey_dir] = give_points(b,idx,latency,spds(i));
    
    
    
    % Set kinematics for animation ----------------------------------------
    

    % Spacing of prey along x-axis
    prey_pos = linspace(0,prey_space*length(prey),length(prey));
    
    % Number of frames in animation
    n_frames = 2 * (abs(max(prey_pos)-start_pred) *  play_rate / spds(i));
        
    % Changes in position of predator's nose
    pos_vals = linspace(start_pred,max(prey_pos),n_frames);
    
    % Limits to the x-axis
    xlims = [-2 prey_space*length(prey)+.5];
    
    % Proximity between predtaor and prey
    prox = abs(repmat(pos_vals,length(prey_pos),1) - ...
               repmat(prey_pos',1,length(pos_vals)));
    
    
    
    % Render the prey -----------------------------------------------------
        
    
    % New Figure window
    hF = figure;
    set(hF,'DoubleBuffer','on','Color','w')
    
    % Render each prey
    for j = 1:length(prey)
        
        % Prey number
        n = prey(j);
        
        % Body length of prey
        blength = norm([tail0(n,1)-rost0(n,1) ...
                        tail0(n,2)-rost0(n,2) ...
                        tail0(n,3)-rost0(n,3)]);
                    
        % Tranformed position of prey (stage 1)
        rostT = [prey_pos(j)+rost1(n,1) rost1(n,2:3)];
        comT  = [prey_pos(j)+com1(n,1)  com1(n,2:3)];
        tailT = [prey_pos(j)+tail1(n,1) tail1(n,2:3)];
        
        % Tranformed position of prey (stage 2)
        rostT2 = [prey_pos(j)+rost2(n,1) rost2(n,2:3)];
        comT2  = [prey_pos(j)+com2(n,1)  com2(n,2:3)];
        tailT2 = [prey_pos(j)+tail2(n,1) tail2(n,2:3)];
        
        % Returns coordinates for prey body in global FOR
        [pXg,pYg,pZg] = prey_global(rostT,comT,tailT,pX,pY,pZ);
        
        % Render the prey
        h_prey1(j) = patch(real(pXg),real(pYg),real(pZg),real(pZg)*0);
        
        % Set properties
        set(h_prey1(j),'FaceLighting','gouraud',...
            'LineStyle','none',...
            'BackFaceLighting','reverselit',...
            'FaceColor',preyColor,...
            'AmbientStrength',.5);
        hold on
        
        % Shadow on x-y plane
        h_prey1S(j) = patch(pXg,pYg,0.*pZg+z_wall,pZg*0);
        set(h_prey1S(j),'LineStyle','none','FaceColor',sh_clr);
        hold on
        
        % Returns coordinates for prey body in global FOR
        [pXg2,pYg2,pZg2] = prey_global(rostT2,comT2,tailT2,pX,pY,pZ);
        
        % Render the prey
        h_prey2(j) = patch(real(pXg2),real(pYg2),real(pZg2),real(pZg2)*0);
        
        % Set properties
        set(h_prey2(j),'FaceLighting','gouraud',...
            'LineStyle','none',...
            'BackFaceLighting','reverselit',...
            'FaceColor','r',...
            'AmbientStrength',.5,...
            'Visible','off');
        
        % Shadow for stage 2 on x-y plane
        h_prey2S(j) = patch(pXg2,pYg2,0.*pZg2+z_wall,pZg2*0);
        set(h_prey2S(j),'LineStyle','none','FaceColor',sh_clr,...
                         'Visible','off');
        
        axis equal
        
        %pause(.1)
    end
    
    %hC = camlight;
    lighting gouraud
    set(gca,'XColor','w','YColor','w','ZColor','w')
    
    %view([165 35])
    %view([180 0])
    hL = light('position',[0 0 20]);
    
    view([115 22])
    
    j = 1;
        
    for k = 1:length(pos_vals)
        
        % Render the predator -------------------------------------------------
        h_pred = patch(pred3DshapeX+pos_vals(k),pred3DshapeY,...
                   pred3DshapeZ,pred3DshapeX*0);
        
        set(h_pred,'FaceLighting','gouraud',...
            'LineStyle','none',...
            'BackFaceLighting','reverselit',...
            'FaceColor',pd_clr,...
            'AmbientStrength',.5);

        % Shadow on x-y plane
        h_predS = patch(pred3DshapeX+pos_vals(k),pred3DshapeY,...
                        0.*pred3DshapeY+z_wall,pred3DshapeY.*0);
        set(h_predS,'LineStyle','none','FaceColor',sh_clr);
        
        
        % Isosurface
        if do_bow
            [faces,verts] = isosurface(cB.x+pos_vals(k),cB.y,cB.z,smooth3(cB.spd),spd_thresh);
            h_bow = patch('Vertices', verts, 'Faces', faces, ...
                'FaceColor','interp', 'edgecolor', 'interp');
            set(h_bow,'FaceColor','b','EdgeColor','none');
            alpha(h_bow,bow_alpha)
            
%             % Shadow for the bow wave causes it to crash:
%             h_bowS = patch(verts(:,1)+pos_vals(k),verts(:,2),...
%                            verts(:,3).*0+z_wall,verts(:,3).*0);
%             set(h_bowS,'LineStyle','none','FaceColor',sh_clr);
%             alpha(h_bowS,bow_alpha)
        end
        
        
        % If bow wave in contact with prey, switch to stage 2 larva
        if prox(j,k)==min(prox(j,:))
            
            set(h_prey1S(j),'Visible','off');
            set(h_prey1(j),'Visible','off');
            set(h_prey2S(j),'Visible','on');
            set(h_prey2(j),'Visible','on');
                        
            if j ~= length(prey)
                j = j + 1;
            end
            
        end
        
        hold off

        %view([113 -2])
        xlim([xlims(1) xlims(2)])
        %axis equal
        daspect([1,1,1])
        axis tight
        
        pause(.1)
        
        if k<length(pos_vals)
            delete(h_pred)
            delete(h_predS)
        end
        
        if do_bow, delete(h_bow);end
    end
end



function [pXg,pYg,pZg] = prey_global(rost,com,tail,pX,pY,pZ)
% Transforms prey points in global FOR

tmp1 = local_to_global(rost,com,tail,...
                        [pX(1,:)' pY(1,:)' pZ(1,:)']);
tmp2 = local_to_global(rost,com,tail,...
                        [pX(2,:)' pY(2,:)' pZ(2,:)']);
tmp3 = local_to_global(rost,com,tail,...
                        [pX(3,:)' pY(3,:)' pZ(3,:)']);
tmp4 = local_to_global(rost,com,tail,...
                        [pX(4,:)' pY(4,:)' pZ(4,:)']);

pXg = [tmp1(:,1)'; tmp2(:,1)'; tmp3(:,1)'; tmp4(:,1)'];
pYg = [tmp1(:,2)'; tmp2(:,2)'; tmp3(:,2)'; tmp4(:,2)'];
pZg = [tmp1(:,3)'; tmp2(:,3)'; tmp3(:,3)'; tmp4(:,3)'];


function  [rost0,com0,tail0,rost1,com1,tail1,rost2,com2,tail2,i_vent,...
           i_wrongLR,i_wrongDV,prey_dir] = give_points(b,idx,latency,spd)
% Returns the points of the prey body for all sequences denoted by 'idx'
       
% Offset in x, due to latency
lat_offset = latency * spd;

% Position (wrt predator) of body points when flow sensed
rost0 = [b.preyx(idx,1)+lat_offset b.preyy(idx,1) b.preyz(idx,1)];
com0  = [b.preyx(idx,2)+lat_offset b.preyy(idx,2) b.preyz(idx,2)];
tail0 = [b.preyx(idx,3)+lat_offset b.preyy(idx,3) b.preyz(idx,3)];

% Position (wrt predator) of body points when larva first moves
rost1 = [b.preyx(idx,1) b.preyy(idx,1) b.preyz(idx,1)];
com1  = [b.preyx(idx,2) b.preyy(idx,2) b.preyz(idx,2)];
tail1 = [b.preyx(idx,3) b.preyy(idx,3) b.preyz(idx,3)];

% Position (wrt predator) of body points at end of stage 2
rost2 = [b.preyx2(idx,1) b.preyy2(idx,1) b.preyz2(idx,1)];
com2  = [b.preyx2(idx,2) b.preyy2(idx,2) b.preyz2(idx,2)];
tail2 = [b.preyx2(idx,3) b.preyy2(idx,3) b.preyz2(idx,3)];

% Find direction of response
prey_dir(:,1) = com2(:,1) - com1(:,1);
prey_dir(:,2) = com2(:,2) - com1(:,2);
prey_dir(:,3) = com2(:,3) - com1(:,3);

% Transform prey direction, assuming mirror symmetry about the predator
prey_dir(com1(:,2)<0,2) = -prey_dir(com1(:,2)<0,2);

% Transform body points, assuming mirror symmetry about the predator
rost0(:,2)  = abs(rost0(:,2));
com0(:,2)   = abs(com0(:,2));
tail0(:,2)  = abs(tail0(:,2));

rost1(:,2)  = abs(rost1(:,2));
com1(:,2)   = abs(com1(:,2));
tail1(:,2)  = abs(tail1(:,2));

rost2(:,2)  = abs(rost2(:,2));
com2(:,2)   = abs(com2(:,2));
tail2(:,2)  = abs(tail2(:,2));

% Index of individuals positioned ventral to predator
i_vent = com0(:,3)<=0;

% Index of L-R responses in the 'wrong' direction
i_wrongLR = prey_dir(:,2)<0;

% Index of D-V responses in the 'wrong' direction
i_wrongDV = (i_vent & (prey_dir(:,3)>0)) |  ...
    (~i_vent & (prey_dir(:,3)<=0));



function ptsT = local_to_global(rost,com,tail,pts)
% Transforms coordinates (pts, in n x 3) from global to local coordinate
% system, assuming the y-axis of larvae runs perpendicular to gravity

% Check dimensions
if size(rost,1)~=1 || size(rost,2)~=3 || size(com,1)~=1 || size(com,2)~=3 ...
   size(tail,1)~=1 || size(tail,2)~=3 || size(pts,2)~=3 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = tail(1) - rost(1);
xaxis(1,2) = tail(2) - rost(2);
xaxis(1,3) = tail(3) - rost(3);

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis' yaxis' zaxis'];

% Rotate points
ptsT = [R * pts']';

% Translate global coordinates wrt rostrum
ptsT(:,1) = ptsT(:,1) + rost(1);
ptsT(:,2) = ptsT(:,2) + rost(2);
ptsT(:,3) = ptsT(:,3) + rost(3);

% Visualize to test
if 0
    
    blength = norm([tail(1)-rost(1) tail(2)-rost(2) tail(3)-rost(3)]);    
    
    figure
    
    subplot(2,2,[1 3])
    plot3([tail(1) rost(1)],[tail(2) rost(2)],[tail(3) rost(3)],'b',...
          rost(1),rost(2),rost(3),'bo');
    hold on
    plot3(ptsT(:,1),ptsT(:,2),ptsT(:,3),'ro')
    xlabel('x'); ylabel('y'); zlabel('z')
    hold off
    grid on;axis equal
    view(3)
    title('global')
    
    subplot(2,2,2)
    plot([0 blength],[0 0],'b',0,0,'ob',pts(:,1),pts(:,2),'ro')
    xlabel('x');ylabel('y')
    grid on; axis equal
    title('local')
    
    subplot(2,2,4)
    plot([0 blength],[0 0],'b',0,0,'ob',pts(:,1),pts(:,3),'ro')
    xlabel('x');ylabel('z')
    grid on; axis equal
end



function [X,Y,Z]= prey_surf(s,h,w,c,numPts)
% Provides 3d coordinates of the body

% Define radial positions along vector
theta = linspace(0,2*pi,numPts)';

% Define empty vectors for coordinates
x=[];y=[];z=[];

% Make mouth cap  _______________________________________
n = numPts/10;
phi = linspace(0,.75*pi/2,n)';
ds = .02.*range(s); %2*s(2)-s(1);
%sC = linspace(s(1)-ds,s(1),n);
hC = h(1) .* sin(phi)./max(sin(phi));
wC = w(1) .* sin(phi)./max(sin(phi));
sC = -(ds.*cos(phi)-ds.*cos(phi(end)));

% Loop down the body length
for i=1:length(sC)-1  
    
  % Draw first ellipse   
    yTemp1 = sC(i)*ones(size(theta));
    xTemp1 = (wC(i)/2) .* cos(theta);
    zTemp1 = (hC(i)/2) .* sin(theta) + c(1);
    
  % Draw second ellipse  
    yTemp2 = sC(i+1)*ones(size(theta));
    xTemp2 = (wC(i+1)/2) .* cos(theta);
    zTemp2 = (hC(i+1)/2) .* sin(theta) + c(1);
    
  % Combine data (works with 'patch')
    x	= [x [xTemp1(1:end-1)';... 
              xTemp2(1:end-1)';... 
              xTemp2(2:end)';... 
              xTemp1(2:end)']];
                      
    y   = [y [yTemp1(1:end-1)';... 
              yTemp2(1:end-1)';...
              yTemp2(2:end)';...
              yTemp1(2:end)']];
                      
    z   = [z [zTemp1(1:end-1)';...
              zTemp2(1:end-1)';...
              zTemp2(2:end)';...
              zTemp1(2:end)']];
end 

clear xTemp1 yTemp1 zTemp1 xTemp2 yTemp2 zTemp2
clear n phi ds sC hC wC


% Make body coordinates _______________________________________

% Loop down the body length
for i=1:length(s)-1  
    
  % Draw first ellipse  
    yTemp1      = s(i)*ones(size(theta));
    xTemp1      = (w(i)/2) .* cos(theta);
    zTemp1      = (h(i)/2) .* sin(theta) + c(i);
    
  % Draw second ellipse    
    yTemp2      = s(i+1)*ones(size(theta));
    xTemp2      = (w(i+1)/2) .* cos(theta);
    zTemp2      = (h(i+1)/2) .* sin(theta) + c(i+1);
    
  % Combine data (works with 'patch')
    x	= [x [xTemp1(1:end-1)';... 
              xTemp2(1:end-1)';... 
              xTemp2(2:end)';... 
              xTemp1(2:end)']];
                      
    y   = [y [yTemp1(1:end-1)';... 
              yTemp2(1:end-1)';...
              yTemp2(2:end)';...
              yTemp1(2:end)']];
                      
    z   = [z [zTemp1(1:end-1)';...
              zTemp2(1:end-1)';...
              zTemp2(2:end)';...
              zTemp1(2:end)']];
end  

% Make tail cap  _______________________________________
n = numPts/10;
phi = linspace(0,0.75*pi/2,n)';
ds = .02.*range(s); %2*s(2)-s(1);
%sC = linspace(s(1)-ds,s(1),n);
hC = h(end) .* sin(phi)./max(sin(phi));
wC = w(end) .* sin(phi)./max(sin(phi));
sC = s(end) + ds.*cos(phi)-+ ds.*cos(phi(end));

% Loop down the body length
for i=1:length(sC)-1  
    
  % Draw first ellipse   
    yTemp1 = sC(i)*ones(size(theta));
    xTemp1 = (wC(i)/2) .* cos(theta);
    zTemp1 = (hC(i)/2) .* sin(theta) + c(end);
    
  % Draw second ellipse  
    yTemp2 = sC(i+1)*ones(size(theta));
    xTemp2 = (wC(i+1)/2) .* cos(theta);
    zTemp2 = (hC(i+1)/2) .* sin(theta) + c(end);
    
  % Combine data (works with 'patch')
    x	= [x [xTemp1(1:end-1)';... 
              xTemp2(1:end-1)';... 
              xTemp2(2:end)';... 
              xTemp1(2:end)']];
                      
    y   = [y [yTemp1(1:end-1)';... 
              yTemp2(1:end-1)';...
              yTemp2(2:end)';...
              yTemp1(2:end)']];
                      
    z   = [z [zTemp1(1:end-1)';...
              zTemp2(1:end-1)';...
              zTemp2(2:end)';...
              zTemp1(2:end)']];
end 

clear xTemp1 yTemp1 zTemp1 xTemp2 yTemp2 zTemp2
clear n phi ds sC hC wC

% Transform for units and wrt rostrum
X  = (y-min(y(:))).*100;
Y  = x.*100;
Z  = -(z-z(1)+max(h(:))/2).*100;