function run_sims(root_path)
% Run simulations of prey responses, given the CFD of the bow wave of
% predators

% Calculate flow around larvae for simulations (takes a long time)
do_flow = 1;

% Make 3D plot of results (note: use simulations with spare sampling)
vis_3D = 0;

% Test the effects of response wrt fixed angle of the body
do_algo_fixed = 0;

% Test the effects of a variety of algorithms
do_algo_rand = 0;

% Test algo where the larva moves in direction of velocity
do_algo_vel = 0;

% Index for approach speed to analyze
iSpeed = 2;

rose_rand = 1;

rose_fixed = 1;


%% Paths

if  nargin < 1
    root_path = uigetdir(pwd,'Select root directory (holds "cfd" & "behavior")');
    
    if isempty(root_path)
        return
    end    
end

% Paths to CFD data arranged in regular grid (for bow wave isosurface)
cfd_path{1}  = [root_path filesep 'cfd' filesep 'flow_02cmps_biggrid.mat'];
cfd_path{2}  = [root_path filesep 'cfd' filesep 'flow_11cmps_biggrid.mat'];
cfd_path{3}  = [root_path filesep 'cfd' filesep 'flow_20cmps_biggrid.mat'];

% Paths to CFD data arranged in regular grid (for local flow)
cfd_pathR{1}  = [root_path filesep 'cfd' filesep 'flow_02cmps_reggrid.mat'];
cfd_pathR{2}  = [root_path filesep 'cfd' filesep 'flow_11cmps_reggrid.mat'];
cfd_pathR{3}  = [root_path filesep 'cfd' filesep 'flow_20cmps_reggrid.mat'];

% Paths for flow simulation data
simflow_path{1} = [root_path filesep 'simulations' filesep 'flow data 2'];
simflow_path{2} = [root_path filesep 'simulations' filesep 'flow data 11'];
simflow_path{3} = [root_path filesep 'simulations' filesep 'flow data 20'];

% Paths for algorithm simulation data (fixed angle)
simalgo_fixed_path{1} = [root_path filesep 'simulations' filesep 'algofixed data 2'];
simalgo_fixed_path{2} = [root_path filesep 'simulations' filesep 'algofixed data 11'];
simalgo_fixed_path{3} = [root_path filesep 'simulations' filesep 'algofixed data 20'];

% Paths for algorithm simulation data (random angle)
simalgo_rand_path{1} = [root_path filesep 'simulations' filesep 'algorand data 2'];
simalgo_rand_path{2} = [root_path filesep 'simulations' filesep 'algorand data 11'];
simalgo_rand_path{3} = [root_path filesep 'simulations' filesep 'algorand data 20'];


%% Parameters

% Speeds
spds = [2 11 20];

% Number of start angles
num_ang = 50;

% Number of prey positions to consider along y-axis 
% (determines spatial sampling density for simulations)
num_y = 20;

% Modeled threshold of speed
spd_thresh = .07;

% Transparency of bow wave
bow_alpha = 0.5;

% Predator surface color
p_clr = .5.*[1 1 1];

% Start angles (rad)
%start_ang = linspace(-pi,pi,36)';
%start_ang = start_ang(1:end-1);

% Position of body points in local coordinates
rostL = [0 0 0];
comL  = [0.07 0 0];
tailL = [0.40 0 0];

% Scale factor for determining subvolume (body lengths)
sclfactr = 1.5;

% Body length of larva (cm)
blength = tailL(1);

% Displacement of fast start distance (cm)
FS_dist = 0.17;

% Series of points to evaluate flow along the body (body lengths)
s = linspace(0,1,20)';

% Fixed angle of FS relative to body (rad)
ang_fix = 100/180*pi;

% Lateral position at which a prey is categorized as 'lateral' (cm)
% (Also used for dorsal and ventral positions)
lat_pos = 0.5;

num_bin = 20;


%% Calculate flow data


if do_flow

tic
    
% Load predator morphology
load([root_path filesep 'morphology' filesep 'Pred3Dbodyshape.mat'])


% Points along the L and R of the body
tmp = ones(length(s),1);
l_ptL = [s.*blength tmp.*(-blength/10) 0.*tmp];
r_ptL = [s.*blength tmp.*(blength/10) 0.*tmp];
clear tmp

% Index for saving data
n = 1;

% Loop through speeds
for i = 1:3
    
    start_ang = 2*pi*(rand(num_ang,1)-.5);
    
    % Find coordinates of larvae at bow wave ------------------------------
    
    % Range of volume to interrogate
    rangeX   = [0 2];
    rangeY   = [-1.5 1.5];
    rangeZ   = [-1.5 1.5];
    subRange = [rangeX(1) rangeX(2) ...
                rangeY(1) rangeY(2) ...
                rangeZ(1) rangeZ(2)];
    
    % Load CFD data in 'cB' strcuture
    load(cfd_path{i});
    
    % Reduce flow field to small volume
    [xS,yS,zS,spdS] = subvolume(cB.x,cB.y,cB.z,cB.spd,subRange);
    
    % Find isosurface of bow wave
    [faces,verts] = isosurface(xS,yS,zS,smooth3(spdS),spd_thresh);
    
    % Y & Z positions to consider
    Ytmp = linspace(min(verts(:,2)),max(verts(:,2)),num_y)';
    Ztmp = [min(verts(:,3)):mean(diff(Ytmp)):max(verts(:,3))]';
    [Yvals,Zvals] = meshgrid(Ytmp,Ztmp);
    
    % Find corresponding X value on bow wave surface
    warning off
    Xvals = griddata(verts(:,2),verts(:,3),verts(:,1),Yvals,Zvals);
    warning on
    
    % Make into vector
    Xvals = Xvals(:);
    Yvals = Yvals(:);
    Zvals = Zvals(:);
    
    % Remove nan points
    i_pts = ~isnan(Xvals);
    Xvals = Xvals(i_pts);
    Yvals = Yvals(i_pts);
    Zvals = Zvals(i_pts);
    
    clear Ytmp Ztmp rangeX rangeY rangeZ subRange xS yS zS spdS 
    
    % Visualize results----------------------------------------------------
    if 0
        figure
            % Render the predator 
            h1 = patch(real(pred3DshapeX),real(pred3DshapeY),...
                real(pred3DshapeZ),real(pred3DshapeX)*0);
            
            % Set properties
            set(h1,'FaceLighting','gouraud',...
                'LineStyle','none',...
                'BackFaceLighting','reverselit',...
                'FaceColor',p_clr,...
                'AmbientStrength',.5);
            hold on
        
        % Render bow wave
        h_bow = patch('Vertices', verts, 'Faces', faces, ...
            'FaceColor','interp', 'edgecolor', 'interp');
        set(h_bow,'FaceColor','b','EdgeColor','none');
        alpha(h_bow,bow_alpha)
        
        % Plot positions of virtual prey
        hold on
        plot3(Xvals,Yvals,Zvals,'ok')
        
        % Adjust appearence
        axis equal       
        axis tight
        lighting gouraud
    end 
    
    
    % Find flow conditions at each position -------------------------------
     
    % Load CFD data in 'cR' structure
    load(cfd_pathR{i})
    
    % Loop through each larva position
    for j = 1:length(Xvals)

        % Range of volume to interrogate
        rangeX   = [Xvals(j)-blength*sclfactr ...
                    Xvals(j)+blength*sclfactr];
        rangeY   = [Yvals(j)-blength*sclfactr ...
                    Yvals(j)+blength*sclfactr];
        rangeZ   = [Zvals(j)-blength*sclfactr/2 ...
                    Zvals(j)+blength*sclfactr/2];
        subRange = [rangeX(1) rangeX(2) ...
                    rangeY(1) rangeY(2) ...
                    rangeZ(1) rangeZ(2)];
        
        % Reduce flow field to small volume
        [xS,yS,zS,spdS] = subvolume(cR.x,cR.y,cR.z,cR.spd,subRange);
        [xS,yS,zS,uS]   = subvolume(cR.x,cR.y,cR.z,cR.u,subRange);
        [xS,yS,zS,vS]   = subvolume(cR.x,cR.y,cR.z,cR.v,subRange);
        [xS,yS,zS,wS]   = subvolume(cR.x,cR.y,cR.z,cR.w,subRange);
            
        % Loop through each angular position
        for k = 1:length(start_ang)
            
            % Rotate body points according to current angle
            ptsL(1,:) = [0 0 0];
            ptsL(2,:) = [comL(1).*cos(start_ang(k)) ...
                         comL(1).*sin(start_ang(k)) 0];
            ptsL(3,:) = [tailL(1).*cos(start_ang(k)) ...
                         tailL(1).*sin(start_ang(k)) 0];
            
            % Translate body points to current position
            ptsG(:,1) = ptsL(:,1) + Xvals(j);
            ptsG(:,2) = ptsL(:,2) + Yvals(j);
            ptsG(:,3) = ptsL(:,3) + Zvals(j);
                      
            % L and R points along body transformed into global cooridnate system
            l_pt = local_to_global(ptsG(1,:),ptsG(2,:),ptsG(3,:),l_ptL);
            r_pt = local_to_global(ptsG(1,:),ptsG(2,:),ptsG(3,:),r_ptL);
  
            % Package points for griddata
            pnts = [l_pt; r_pt];
            
            % Interpolate CFD data for speed/velocity 
            spd_vals = griddata(xS,yS,zS,spdS,pnts(:,1),pnts(:,2),pnts(:,3));
            u_vals   = griddata(xS,yS,zS,uS,pnts(:,1),pnts(:,2),pnts(:,3));
            v_vals   = griddata(xS,yS,zS,vS,pnts(:,1),pnts(:,2),pnts(:,3));
            w_vals   = griddata(xS,yS,zS,wS,pnts(:,1),pnts(:,2),pnts(:,3));
            
            % Check output
            if max(isnan(spd_vals))
                error('nan values generated by griddata -- enlarge subRange');
            end

            % Speed values at different spots on the body
            spd_left  = spd_vals(1:length(l_pt));
            spd_right = spd_vals((length(l_pt)+1):(length(l_pt)+length(r_pt)));  
            u_left    = u_vals(1:length(l_pt));
            u_right   = u_vals((length(l_pt)+1):(length(l_pt)+length(r_pt)));  
            v_left    = v_vals(1:length(l_pt));
            v_right   = v_vals((length(l_pt)+1):(length(l_pt)+length(r_pt)));  
            w_left    = w_vals(1:length(l_pt));
            w_right   = w_vals((length(l_pt)+1):(length(l_pt)+length(r_pt)));  
           
            % Store data           
            sim.spd(n,1)    = spds(i);
            sim.rostG(n,:)  = [Xvals(j) Yvals(j) Zvals(j)];
            sim.comG(n,:)   = ptsG(2,:);
            sim.tailG(n,:)  = ptsG(3,:); 
            sim.ang(n,1)    = start_ang(k);
            sim.s{n}        = s;
            sim.spd_r{n}    = spd_right;
            sim.spd_l{n}    = spd_left;
            sim.u_r{n}      = u_right;
            sim.v_r{n}      = v_right;
            sim.w_r{n}      = w_right;
            sim.u_l{n}      = u_left;
            sim.v_l{n}      = v_left;
            sim.w_l{n}      = w_left;
            
            % Increment index
            n = n + 1;
            
            % Visualize to check
            if 0
                plot3([ptsG(1,1) ptsG(3,1)],[ptsG(1,2) ptsG(3,2)],...
                      [ptsG(1,3) ptsG(3,3)],'-',...
                      ptsG(1,1),ptsG(1,2),ptsG(1,3),'o');
            end
            
            clear spd_vals pnts spd_left spd_right pnts l_pt r_pt ptsL ptsG
        end
         
        clear rangeX rangeY rangeZ subRange
        
        % Update
        disp(['         Done ' num2str(j) ' of ' ...
              num2str(length(Xvals)) ' positions']) 
        
    end
    
    % Update user
    disp(' '); disp(['Done ' num2str(spds(i)) ' cm/s']); disp(' ');beep
    
    disp('Saving results')
    save([root_path filesep 'simulations' filesep 'flow data ' ...
          num2str(spds(i))],'sim')
end

toc

% else
%     disp(' '); disp('Loading flow data . . .')
%     
%     % load flow data ('sim')
%     %load([root_path filesep 'simulations' filesep 'flow data 2']);
%     %sim2 = sim;
%     
%     load([root_path filesep 'simulations' filesep 'flow data 11']);
%     %sim11 = sim;
%     
%     %load([root_path filesep 'simulations' filesep 'flow data 20']);
%     %sim20 = sim;
%     
%     %clear sim
end


%% Simulate response for different behavioral algorithms (fixed angle)

if do_algo_fixed

    % Loop through speeds
    for i = 1:3
        
        % Load flow data 'sim'
        load(simflow_path{i})
        
        algo.away_fix = [];
        algo.tow_fix = [];
        algo.rand_fix = [];
        
        % Loop thru simulations
        for j = 1:length(sim.spd)
            
            % Index of body position with max flow
            iBod = find(sim.spd_r{j} == max(sim.spd_r{j}),1,'first');
            
            % Adjust sign of algorithm, according to side with higher flow
            if sim.spd_r{j}(iBod) > sim.spd_l{j}(iBod)
                ang_away = - (pi - ang_fix);
                ang_tow  = (pi - ang_fix);
                
            elseif max(isnan(sim.spd_r{j}(iBod))) || max(isnan(sim.spd_r{j}(iBod)))
                error('There are nans in the flow data');
                
            else
                ang_away = (pi - ang_fix);
                ang_tow  = - (pi - ang_fix);   
            end
            
            % Random selection of toward or away
            if rand(1)>.5
                ang_rand = ang_tow;
            else
                ang_rand = ang_away;
            end            
            
            % Algorithm: Side away, fixed angle ------------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_away) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_away) + comL(2);
            pos2L(1,3) = 0;
            
            % Store results
            algo.away_fix = give_algo_res(sim,algo.away_fix,pos2L,ang_away,lat_pos,j);
            
            clear pos2L
            
            % Algorithm: Side toward, fixed angle ----------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_tow) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_tow) + comL(2);
            pos2L(1,3) = 0;
            
            algo.tow_fix = give_algo_res(sim,algo.tow_fix,pos2L,ang_tow,lat_pos,j);
            
            clear pos2L
            
            % Algorithm: Side rand, fixed angle ----------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_rand) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_rand) + comL(2);
            pos2L(1,3) = 0;
            
            algo.rand_fix = give_algo_res(sim,algo.rand_fix,pos2L,ang_rand,lat_pos,j);
            clear pos2L
            
            % Visualize to check calculations
            if 0

               subplot(2,3,1)
               plot([comL(1) FS_dist*cos(ang_away) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_away) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('side away, fixed, local')
               axis equal
               
               subplot(2,3,2)
               plot([comL(1) FS_dist*cos(ang_tow) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_tow) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('size toward, fixed, local')
               axis equal
               
               subplot(2,3,3)
               plot([comL(1) FS_dist*cos(ang_rand) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_rand) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('side rand, fixed, local')
               axis equal
               
               subplot(2,3,4)
               plot([algo.away_fix.rost1(end,1) algo.away_fix.tail1(end,1)],...
                    [algo.away_fix.rost1(end,2) algo.away_fix.tail1(end,2)],'r-',...
                    algo.away_fix.rost1(end,1),algo.away_fix.rost1(end,2),'or',...
                    [algo.away_fix.com1(end,1) algo.away_fix.com2(end,1)],...
                    [algo.away_fix.com1(end,2) algo.away_fix.com2(end,2)],'k-')
                title('side away, fixed, global')
               axis equal
               
               subplot(2,3,5)
               plot([algo.tow_fix.rost1(end,1) algo.tow_fix.tail1(end,1)],...
                    [algo.tow_fix.rost1(end,2) algo.tow_fix.tail1(end,2)],'r-',...
                    algo.tow_fix.rost1(end,1),algo.tow_fix.rost1(end,2),'or',...
                    [algo.tow_fix.com1(end,1) algo.tow_fix.com2(end,1)],...
                    [algo.tow_fix.com1(end,2) algo.tow_fix.com2(end,2)],'k-')
                title('size toward, fixed, global')
               axis equal
               
               subplot(2,3,6)
               plot([algo.rand_fix.rost1(end,1) algo.rand_fix.tail1(end,1)],...
                    [algo.rand_fix.rost1(end,2) algo.rand_fix.tail1(end,2)],'r-',...
                     algo.rand_fix.rost1(end,1),algo.rand_fix.rost1(end,2),'or',...
                    [algo.rand_fix.com1(end,1) algo.rand_fix.com2(end,1)],...
                    [algo.rand_fix.com1(end,2) algo.rand_fix.com2(end,2)],'k-')
                title('side rand, fixed, global')
               axis equal
            end       
          end
        
        % Save
        disp(['Saving results for ' num2str(spds(i)) ' cm/s'])
        save(simalgo_fixed_path{i},'algo')
        
        clear sim algo
    end      
end


%% Simulate response for different behavioral algorithms (random angle)

if do_algo_rand

    % Loop through speeds
    for i = 1:3
        
        % Load flow data 'sim'
        load(simflow_path{i})
        
        algo.away_rand = [];
        algo.tow_rand = [];
        algo.rand_rand = [];
        
        % Loop thru simulations
        for j = 1:length(sim.spd)
            
            % Random angle
            ang_rand = rand(1)*pi;
            
            % Index of body position with max flow
            iBod = find(sim.spd_r{j} == max(sim.spd_r{j}),1,'first');
            
            % Adjust sign of algorithm, according to side with higher flow
            if sim.spd_r{j}(iBod) > sim.spd_l{j}(iBod)
                ang_away = - (pi - ang_rand);
                ang_tow  = (pi - ang_rand);
                
            elseif max(isnan(sim.spd_r{j}(iBod))) || max(isnan(sim.spd_r{j}(iBod)))
                error('There are nans in the flow data');
                
            else
                ang_away = (pi - ang_rand);
                ang_tow  = - (pi - ang_rand);   
            end
            
            % Random selection of toward or away
            if rand(1)>.5
                ang_rand = ang_tow;
            else
                ang_rand = ang_away;
            end            
            
            % Algorithm: Side away, fixed angle ------------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_away) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_away) + comL(2);
            pos2L(1,3) = 0;
            
            % Store results
            algo.away_rand = give_algo_res(sim,algo.away_rand,pos2L,ang_away,lat_pos,j);
            
            clear pos2L
            
            % Algorithm: Side toward, fixed angle ----------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_tow) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_tow) + comL(2);
            pos2L(1,3) = 0;
            
            algo.tow_rand = give_algo_res(sim,algo.tow_rand,pos2L,ang_tow,lat_pos,j);
            
            clear pos2L
            
            % Algorithm: Side rand, fixed angle ----------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_rand) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_rand) + comL(2);
            pos2L(1,3) = 0;
            
            algo.rand_rand = give_algo_res(sim,algo.rand_rand,pos2L,ang_rand,lat_pos,j);
            clear pos2L
            
            % Visualize to check calculations
            if 0

               subplot(2,3,1)
               plot([comL(1) FS_dist*cos(ang_away) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_away) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('side away, fixed, local')
               axis equal
               
               subplot(2,3,2)
               plot([comL(1) FS_dist*cos(ang_tow) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_tow) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('size toward, fixed, local')
               axis equal
               
               subplot(2,3,3)
               plot([comL(1) FS_dist*cos(ang_rand) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_rand) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('side rand, fixed, local')
               axis equal
               
               subplot(2,3,4)
               plot([algo.away_rand.rost1(end,1) algo.away_rand.tail1(end,1)],...
                    [algo.away_rand.rost1(end,2) algo.away_rand.tail1(end,2)],'r-',...
                    algo.away_rand.rost1(end,1),algo.away_rand.rost1(end,2),'or',...
                    [algo.away_rand.com1(end,1) algo.away_rand.com2(end,1)],...
                    [algo.away_rand.com1(end,2) algo.away_rand.com2(end,2)],'k-')
                title('side away, fixed, global')
               axis equal
               
               subplot(2,3,5)
               plot([algo.tow_rand.rost1(end,1) algo.tow_rand.tail1(end,1)],...
                    [algo.tow_rand.rost1(end,2) algo.tow_rand.tail1(end,2)],'r-',...
                    algo.tow_rand.rost1(end,1),algo.tow_rand.rost1(end,2),'or',...
                    [algo.tow_rand.com1(end,1) algo.tow_rand.com2(end,1)],...
                    [algo.tow_rand.com1(end,2) algo.tow_rand.com2(end,2)],'k-')
                title('size toward, fixed, global')
               axis equal
               
               subplot(2,3,6)
               plot([algo.rand_rand.rost1(end,1) algo.rand_rand.tail1(end,1)],...
                    [algo.rand_rand.rost1(end,2) algo.rand_rand.tail1(end,2)],'r-',...
                     algo.rand_rand.rost1(end,1),algo.rand_rand.rost1(end,2),'or',...
                    [algo.rand_rand.com1(end,1) algo.rand_rand.com2(end,1)],...
                    [algo.rand_rand.com1(end,2) algo.rand_rand.com2(end,2)],'k-')
                title('side rand, fixed, global')
               axis equal
            end       
          end
        
        % Save
        disp(['Saving results for ' num2str(spds(i)) ' cm/s'])
        save(simalgo_rand_path{i},'algo')
        
        clear sim algo
    end         
end


%% Simulate response for different behavioral algorithms (velocity)

if do_algo_vel

    % Loop through speeds
    for i = 1:3
        
        % Load flow data 'sim'
        load(simflow_path{i})
        
        algo.away_rand = [];
        algo.tow_rand = [];
        algo.rand_rand = [];
        
        % Loop thru simulations
        for j = 1:length(sim.spd)
            
            % Random angle
            ang_rand = rand(1)*pi;
            
            % Index of body position with max flow
            iBod = find(sim.spd_r{j} == max(sim.spd_r{j}),1,'first');
            
            % Adjust sign of algorithm, according to side with higher flow
            if sim.spd_r{j}(iBod) > sim.spd_l{j}(iBod)
                ang_away = - (pi - ang_rand);
                ang_tow  = (pi - ang_rand);
                
            elseif max(isnan(sim.spd_r{j}(iBod))) || max(isnan(sim.spd_r{j}(iBod)))
                error('There are nans in the flow data');
                
            else
                ang_away = (pi - ang_rand);
                ang_tow  = - (pi - ang_rand);   
            end
            
            % Random selection of toward or away
            if rand(1)>.5
                ang_rand = ang_tow;
            else
                ang_rand = ang_away;
            end            
            
            % Algorithm: Side away, fixed angle ------------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_away) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_away) + comL(2);
            pos2L(1,3) = 0;
            
            % Store results
            algo.away_rand = give_algo_res(sim,algo.away_rand,pos2L,ang_away,lat_pos,j);
            
            clear pos2L
            
            % Algorithm: Side toward, fixed angle ----------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_tow) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_tow) + comL(2);
            pos2L(1,3) = 0;
            
            algo.tow_rand = give_algo_res(sim,algo.tow_rand,pos2L,ang_tow,lat_pos,j);
            
            clear pos2L
            
            % Algorithm: Side rand, fixed angle ----------------------------------
            
            % Displacement in local system
            pos2L(1,1) =  FS_dist*cos(ang_rand) + comL(1);
            pos2L(1,2) =  FS_dist*sin(ang_rand) + comL(2);
            pos2L(1,3) = 0;
            
            algo.rand_rand = give_algo_res(sim,algo.rand_rand,pos2L,ang_rand,lat_pos,j);
            clear pos2L
            
            % Visualize to check calculations
            if 0

               subplot(2,3,1)
               plot([comL(1) FS_dist*cos(ang_away) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_away) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('side away, fixed, local')
               axis equal
               
               subplot(2,3,2)
               plot([comL(1) FS_dist*cos(ang_tow) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_tow) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('size toward, fixed, local')
               axis equal
               
               subplot(2,3,3)
               plot([comL(1) FS_dist*cos(ang_rand) + comL(1)],...
                    [comL(2) FS_dist*sin(ang_rand) + comL(2)],'k-',...
                    [0 4e-1],[0 0],'r-',0,0,'or')
                title('side rand, fixed, local')
               axis equal
               
               subplot(2,3,4)
               plot([algo.away_rand.rost1(end,1) algo.away_rand.tail1(end,1)],...
                    [algo.away_rand.rost1(end,2) algo.away_rand.tail1(end,2)],'r-',...
                    algo.away_rand.rost1(end,1),algo.away_rand.rost1(end,2),'or',...
                    [algo.away_rand.com1(end,1) algo.away_rand.com2(end,1)],...
                    [algo.away_rand.com1(end,2) algo.away_rand.com2(end,2)],'k-')
                title('side away, fixed, global')
               axis equal
               
               subplot(2,3,5)
               plot([algo.tow_rand.rost1(end,1) algo.tow_rand.tail1(end,1)],...
                    [algo.tow_rand.rost1(end,2) algo.tow_rand.tail1(end,2)],'r-',...
                    algo.tow_rand.rost1(end,1),algo.tow_rand.rost1(end,2),'or',...
                    [algo.tow_rand.com1(end,1) algo.tow_rand.com2(end,1)],...
                    [algo.tow_rand.com1(end,2) algo.tow_rand.com2(end,2)],'k-')
                title('size toward, fixed, global')
               axis equal
               
               subplot(2,3,6)
               plot([algo.rand_rand.rost1(end,1) algo.rand_rand.tail1(end,1)],...
                    [algo.rand_rand.rost1(end,2) algo.rand_rand.tail1(end,2)],'r-',...
                     algo.rand_rand.rost1(end,1),algo.rand_rand.rost1(end,2),'or',...
                    [algo.rand_rand.com1(end,1) algo.rand_rand.com2(end,1)],...
                    [algo.rand_rand.com1(end,2) algo.rand_rand.com2(end,2)],'k-')
                title('side rand, fixed, global')
               axis equal
            end       
          end
        
        % Save
        disp(['Saving results for ' num2str(spds(i)) ' cm/s'])
        save(simalgo_rand_path{i},'algo')
        
        clear sim algo
    end         
end


%% Visualize results

if vis_3D

    % Load predator morphology
    load([root_path filesep 'morphology' filesep 'Pred3Dbodyshape.mat'])
    
    % Load CFD data in 'cB' strcuture
    load(cfd_path{iSpeed});
    
    % Range of volume to interrogate
    rangeX   = [0 2];
    rangeY   = [-1.5 1.5];
    rangeZ   = [-1.5 1.5];
    subRange = [rangeX(1) rangeX(2) ...
                rangeY(1) rangeY(2) ...
                rangeZ(1) rangeZ(2)];

    % Reduce flow field to small volume
    [xS,yS,zS,spdS] = subvolume(cB.x,cB.y,cB.z,cB.spd,subRange);
    
    % Loop through each algorithm
    
    curr = sim.away_fix;
    t_txt = 'Away-Fix';
    
    render_results(xS,yS,zS,spdS,spd_thresh,pred3DshapeX,pred3DshapeY,...
                        pred3DshapeZ,curr,t_txt,p_clr)

    clear cB pred3DshapeX pred3DshapeY pred3DshapeZ curr rangeX rangeY rangeZ 
    clear xS yS zS spdS t_test
end


%% Analysis: Percentage correct (fixed)

for i = 1:3
    % Load 'algo'
    load(simalgo_fixed_path{i})
    
    disp(' ');disp(['Percentage away from predator ' num2str(spds(i)) ...
        ' (cm/s) -----------------']);disp(' ')
    
    disp(['   Away side, fixed angle: ' num2str(100*...
        sum(~algo.away_fix.wrong)/length(algo.away_fix.wrong)) '% '])
    disp(' ')
    disp(['   Toward side, fixed angle: ' num2str(100*...
        sum(~algo.tow_fix.wrong)/length(algo.tow_fix.wrong)) '% '])
    disp(' ')
    disp(['   Random side, fixed angle: ' num2str(100*...
        sum(~algo.rand_fix.wrong)/length(algo.rand_fix.wrong)) '% '])
    disp(' '); disp(' ')
    
    clear algo
end


%% Analysis: Percentage correct (rand)

for i = 1:3
    % Load 'algo'
    load(simalgo_rand_path{i})
    
    disp(' ');disp(['Percentage away from predator ' num2str(spds(i)) ...
        ' (cm/s) -----------------']);disp(' ')
    
    disp(['   Away side, rand angle: ' num2str(100*...
        sum(~algo.away_rand.wrong)/length(algo.away_rand.wrong)) '% '])
    disp(' ')
    disp(['   Toward side, rand angle: ' num2str(100*...
        sum(~algo.tow_rand.wrong)/length(algo.tow_rand.wrong)) '% '])
    disp(' ')
    disp(['   Random side, rand angle: ' num2str(100*...
        sum(~algo.rand_rand.wrong)/length(algo.rand_rand.wrong)) '% '])
    disp(' '); disp(' ')
    
    clear algo
end


%% Analysis: Global rose plots (fixed ang)

if rose_fixed
% Loop thru speeds
for i = 2 %1:3
    
    % Load 'algo'
    load(simalgo_fixed_path{i})

    % Algorithm: Side away, fixed angle ------------------------------------

    rosey(algo.away_fix,num_bin,lat_pos,'Away-fixed',spds(i))
    
     % Algorithm: Side away, fixed angle ------------------------------------

    rosey(algo.tow_fix,num_bin,lat_pos,'Toward-fixed',spds(i))
    
     % Algorithm: Side away, fixed angle ------------------------------------

    rosey(algo.rand_fix,num_bin,lat_pos,'Rand-fixed',spds(i))
    
end 

end


%% Analysis: Global rose plots (rand ang)

if rose_rand
% Loop thru speeds
for i = 2 %1:3
    
    % Load 'algo'
    load(simalgo_rand_path{i})

    % Algorithm: Side away, fixed angle ------------------------------------

    rosey(algo.away_rand,num_bin,lat_pos,'Away-rand',spds(i))
    
     % Algorithm: Side away, fixed angle ------------------------------------

    rosey(algo.tow_rand,num_bin,lat_pos,'Toward-rand',spds(i))
    
     % Algorithm: Side away, fixed angle ------------------------------------

    rosey(algo.rand_rand,num_bin,lat_pos,'Rand-rand',spds(i))
    
end   
    
end


%% Functions ----------------------------------------------------------------
    
function render_results(xS,yS,zS,spdS,spd_thresh,pred3DshapeX,pred3DshapeY,...
                        pred3DshapeZ,c,t_txt,p_clr)
 bow_alpha=.5;
% Find isosurface of bow wave
[faces,verts] = isosurface(xS,yS,zS,smooth3(spdS),spd_thresh);

figure

% Render the predator
h1 = patch(real(pred3DshapeX),real(pred3DshapeY),...
    real(pred3DshapeZ),real(pred3DshapeX)*0);

% Set properties
set(h1,'FaceLighting','gouraud',...
    'LineStyle','none',...
    'BackFaceLighting','reverselit',...
    'FaceColor',p_clr,...
    'AmbientStrength',.5);
hold on

% % Render bow wave
% h_bow = patch('Vertices', verts, 'Faces', faces, ...
%     'FaceColor','interp', 'edgecolor', 'interp');
% set(h_bow,'FaceColor','b','EdgeColor','none');
% alpha(h_bow,bow_alpha)

% Adjust appearence

%axis tight
lighting gouraud

% Plot positions of virtual prey
hold on
for i = 1:size(c.com1,1)
    h1 = plot3(c.rost1(i,1),c.rost1(i,2),c.rost1(i,3),'o');
    hold on
    h2 = plot3([c.rost1(i,1) c.tail1(i,1)],[c.rost1(i,2) c.tail1(i,2)],...
               [c.rost1(i,3) c.tail1(i,3)],'-k');
    h3 = plot3([c.com1(i,1) c.com2(i,1)],[c.com1(i,2) c.com2(i,2)],...
               [c.com1(i,3) c.com2(i,3)],'-m');
    set(h1,'Color',.7.*[1 1 1]);
    set(h2,'Color',.7.*[1 1 1]);
    if c.wrong(i)
        set(h3,'Color','r')
    else
        set(h3,'Color','g')
    end
end

xl = xlim;

hL = plot3([0 max(xlim)],[0 0],[0 0],'k-');

axis equal
view([68 18])
title(t_txt)

function res = give_algo_res(sim,res,pos2L,angl,lat_pos,i)
% Returns the results of a simulation, given particular angular response
   
% COM for times 1 & 2
com1 = sim.comG(i,:);
com2 = local_to_global(sim.rostG(i,:),sim.comG(i,:),sim.tailG(i,:),pos2L); 

% Calculate dist from midline at t1 and t2
dist1 = sqrt(com1(2)^2 + com1(3)^2);
dist2 = sqrt(com2(2).^2 + com2(3).^2);
  
% Azimuth of response
[az,rad] = cart2pol(com2(1)-com1(1),com2(2)-com1(2));

% Fast start displacement
disp = sqrt((com2(1)-com1(1))^2 + (com2(2)-com1(2))^2 + (com2(3)-com1(3))^2 ) ;

% Record results
res.rost1(i,:) = sim.rostG(i,:);
res.tail1(i,:) = sim.tailG(i,:);
res.com1(i,:)  = com1;
res.com2(i,:)  = com2;
res.wrong(i,1) = dist2 < dist1;
res.azG(i,1)   = az;
res.azL(i,1)   = angl;
res.vent(i,1)  = com1(3)<0;
res.lat(i,1)   = abs(com1(2))>lat_pos;


function rosey(curr,num_bin,lat_pos,t_txt,spds)

figure

% Lateral --------------------------------------------------
idx = (abs(curr.com1(:,2))>lat_pos);

% Define current azimuth
c_az   = curr.azG(idx);
c_com1 = curr.com1(idx,:);

% Flip sign of azimuth, if negative y coordinate
c_az(c_com1(:,2)<0) = -c_az(c_com1(:,2)<0);

subplot(2,1,1)
rose_plot(c_az,curr.wrong(idx),num_bin)
title([num2str(spds) ' cm/s Dors & Lat ' t_txt])

% Medial --------------------------------------------------
idx = (abs(curr.com1(:,2))<lat_pos);

subplot(2,1,2)
rose_plot(curr.azG(idx),curr.wrong(idx),num_bin)
title([num2str(spds) ' cm/s Dors & Medial ' t_txt])



function rose_plot(az,wrong,num_bin)
    
% Rose plot in correct direction
h = rose(az(~wrong),num_bin);
set(h,'Color','g')
hold on

% Rose plot of wrong direction
h = rose(az(wrong),num_bin);
set(h,'Color','r')

% Calculate and plot mean and CIs
warning off
[mu,l1,l2] = circ_mean(az);
warning on

h = polar([mu mu],[0 max(abs(ylim))]);
set(h,'Color','k')
h = polar([l1 l1],[0 max(abs(ylim))]);
set(h,'Color','k')
set(h,'LineStyle','--')
h = polar([l2 l2],[0 max(abs(ylim))]);
set(h,'Color','k')
set(h,'LineStyle','--')


function rose_wrt_body(b,r,num_bin,index,spd1,spd2,txt1,txt2,angl)
% Creates rose plots of response diretcion wrt the prey body for different
% categories that are defined by differences in flow along two regions of
% the body
      
% List of all usable sequences
seqs = find(~isnan(b.preyx2(:,2)));

% Speed values
spds = [2 11 20];
    
if strcmp(angl,'both') %Bill's addition
    dir = r.az_dir;
    wrong = r.wrong;
    
else
    error('unrecognized "angle" input -- should be "az" or "el" or "both"')
end


% Vector to keep track of which have 1 being greater than 2
one_more = nan(length(spd1),1);
    
% Loop through squences, mark where one is greater
for i = 1:length(seqs)

        % Difference in speed between 1 and 2
        Dspd = spd1{seqs(i)} - spd2{seqs(i)};
        
        % Idenx of position where greatest
        iBod = find(abs(Dspd)==max(abs(Dspd)),1,'first');
        
        % Update 'one_more'
        if Dspd(iBod)<0
            one_more(seqs(i)) = 0;
        else
            one_more(seqs(i)) = 1;
        end 
end

figure

% Step through each speed
for i = 1:3
    
    
    % Index for spd1>spd2 & 'correct' direction ---------------------------
    idx = index{i} & ~isnan(b.preyx2(:,2)) & (one_more==1) & ~wrong;
    
    % Plot 
    subplot(2,3,i)
    rose(dir(idx),num_bin)
    title([num2str(spds(i)) ' cm/s ' txt1 ' > ' txt2])
    hold on
    xlabel(angl)
    
    
    % Index for spd1>spd2 & 'wrong' direction -----------------------------
    idx = index{i} & ~isnan(b.preyx2(:,2)) & (one_more==1) & wrong;
    
    % Plot
    hR = rose(dir(idx),num_bin);
    set(hR,'Color','r')
    
    
    % Index for spd1<spd2 & 'correct' direction ---------------------------
    idx = index{i} & ~isnan(b.preyx2(:,2)) & (one_more==0) & ~wrong;
    
    % Plot 
    subplot(2,3,3+i)
    rose(dir(idx),num_bin)
    title([num2str(spds(i)) ' cm/s ' txt1 ' < ' txt2])
    xlabel(angl)
    hold on
    
    
    % Index for spd1<spd2 & 'wrong' direction -----------------------------
    idx = index{i} & ~isnan(b.preyx2(:,2)) & (one_more==0) & wrong;
    
    % Plot
    hR = rose(dir(idx),num_bin);
    set(hR,'Color','r')
    
end


function ptsT = global_to_local(rost,com,tail,pts)
% Transforms coordinates (pts, in n x 3) from global to local coordinate
% system, assuming the y-axis of larvae runs perpendicular to gravity

% Check dimensions
if size(rost,1)~=1 || size(rost,2)~=3 || size(com,1)~=1 | size(com,2)~=3 ...
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

% Translate global coordinates wrt rostrum
ptsT(:,1) = pts(:,1) - rost(1);
ptsT(:,2) = pts(:,2) - rost(2);
ptsT(:,3) = pts(:,3) - rost(3);

% Rotate points
ptsT = [inv(R) * ptsT']';

% Visualize to test
if 0
    
    blength = norm([tail(1)-rost(1) tail(2)-rost(2) tail(3)-rost(3)]);    
    
    figure
    
    subplot(2,2,[1 3])
    plot3([tail(1) rost(1)],[tail(2) rost(2)],[tail(3) rost(3)],'b',...
          rost(1),rost(2),rost(3),'bo');
    hold on
    plot3(pts(:,1),pts(:,2),pts(:,3),'ro')
    xlabel('x'); ylabel('y'); zlabel('z')
    hold off
    grid on;axis equal
    view(3)
    title('global')
    
    subplot(2,2,2)
    plot([0 blength],[0 0],'b',0,0,'ob',ptsT(:,1),ptsT(:,2),'ro')
    xlabel('x');ylabel('y')
    grid on; axis equal
    title('local')
    
    subplot(2,2,4)
    plot([0 blength],[0 0],'b',0,0,'ob',ptsT(:,1),ptsT(:,3),'ro')
    xlabel('x');ylabel('z')
    grid on; axis equal
end


function ptsT = local_to_global(rost,com,tail,pts)
% Transforms coordinates (pts, in n x 3) from global to local coordinate
% system, assuming the y-axis of larvae runs perpendicular to gravity

% Check dimensions of landmark coordinates
if size(rost,1)~=1 || size(rost,2)~=3 || size(com,1)~=1 | size(com,2)~=3 ...
   size(tail,1)~=1 || size(tail,2)~=3
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

% If points are not meshgridded
if size(pts,2)==3
    % Rotate points
    ptsT = [R * pts']';
    
    % Translate global coordinates wrt rostrum
    ptsT(:,1) = ptsT(:,1) + rost(1);
    ptsT(:,2) = ptsT(:,2) + rost(2);
    ptsT(:,3) = ptsT(:,3) + rost(3);
    
else
    error('points need to be arranged in 3 columns')
    
end

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


function [xT,yT,zT] = global_to_local_matrix(rost,com,tail,xpts,ypts,zpts)
% Transforms coordinates (pts, in n x 3) from global to local coordinate
% system, assuming the y-axis of larvae runs perpendicular to gravity

% Check dimensions
if size(rost,1)~=1 || size(rost,2)~=3 || size(com,1)~=1 || size(com,2)~=3 ...
   size(tail,1)~=1 || size(tail,2)~=3 
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

% Translate global coordinates wrt rostrum
xptsT = xpts - rost(1);
yptsT = ypts - rost(2);
zptsT = zpts - rost(3);

% Rotate points
%ptsT = [inv(R) * ptsT']';

% Transformation
for i = 1:size(xptsT,2)
    
    for j = 1:size(xptsT,3)
        ptsT = [xptsT(:,i,j) yptsT(:,i,j) zptsT(:,i,j)];
        
        % Rotate points
        ptsT = [inv(R) * ptsT']';
        
        % Store
        xT(:,i,j) = ptsT(:,1);
        yT(:,i,j) = ptsT(:,2);
        zT(:,i,j) = ptsT(:,3);    
    end    
end


function [xptsT,yptsT,zptsT] = local_to_global_matrix(rost,com,tail,xpts,ypts,zpts)
% Transforms coordinates (pts, in n x 3) from global to local coordinate
% system, assuming the y-axis of larvae runs perpendicular to gravity

% Check dimensions of landmark coordinates
if size(rost,1)~=1 || size(rost,2)~=3 || size(com,1)~=1 | size(com,2)~=3 ...
   size(tail,1)~=1 || size(tail,2)~=3
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

% Transformation
for i = 1:size(xpts,2)
    
    for j = 1:size(xpts,3)
        pts = [xpts(:,i,j) ypts(:,i,j) zpts(:,i,j)];
        
        % Rotate points
        ptsT = [R * pts']';
        
        % Store
        xptsT(:,i,j) = ptsT(:,1) + rost(1);
        yptsT(:,i,j) = ptsT(:,2) + rost(2);
        zptsT(:,i,j) = ptsT(:,3) + rost(3);    
    end    
end


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

% Index of individuals positioned ventral to predator
i_vent = com0(:,3)<=0;

% Index of L-R responses in the 'wrong' direction
i_wrongLR = prey_dir(:,2)<0;

% Index of D-V responses in the 'wrong' direction
i_wrongDV = (i_vent & (prey_dir(:,3)>0)) |  ...
    (~i_vent & (prey_dir(:,3)<=0));