function compare_piv_cfd(root)
% Relates cfd data to piv data for a particular predator speed
% 
% Note that the cfd data are massive.  We therefore makes crop the
% data prior to analysis

% Speeds of predator approach to compare
pred_spds = [2 11 20];

% Scale factor for expanding cfd data beyond piv volume 
scl_fctr = 1.3;

% Number of points for interpolating the data
num_pts = 100;


%% Define paths 

% Root paths for piv & cfd data
if nargin < 1
    root = uigetdir(pwd,'Select root');
    if root==0
        disp('Canceled')
        return
    end
end

% Check for contents
if isempty(dir([root filesep 'piv']))
    error('Directory "piv" is missing in root')
elseif isempty(dir([root filesep 'cfd']))
    error('Directory "piv" is missing in root')
elseif isempty(dir([root filesep 'morphology']))
    error('Directory "morphology" is missing in root')
end

% Morphology path
morph_path  = [root filesep 'morphology' filesep 'adult morphology.mat'];

% Load morphology data ('m')
load(morph_path)


%% Step through comparisons

if ~isempty(dir([root filesep 'piv cfd combined' filesep 'pivcfd DV.mat']))
    button = questdlg('Rerun analysis?','???','Yes','No, just plot',...
                       'Cancel','No, just plot');
else
    button = 'Yes';
end

if strcmp(button,'Cancel')
    return
    
elseif strcmp(button,'Yes')
    
    for i = 1:length(pred_spds)
        
        % Data on the frontal plane ------------------------------------------
        perspective = 'DV';
        
        % Load data
        [p,c] = load_data(root,pred_spds(i),scl_fctr,perspective);
        
        % Calculate flow fields
        [x,y,u_piv,v_piv,u_cfd,v_cfd,x_m,T_m,B_m] = ...
            calc_field(m,p,c,pred_spds(i),perspective);
        
        % Interpolation comparison
        [x_lin,y_lin,spd_cfd,spd_piv] = ...
            line_interp(num_pts,x,y,u_piv,v_piv,c,perspective);
        
        % Title for the field plots
        tit_txt = ['Frontal, ' num2str(pred_spds(i)) ' cm/s '];
        
        % Store data
        D(i).spd = pred_spds(i);
        D(i).x = x;
        D(i).y = y;
        D(i).u_piv = u_piv;
        D(i).v_piv = v_piv;
        D(i).u_cfd = u_cfd;
        D(i).v_cfd = v_cfd;
        D(i).tit_txt = tit_txt;
        D(i).x_m = x_m;
        D(i).T_m = T_m;
        D(i).B_m = B_m;
        D(i).x_lin = x_lin;
        D(i).y_lin = y_lin;
        D(i).spd_piv = spd_piv;
        D(i).spd_cfd = spd_cfd;
        
        clear p c x y u_piv v_piv u_cfd v_cfd tit_txt x_m T_m B_m x_lin y_lin
        clear spd_piv spd_cfd
        
        
        % Data on the saggital plane -------------------------------------------
        perspective = 'LR';
        
        % Load data
        [p,c] = load_data(root,pred_spds(i),scl_fctr,perspective);
        
        % Calculate flow fields
        [x,y,u_piv,v_piv,u_cfd,v_cfd,x_m,T_m,B_m] = ...
            calc_field(m,p,c,pred_spds(i),perspective);
        
        % Interpolation comparison
        [x_lin,y_lin,spd_cfd,spd_piv] = ...
            line_interp(num_pts,x,y,u_piv,v_piv,c,perspective);
        
        % Title for the field plots
        tit_txt = ['Saggital, ' num2str(pred_spds(i)) ' cm/s '];
        
        % Store data
        L(i).spd = pred_spds(i);
        L(i).x = x;
        L(i).y = y;
        L(i).u_piv = u_piv;
        L(i).v_piv = v_piv;
        L(i).u_cfd = u_cfd;
        L(i).v_cfd = v_cfd;
        L(i).tit_txt = tit_txt;
        L(i).x_m = x_m;
        L(i).T_m = T_m;
        L(i).B_m = B_m;
        L(i).x_lin = x_lin;
        L(i).y_lin = y_lin;
        L(i).spd_piv = spd_piv;
        L(i).spd_cfd = spd_cfd;
        
        clear p c x y u_piv v_piv u_cfd v_cfd  tit_txt
        
    end
    
    save([root filesep 'piv cfd combined' filesep 'pivcfd DV.mat'],'D')
    save([root filesep 'piv cfd combined' filesep 'pivcfd LR.mat'],'L')
    
else
    % Load 'D' and 'L'
    load([root filesep 'piv cfd combined' filesep 'pivcfd DV.mat'])
    load([root filesep 'piv cfd combined' filesep 'pivcfd LR.mat'])
end


%% Plot results

% Plots of spd vs. distance
plot_comparison(L,D)

max_spds = 2.*[1 1 1];

% Flow fields
for i = 1:length(L)
    % L/R perspective
    plot_fields(L(i).x,L(i).y,L(i).u_piv,L(i).v_piv,L(i).u_cfd,L(i).v_cfd,...
                L(i).x_m,L(i).T_m,L(i).B_m,L(i).spd,L(i).tit_txt,max_spds(i));
    
    % D/V perspective
    plot_fields(D(i).x,D(i).y,D(i).u_piv,D(i).v_piv,D(i).u_cfd,D(i).v_cfd,...
                D(i).x_m,D(i).T_m,D(i).B_m,D(i).spd,D(i).tit_txt,max_spds(i));
end



function plot_comparison(L,D)

figure

clr1 = {'r-','b-','g-'};
clr2 = {'r--','b--','g--'};

h = [];

for i = 1:length(L)
    
    for j = 1:3
        if i==1
            subplot(4,3,j)
            plot(L(i).x_lin(j,:),L(i).y_lin(j,:));
            axis square
            axis([0 max(L(1).x_lin(1,:)) 0 max(L(1).x_lin(1,:))])
            
            if j==2
                title('Saggital plane')
            end
        end
        
        subplot(4,3,[3+j 6+j 9+j])
        dist = sqrt(L(i).x_lin(j,:).^2 + L(i).y_lin(j,:).^2);
        h_tmp = plot(dist,L(i).spd_cfd(j,:),clr1{i},...
                      dist,L(i).spd_piv(j,:),clr2{i});
        set(h_tmp(1),'LineWidth',.5)
        set(h_tmp(2),'LineWidth',1.5)
        
        h = [h; h_tmp];
        hold on
        
        ylim([0 10])
    end
       
end

subplot(4,3,[4 7 10])
legend([h(1);h(7);h(13)],{'2','11','20'})
xlabel('Distance from mouth (cm)')
ylabel('Flow speed (cm/s)')

figure

h = [];

for i = 1:length(D)
    
    for j = 1:3
        if i==1
            subplot(4,3,j)
            plot(D(i).x_lin(j,:),D(i).y_lin(j,:));
            axis square
            axis([0 max(D(1).x_lin(1,:)) 0 max(D(1).x_lin(1,:))])
            
            if j==2
                title('Frontal plane')
            end
        end
        
        subplot(4,3,[3+j 6+j 9+j])
        dist = sqrt(D(i).x_lin(j,:).^2 + D(i).y_lin(j,:).^2);
        h_tmp = plot(dist,D(i).spd_cfd(j,:),clr1{i},...
                      dist,D(i).spd_piv(j,:),clr2{i});
        set(h_tmp(1),'LineWidth',.5)
        set(h_tmp(2),'LineWidth',1.5)
        
        h = [h; h_tmp];
        hold on
        
        ylim([0 10])
    end
       
end

subplot(4,3,[4 7 10])
legend([h(1);h(7);h(13)],{'2','11','20'})
xlabel('Distance from mouth (cm)')
ylabel('Flow speed (cm/s)')


function [x_lin,y_lin,spd_cfd,spd_piv] = line_interp(num_pts,x,y,...
                                                 u_piv,v_piv,c,perspective)

% Calculate PIV speed
U_piv = sqrt(u_piv.^2 + v_piv.^2);


% Coordinates for lines in front of the predator
x_lin = [linspace(0,max(x(:)),num_pts); ...
         sqrt(2).*linspace(0,max(x(:)),num_pts).*cos(pi/4); ...
         zeros(1,num_pts)];
     
% Step through each line
for i = 1:size(x_lin,1)
    
    % Calculate CFD speed
    if strcmp(perspective,'DV')
        y_lin = [x_lin(3,:); x_lin(2,:); x_lin(1,:)];
        z_lin = y_lin.*0;
        
        % Interpolate for PIV speed
        spd_piv(i,:) = interp2(x,y,U_piv,x_lin(i,:),y_lin(i,:));
        
        % Interpolate sparse CFD data
        warning off
        u_cfd = griddata(c.x,c.y,c.z,c.u,x_lin(i,:),y_lin(i,:),z_lin(i,:));
        v_cfd = griddata(c.x,c.y,c.z,c.v,x_lin(i,:),y_lin(i,:),z_lin(i,:));
        warning on
        
        % CFD speed
        spd_cfd(i,:) = sqrt(u_cfd.^2 + v_cfd.^2);
        
    elseif strcmp(perspective,'LR')
        z_lin = [x_lin(3,:); x_lin(2,:); x_lin(1,:)];
        y_lin = z_lin.*0;
        
        % Interpolate for PIV speed
        spd_piv(i,:) = interp2(x,y,U_piv,x_lin(i,:),z_lin(i,:));
        
        % Interpolate sparse CFD data
        warning off
        u_cfd = griddata(c.x,c.y,c.z,c.u,x_lin(i,:),y_lin(i,:),z_lin(i,:));
        w_cfd = griddata(c.x,c.y,c.z,c.w,x_lin(i,:),y_lin(i,:),z_lin(i,:));
        warning on
        
        % CFD Speed
        spd_cfd(i,:) = sqrt(u_cfd.^2 + w_cfd.^2);
        
        % Convert z values to y values
        y_lin = z_lin;
    else
        error('Requested "persepctive" not recognized')
    end
end
     

function [p,c] = load_data(root,pred_spd,scl_fctr,perspective)
% Returns the flow field structures for piv and cfd data


% Define subdirectories 
if pred_spd == 2
    
    cfd_path = [root filesep 'cfd' filesep 'flow_02cmps_around_zebrafish.ascii'];
    piv_path = [root filesep 'piv' filesep 'compiled data' filesep ...
                perspective ' view S02_piv.mat'];
      
elseif pred_spd == 11
    
    cfd_path = [root filesep 'cfd' filesep 'flow_11cmps_around_zebrafish.ascii'];
    piv_path = [root filesep 'piv' filesep 'compiled data' filesep ...
                perspective ' view S11_piv.mat'];  
            
elseif pred_spd == 20
    
    cfd_path = [root filesep 'cfd' filesep 'flow_20cmps_around_zebrafish.ascii'];
    piv_path = [root filesep 'piv' filesep 'compiled data' filesep ...
                perspective ' view S20_piv.mat'];   
else
    
    error('predator speed must be defined as 2, 11 or 20')
    
end

% Load piv data ('p')
load(piv_path)

% Load cfd velocity data ('C')
c = cfd_import_vel(cfd_path);

% Index that defines domain covered by piv data
i_piv = (c.x>min(p.x(:))*scl_fctr) & (c.x<max(p.x(:))*scl_fctr) & ...
        (c.y>min(p.y(:))*scl_fctr) & (c.y<max(p.y(:))*scl_fctr) & ...
        (c.z>min(p.z(:))*scl_fctr) & (c.z<max(p.z(:))*scl_fctr);
    
% Trim the cfd data to the piv domain
x_cfd = c.x(i_piv);
y_cfd = c.y(i_piv);
z_cfd = c.z(i_piv);
u_cfd = c.u(i_piv);
v_cfd = c.v(i_piv);
w_cfd = c.w(i_piv);


function [x,y,u_piv,v_piv,u_cfd,v_cfd,x_m,T_m,B_m] = calc_field(...
                                                m,p,c,pred_spd,perspective)
% Finds the PIV and CFD data to be compared on either the saggital or
% frontal planes

% Flow on the saggital plane
if strcmp(perspective,'LR')
    
    % Index for saggital plane values
    y_idx = p.y==0;
    
    if sum(y_idx) > 1
        error('More than one saggital plane in the data')
    end
    
    % PIV coordinates and vectors for the saggital plane
    x_piv = reshape(p.x,size(p.x,1),size(p.x,3));
    y_piv = x_piv.*0;
    z_piv = reshape(p.z,size(p.z,1),size(p.z,3));
    
    u_piv = p.u(:,y_idx,:);
    w_piv = p.w(:,y_idx,:);
    
    u_piv = reshape(u_piv,size(u_piv,1),size(u_piv,3));
    w_piv = reshape(w_piv,size(w_piv,1),size(w_piv,3));
    
    % Subtract freestream
    u_piv = u_piv + pred_spd;
    
    % Margin of the body extracted (in cm)
    x_m = (m.s-max(m.s)).*100;
    T_m = -(m.h./2 + (m.c-m.c(end))).*100;
    B_m = -((m.c-m.c(end)) - m.h./2).*100;
    
    % Spline fit to the body profile
    T_cs = csapi(x_m,T_m);
    B_cs = csapi(x_m,B_m);
    
% Flow on the frontal plane
elseif strcmp(perspective,'DV')
    
    % Index for saggital plane values
    z_idx = p.z==0;
    
    if sum(z_idx) > 1
        error('More than one frontal plane in the data')
    end
    
    % PIV coordinates and vectors for the saggital plane
    x_piv = p.x;
    y_piv = p.y;
    z_piv = p.x.*0;
    
    u_piv = p.u(:,:,z_idx);
    v_piv = p.v(:,:,z_idx);
    
    % Subtract freestream
    u_piv = u_piv + pred_spd;
   
    % Margin of the body extracted (in cm)
    x_m = (m.s-max(m.s)).*100;
    T_m = m.w./2.*100;
    B_m = -m.w./2.*100;
    
    % Spline fit to the body profile
    T_cs = csapi(x_m,T_m);
    B_cs = csapi(x_m,B_m);
    
else
    error('Requested "perspective" not identified')
end

% Interpolation of sparse CFD data
warning off
u_cfd = griddata(c.x,c.y,c.z,c.u,x_piv,y_piv,z_piv);
v_cfd = griddata(c.x,c.y,c.z,c.v,x_piv,y_piv,z_piv);
w_cfd = griddata(c.x,c.y,c.z,c.w,x_piv,y_piv,z_piv);
warning on

% Translate output to 2D for saggital and frontal planes
if strcmp(perspective,'LR')
    x     = x_piv;
    y     = z_piv;
    v_piv = w_piv;
    v_cfd = w_cfd;
else
    x     = x_piv;
    y     = y_piv;  
end
clear w_piv w_cfd x_piv y_piv z_piv  
    
% % Identify the nodes within the body
% i_morph = logical(y.*0);
% for i = 1:size(x,2)
%     if x(1,i) <= max(x_m)
%         i_morph(:,i) = (y(:,i)<fnval(T_cs,x(1,i))) & ...
%             (y(:,i)>fnval(B_cs,x(1,i)));
%     end
% end
% clear T_cs B_cs

% Zero out values within the body
%u_piv(i_morph) = nan;
%v_piv(i_morph) = nan;
%u_cfd(i_morph) = nan;
%v_cfd(i_morph) = nan;


function plot_fields(x,y,u_piv,v_piv,u_cfd,v_cfd,x_m,T_m,B_m,pred_spd,tit_txt,max_spd)

% Color to render the body color
bod_clr = .5.*[1 1 1];

% Calculate speed
spd_piv = sqrt(u_piv.^2 + v_piv.^2);
spd_cfd = sqrt(u_cfd.^2 + v_cfd.^2);

% Offset base of arrows by half a grid
x_off = mean(diff(x(1,:)))/2;
y_off = mean(diff(y(:,1)))/2;

% New figure window
figure

% Plot PIV field
subplot(1,2,1)
h1 = pcolor(x,y,log10(spd_piv));
caxis([0 log10(max_spd)])
colorbar('South')
set(h1,'EdgeColor','none')
hold on
h2 = quiver(x+x_off,y+y_off,u_piv,v_piv,'k');
axis square
axis equal
xLim = xlim;
yLim = ylim;
h3 = patch([x_m; x_m(end:-1:1)],[B_m; T_m(end:-1:1)],bod_clr);
set(h3,'EdgeColor','none')
axis([xLim yLim])
hold off
title([tit_txt '  PIV'])

clear h1 h2 h3

% Plot CFD field
subplot(1,2,2)
h1 = pcolor(x,y,log10(spd_cfd));
shading interp
caxis([0 log10(max_spd)])
colorbar('South')
set(h1,'EdgeColor','none')
hold on
h2 = quiver(x+x_off,y+y_off,u_cfd,v_cfd,'k');
axis square
axis equal
hold on 
h3 = patch([x_m; x_m(end:-1:1)],[T_m; B_m(end:-1:1)],bod_clr);
set(h3,'EdgeColor','none')
axis([xLim yLim])
title('CFD')


function comparison_plots(x_piv,y_piv,spd_piv,x_cfd,y_cfd,spd_cfd)

idx = 12;
figure;
plot(x(1,:),spd_piv(idx,:)./pred_spd,'k-',x(1,:),spd_cfd(idx,:)./pred_spd,'k--')
xlabel('X (cm)')
ylabel('Flow speed (Normalized to freestream)')
legend('piv','cfd')
xlim([0 max(x(1,:))])


