function CH4_circular

%% load b


if ~isempty(dir('/Users/williamstewart/Dropbox/Robopredator/behavior/Transformed_Prey_Coords.mat'))
    load('/Users/williamstewart/Dropbox/Robopredator/behavior/Transformed_Prey_Coords.mat');
    
else
    load('/Users/mmchenry/Dropbox/Projects/Robopredator/behavior/Transformed_Prey_Coords.mat')
end
    
    
%%


x = b.preyx(:,2);
y = b.preyy(:,2);
z = b.preyz(:,2);

x2 = b.preyx2(:,2);
y2 = b.preyy2(:,2);
z2 = b.preyz2(:,2);

%mirror responses that are at position y < -0.5;
who = y < -0.5;
y(who) = abs(y(who));
y2(who) = abs(y2(who));


%% calculate azimuth

speed = [2 11 20];

for i = 1:3

    %w = w = whose(speed(i)) & ~isnan(b.preyx2(:,2)) & (b.preyy(:,2) <= 0.5);
    w = whose(b,speed(i)) & ~isnan(b.preyx2(:,2)) & (b.preyy(:,2) > 0.5);

    dx = x2(w)- x(w);
    dy = y2(w) - y(w);
    dz = z2(w) - z(w);


    [az,r] = cart2pol(dx,dy);

    theta = az;
    
    [mu,l1,l2] = circ_mean(theta);

    disp(['speed = ' num2str(speed(i))]);
    l1*(180/pi)
    mu*(180/pi)
    l2*(180/pi)
    
    
    
    
end


%% calculate elevation

speed = [2 11 20];

for i = 1:3

    %w = whose(speed(i)) & ~isnan(b.preyx2(:,2)) & (b.preyz(:,2) > 0.5);
    w = whose(b,speed(i)) & ~isnan(b.preyx2(:,2)) & ...
        (b.preyz(:,2) <= 0.5) & (b.preyz(:,2) >= -0.5);
    %w = whose(speed(i)) & ~isnan(b.preyx2(:,2)) & (b.preyz(:,2) < -0.5);

    dx = x2(w)- x(w);
    dy = y2(w) - y(w);
    dz = z2(w) - z(w);


    [el,r] = cart2pol(dx,dz);

    theta = el;
    
    [mu,l1,l2] = circ_mean(theta);

    disp(['speed = ' num2str(speed(i))]);
    l1*(180/pi)
    mu*(180/pi)
    l2*(180/pi)
    
    
    
    
end


function [w_out] = whose(b,speed,light,LL,respond)
% function [who] = whose(speed,light,LL,respond);
%
%creates logical for different treatments of robopredator data
%
%only specify respond if interested in non-responders


%load('/Users/williamstewart/Dropbox/Robopredator/behavior/Transformed_Prey_Coords.mat')


commandwindow;

% if nargin == 3
%     
%     who = (isfinite(b.preyx) & (b.speed == speed) & ...
%         (b.lit == light) & (s.LL == LL));
%     
% end
% 
% if (nargin == 4) & (respond == 0) 
%     
%     who = (isnan(s.DorFSframe) & (s.speed == speed) & ...
%     (s.lit == light) & (s.LL == LL));
%     
% end

if nargin == 2
    
    w_out = (isfinite(b.preyx(:,1))) & (b.speed == speed) & ...
    (b.lit == 0) & (b.LL == 1);

end

%if no input arguments just specify the prey that have stage 2 data
if nargin == 1
    
    w_out = (isfinite(b.preyx2(:,1))) & (b.lit == 0);

end