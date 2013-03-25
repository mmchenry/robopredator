function [who] = whose(speed,light,LL,respond)
% function [who] = whose(speed,light,LL,respond);
%
%creates logical for different treatments of robopredator data
%
%only specify respond if interested in non-responders


load('/Users/williamstewart/Dropbox/Robopredator/behavior/Transformed_Prey_Coords.mat')


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

if nargin == 1
    
    who = (isfinite(b.preyx(:,1))) & (b.speed == speed) & ...
    (b.lit == 0) & (b.LL == 1);

end

%if no input arguments just specify the prey that have stage 2 data
if nargin == 0
    
    who = (isfinite(b.preyx2(:,1))) & (b.lit == 0);

end