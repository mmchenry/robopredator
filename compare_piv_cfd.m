function compare_piv_cfd(pred_spd)
% Relates cfd data to piv data for a particular predator speed
%
% pred_spd - integer with value of 2, 11 or 20
% 
% Note that the cfd data are massive.  It therefore makes sense to crop the
% data and/or downsample prior to many analyses or visualizations.



%% Define paths 

% Root paths for piv & cfd data
cfd_root = '/Users/mmchenry/Dropbox/Projects/Robopredator/cfd';
piv_root = '/Users/mmchenry/Dropbox/Projects/Robopredator/piv';

% Define subdirectories 
if pred_spd == 2
    
    cfd_path = [cfd_root filesep 'flow_02cmps_around_zebrafish.ascii'];
    
elseif pred_spd == 11
    
    cfd_path = [cfd_root filesep 'flow_11cmps_around_zebrafish.ascii'];
    
elseif pred_spd == 20
    
    cfd_path = [cfd_root filesep 'flow_20cmps_around_zebrafish.ascii'];
    
else
    
    error('predator speed must be defined as 2, 11 or 20')
    
end


% Load cfd velocity data
d = cfd_import_vel(cfd_path);




