function d = cfd_import_vel(f_path,type)
% Imports velocity data from an ascii file of cfd data
% 'type' - defines whether the file includes all partial derivatives, or not
%          possibile inputs: 'short' or 'full'

%% Get file ID, check inputs

fid = fopen(f_path);

if fid==-1
    error('file cannot be opened');
end

if nargin < 2
    error('you need to specify either "short" or "full" for type')
    
elseif strcmp(type,'short')
    type_val = 1;
    
elseif strcmp(type,'full')
    type_val = 2;
    
else
    error('"type" input not recognized')
    
end

%% Read data

if type_val == 1
    
    % Read header
    H = textscan(fid, ...
        '%s %s %s %s %s %s %s %s',1,'delimiter',',');
    
    % Read data
    N = textscan(fid, ...
        '%d %f %f %f %f %f %f %f','delimiter',',');
    
    % Check size of structues
    if (length(H)~=8) || (length(N)~=8)
        error('The ASCII file was expected to have 8 columns');
    end
    
else
    
    % Read header
    H = textscan(fid, ...
        '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',1,'delimiter',',');
    
    % Read data
    N = textscan(fid, ...
        '%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');
    
    % Check size of structues
    if (length(H)~=16) || (length(N)~=16)
        error('The ASCII file was expected to have 16 columns');
    end
    
end

% Close file identifier
fclose(fid);



%% Translate data

if type_val == 1
    
% Check column headings
check_str(H,1,'nodenumber');
check_str(H,2,'x-coordinate');
check_str(H,3,'y-coordinate');
check_str(H,4,'z-coordinate');
check_str(H,5,'standard-deformation-in-time-inverse');
check_str(H,6,'u-in-stationary-rf');
check_str(H,7,'y-velocity');
check_str(H,8,'z-velocity');

% Store data (in cm, cm/s)
d.node = N{1};
d.x    = N{2}.*100;
d.y    = N{3}.*100;
d.z    = N{4}.*100;
d.u    = N{6}.*100;
d.v    = N{7}.*100;
d.w    = N{8}.*100;

else
    % Check column headings
    check_str(H,1,'nodenumber');
    check_str(H,2,'x-coordinate');
    check_str(H,3,'y-coordinate');
    check_str(H,4,'z-coordinate');
    check_str(H,5,'u-in-stationary-rf');
    check_str(H,6,'y-velocity');
    check_str(H,7,'z-velocity');
    check_str(H,8,'dx-velocity-dx');
    check_str(H,9,'dy-velocity-dx');
    check_str(H,10,'dz-velocity-dx');
    check_str(H,11,'dx-velocity-dy');
    check_str(H,12,'dy-velocity-dy');
    check_str(H,13,'dz-velocity-dy');
    check_str(H,14,'dx-velocity-dz');
    check_str(H,15,'dy-velocity-dz');
    check_str(H,16,'dz-velocity-dz');
    
    % Store data (in cm, cm/s)
    d.node = N{1};
    d.x    = N{2}.*100;
    d.y    = N{3}.*100;
    d.z    = N{4}.*100;
    d.u    = N{5}.*100;
    d.v    = N{6}.*100;
    d.w    = N{7}.*100;
    d.du_dx = N{8};
    d.dv_dx = N{9};
    d.dw_dx = N{10};
    d.du_dy = N{11};
    d.dv_dy = N{12};
    d.dw_dy = N{13};
    d.du_dz = N{14};
    d.dv_dz = N{15};
    d.dw_dz = N{16};

end
    





function check_str(H,i,str)
if ~strcmp(H{i},str)
    error(['Column ' num2str(i) ' expected to be ' str]);
end
