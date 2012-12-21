function d = cfd_import_vel_mod(f_path)
% Imports velocity data from an ascii file of cfd data that includes all
% partial derivatives.

fid = fopen(f_path);

if fid==-1
    error('file cannot be opened');
end

% Read header
H = textscan(fid, ...
    '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',1,'delimiter',',');

% Read data
N = textscan(fid, ...
    '%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',',');

% Close file identifier
fclose(fid);

% Check size of structues
if (length(H)~=16) || (length(N)~=16)
    error('The ASCII file was expected to have 16 columns');
end

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


function check_str(H,i,str)
if ~strcmp(H{i},str)
    error(['Column ' num2str(i) ' expected to be ' str]);
end
