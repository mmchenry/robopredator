function d = cfd_import_vel(f_path)
% Imports velocity data from an ascii file of cfd data

fid = fopen(f_path);

if fid==-1
    error('file cannot be opened');
end

% Read header
H = textscan(fid, ...
    '%s %s %s %s %s %s %s %s',1,'delimiter',',');

% Read data
N = textscan(fid, ...
    '%d %f %f %f %f %f %f %f','delimiter',',');

% Close file identifier
fclose(fid)

% Check size of structues
if (length(H)~=8) || (length(N)~=8)
    error('The ASCII file was expected to have 8 columns');
end

% Check column headings
check_str(H,1,'nodenumber');
check_str(H,2,'x-coordinate');
check_str(H,3,'y-coordinate');
check_str(H,4,'z-coordinate');
check_str(H,5,'standard-deformation-in-time-inverse');
check_str(H,6,'u-in-stationary-rf');
check_str(H,7,'y-velocity');
check_str(H,8,'z-velocity');

% Store data
d.node = N{1};
d.x    = N{2};
d.y    = N{3};
d.z    = N{4};
d.u    = N{6};
d.v    = N{7};
d.w    = N{8};


function check_str(H,i,str)
if ~strcmp(H{i},str)
    error(['Column ' num2str(i) ' expected to be ' str]);
end
