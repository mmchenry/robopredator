function [velocity_mag, shear_deform, out] =...
    CH4flowquery(path,speed,input_x,input_y,input_z)


%path = '/Users/williamstewart/Dropbox/Robopredator/cfd/';
%path = '/Users/mmchenry/Dropbox/Projects/Robopredator/cfd/';

%% GATHERING DATA

%Input values
%Input data needs to be in vector form (1xn or nx1 size)
%input_x = [1];
%input_y = [1];
%input_z = [1];

%Output arrays
%The output data will be formatted in a 1x3 structure array.
%The indices correspond to the following:
%   index = 1 => speed = 2
%   index = 2 => speed = 11
%   index = 3 => speed = 20

% Each index of the structure array has their own x, y and z fields for the
% output data.
values = {[]};
out = struct('u',values,'v',values,'w',values,'du_dx',values,...
    'du_dy',values,'du_dz',values,'dv_dx',values,'dv_dy',values,'dv_dz',values,...
    'dw_dx',values,'dw_dy',values,'dw_dz',values);


if speed == 2
    %Speed = 2 cm/s
    cfd_file = [path 'flow_02cmps_around_zebrafish.ascii'];
    c = cfd_import_vel_mod(cfd_file);
elseif speed == 11
    %Speed = 11 cm/s
    cfd_file = [path 'flow_11cmps_around_zebrafish.ascii'];
    c = cfd_import_vel_mod(cfd_file);
else
    %Speed = 20 cm/s
    cfd_file = [path 'flow_20cmps_around_zebrafish.ascii'];
    c = cfd_import_vel_mod(cfd_file);
end

warning off
%Velocity
out(1).u = griddata(c.x,c.y,c.z,c.u,input_x,input_y,input_z);
out(1).v = griddata(c.x,c.y,c.z,c.v,input_x,input_y,input_z);
out(1).w = griddata(c.x,c.y,c.z,c.w,input_x,input_y,input_z);
%Velocity gradients
out(1).du_dx = griddata(c.x,c.y,c.z,c.du_dx,input_x,input_y,input_z);
out(1).du_dy = griddata(c.x,c.y,c.z,c.du_dy,input_x,input_y,input_z);
out(1).du_dz = griddata(c.x,c.y,c.z,c.du_dz,input_x,input_y,input_z);
out(1).dv_dx = griddata(c.x,c.y,c.z,c.dv_dx,input_x,input_y,input_z);
out(1).dv_dy = griddata(c.x,c.y,c.z,c.dv_dy,input_x,input_y,input_z);
out(1).dv_dz = griddata(c.x,c.y,c.z,c.dv_dz,input_x,input_y,input_z);
out(1).dw_dx = griddata(c.x,c.y,c.z,c.dw_dx,input_x,input_y,input_z);
out(1).dw_dy = griddata(c.x,c.y,c.z,c.dw_dy,input_x,input_y,input_z);
out(1).dw_dz = griddata(c.x,c.y,c.z,c.dw_dz,input_x,input_y,input_z);
warning on




%% CALCULATE VELOCITY MAGNITUDE AND SHEAR DEFORMATION


%Velocity Magnitude
velocity_mag = nan(1,length(out(1).u));

velocity_mag(1,:) = (out(1).u.^2 + out(1).v.^2 + out(1).w.^2).^0.5;


%Shear Deformation
shear_deform = nan(1,length(out(1).u));

shear_deform(1,:) = (2*(out(1).du_dx).^2 + 2*(out(1).dv_dy).^2 ...
    + 2*(out(1).dw_dz).^2 + (out(1).du_dy+out(1).dv_dx).^2 ...
    + (out(1).du_dz+out(1).dw_dx).^2 + (out(1).dv_dz+out(1).dw_dy).^2).^0.5;


