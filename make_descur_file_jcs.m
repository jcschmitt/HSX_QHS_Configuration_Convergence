function [] = make_descur_file_jcs(filename1,pathname)
%Modified by JL 11/07
% automation added by JCS 4/09

try
    if nargin < 1
        %get the cart. coords. of the lcfs (from make_lcfs_jl.m)
        [filename1,pathname] = uigetfile('*.mat', 'Indicate the file containing the x,y,z coordinates of the LCFS');
    elseif nargin < 2
        pathname = './';
    end
    cur_pwd = pwd;
    cd(pathname);
    load(filename1, ...
        'surf_x', ...
        'surf_y', ...
        'surf_z');
%     load(filename1, ...
%         'surf_x', ...
%         'surf_y', ...
%         'surf_z', ...
%         'ntheta', ...
%         'nphi' ...
%         );
    
    outfilename = ['rzdata_' filename1(6:(end-4)) '.dat'];
    %these values assume you used 160 toridal cuts and 100 pts for each cut
    %(nphi and theta respectively), in make_lcfs_jl.  Using stell. sym. you
    %actually have these values for HSX.
    ntheta = 800;
    nphi = 50;
    % nphi = nphi; % nphi is in the file
    % ntheta = ntheta; % need to take account of HSX having 4-periods
    % nphi =
    % ntheta = nphi * 4;
    % nphi = (length(surf_x) - 1 ) / ntheta;
    % for 256 cuts, and 400 pts/cut
    % ntheta = 800;
    % nphi = 50;
    % ntheta = 1600;
    % nphi = 64;
    % for  128 cuts, and 200 pts/cut
    % ntheta = 800;
    % nphi = 100;
    
    %open output file for writing
    fid = fopen(outfilename,'w');
    fprintf(fid,'%s %s %s \n',num2str(ntheta),num2str(nphi),'4');
    % R = sqrt(surf_x.^2 + surf_y.^2);
    % PHI = atan2(surf_y,surf_x);
    % Z = surf_z;
    
    for phi_ind = 1:nphi %loop over cuts and calc R,phi,Z
        Rsurf = sqrt(surf_x(phi_ind:nphi:end).^2 + surf_y(phi_ind:nphi:end).^2);
        PHI = atan2(surf_y(phi_ind:nphi:end),surf_x(phi_ind:nphi:end));
        Zsurf = surf_z(phi_ind:nphi:end);
        for theta_ind = 1:ntheta %print each point on cut
            fprintf(fid,'%f %f %f \n',Rsurf(theta_ind),PHI(theta_ind),Zsurf(theta_ind));
        end
        figure(400+phi_ind);hold on;box on;hold on;
        plot(Rsurf, Zsurf,'.');
    end
    
    %shouldn't need this if you actually close the file
    %
    % for i = 1:100
    %     fprintf(fid,'%s \n','Finish This');
    % end
    
    fclose(fid);
    cd(cur_pwd);
catch
    warning('error in make_descur_file_jcs.m');
    filename1
    pathname
    cd(cur_pwd);
end