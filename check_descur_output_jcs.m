function [] = check_descur_output_jcs()
%JC's file modified by JL 11/07 to automatically parse the outcurve file
%instead of having to copy and paste lines in.  It assumes that the
%outcurve file is in the same directory as the LCFS file.


%get the R,Z coords from field line following
[filename,pathname] = uigetfile('*.mat', 'Indicate the file containing the x,y,z coordinates of the LCFS');
load([pathname '\' filename], 'surf_x', 'surf_y', 'surf_z');
rsurf = sqrt(surf_x.^2 + surf_y.^2);
zsurf = surf_z;
%
numphi = 50;

% numtormodes_first_pass = 49;
numtormodes_first_pass = 25;

%phi angle of each cut
phi = 0:2*pi/(4*numphi):pi/2;
theta = linspace(0,2*pi,1000);

[fid_outcurve]=fopen(uigetfile('*.*', 'Indicate the descur outcurve file'),'r');
%parse the file to find the RBC lines
try
    tline = fgetl(fid_outcurve)
    if length(tline) < 5
        cmpSTR = ' ';
    else
        cmpSTR = tline(4:5);
    end
    while((~strcmp(cmpSTR,'MB')))
        tline = fgetl(fid_outcurve);
        if length(tline) < 5
            cmpSTR = ' ';
        else
            cmpSTR = tline(4:5);
        end
    end
catch
    error('Something messed up when trying to read outcurve file')
end

% for i = 1:number toroidal modes used in the expansion
for ii=1:numtormodes_first_pass
% for i=1:64
    temp_data=fgetl(fid_outcurve);
    data(ii,:)=temp_data(1:57);
end

keep_going = 1;
counter = ii;
while keep_going % check the outcurve file, this ??? changes.
    temp_data=fgetl(fid_outcurve)
    if max(size(temp_data)) > 1
        counter = counter+1;
        data(counter,:)=temp_data;        
    else
        keep_going = 0;
    end    
end
data=str2num(data);
numtotalmodes = counter

m = data(:,1);
n = data(:,2);
rbc = data(:,3); %cos component of R
rbs = data(:,4); %sin component of R
zbc = data(:,5); %cos component of Z
zbs = data(:,6); %sin component of Z
% figure
% hold on
for kk = 1 % just the first surface
% for kk = 1:length(phi) %for each cut
% for k = length(phi) %for each cut
    clear r_rec z_rec
    for ii = 1:length(theta) %sum up the fourier components of R,Z for each m,n for all angles theta
        r_rec(ii) = sum(rbc.*cos(m*theta(ii) - 4*n*phi(kk)))+sum(rbs.*sin(m*theta(ii) - 4*n*phi(kk)));
        z_rec(ii) = sum(zbc.*cos(m*theta(ii) - 4*n*phi(kk)))+sum(zbs.*sin(m*theta(ii) - 4*n*phi(kk)));
    end
    figure(kk);box on; hold on;
    %plot field line follow points
    plot(rsurf(kk:numphi:end),zsurf(kk:numphi:end),'k.')
    hold on
    %plot fourier points
%     plot(r_rec,z_rec,'r.','MarkerSize',.5)
    plot(r_rec,z_rec,'r-','LineWidth',1)
    axis equal
%     pause
end

fclose(fid_outcurve)