function [flux] = calculateFlux(current, taper, R, z)
% Returns the flux through a surface defined by [phi, coords].  Assumes the
% surface is at the boxport.  The values in 'phi' are used to determine
% which boxport.

TOL_FLUX = 1e-9;
POLYFIT_ORDER = 40;

% Step 1: Separate the points into two regions of integration
%       Identify the apex of the surface, in terms of z.  
%       Identify the points to the each side of the apex (inboard and outboard w.r.t. the apex).
%       Each set of points is fit with a high-order polynomial, and this
%       polynomial is used to define the flux boundary during surface
%       integration.

[z_apex, index_apex] = max(z); R_apex = R(index_apex);
index_inboard = find(R <= R_apex); index_outboard = find(R >= R_apex);
R_inboard = R(index_inboard);  R_outboard = R(index_outboard);
z_inboard = z(index_inboard);  z_outboard = z(index_outboard);

poly_inboard = polyfit(z_inboard, R_inboard, POLYFIT_ORDER);
poly_outboard = polyfit(z_outboard, R_outboard, POLYFIT_ORDER);

% check the fit properties
if ( min([length(R_inboard) length(R_outboard)]) < POLYFIT_ORDER)
    warning('Polynomial order is too high in calculateFlux');
    flux = -1;
    return;
end

if 1
    z_dense = linspace(-.25, .25, 1e4);
    figure; hold on;
    plot(R_inboard, z_inboard, 'r.');
    plot(R_outboard, z_outboard, 'b.');
    plot(polyval(poly_inboard, z_dense), z_dense,'r');
    plot(polyval(poly_outboard, z_dense), z_dense,'b');
    axis([1.3 1.55 -.3 .3]);
    pause(1)
end

% Step 2: Calculate the surface integral to determine the flux through the surface
flux = dblquad(@flux_integrand, min(z) - 5e-4, max(z) + 5e-4, min(R) - 5e-4, max(R) + 5e-4, TOL_FLUX, [], poly_inboard, poly_outboard, R_apex, current, taper);

%//////////////////////////////////////////////////////////////////////////
% End of calculateFlux
%//////////////////////////////////////////////////////////////////////////

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Auxilliary functions
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function flux = flux_integrand(z, R, poly_inboard, poly_outboard, R_apex, current, taper)

flux = zeros(size(z)); 
phi = 0; %

if (R < R_apex)  % Consider the inboard points
    R_fitline = polyval(poly_inboard, z);    
    index_in = find(R >= R_fitline);   % points inside the flux surface
    Z_in = z(index_in);
    if ( length(index_in) >=1 )
        for ii = 1:length(index_in)
            [B_x, B_y, B_z] = bs_derivs_aux_mex(R, phi, z(ii), current, taper);
%             [B_x, B_y, B_z] = bs_derivs_aux_mex(R, phi, z(ii), earthField, current, taper);
% 		    B_phi = -B_x * sin(phi) + B_y * cos(phi);
		    B_phi = B_y;
            flux(index_in(ii)) = B_phi;
        end
    end
else  % Consider the outboard points
    R_fitline = polyval(poly_outboard, z);    
    index_in = find(R <= R_fitline); % points inside the flux surface
    Z_in = z(index_in);
    if length(index_in) >= 1
        for ii = 1:length(index_in)
            [B_x, B_y, B_z] = bs_derivs_aux_mex(R, phi, z(ii), current, taper);
%             [B_x, B_y, B_z] = bs_derivs_aux_mex(R, phi, z(ii), earthField, current, taper);
% 		    B_phi = -B_x * sin(phi) + B_y * cos(phi);
		    B_phi = B_y;
            flux(index_in(ii)) = B_phi;
        end
    end
end

