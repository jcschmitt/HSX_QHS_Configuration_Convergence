function [] = calc_toroidal_flux()

load LCFS_QHS_Rstart_eq_1_509;
RR = sqrt(surf_x.^2 + surf_y.^2);
ZZ = surf_z;
flux(1) = calculateFlux(10722, [0 0 0 0 0 0], RR(1:50:end), ZZ(1:50:end));

load LCFS_QHS_Rstart_eq_1_510;
RR = sqrt(surf_x.^2 + surf_y.^2);
ZZ = surf_z;
flux(2) = calculateFlux(10722, [0 0 0 0 0 0], RR(1:50:end), ZZ(1:50:end));

load LCFS_QHS_Rstart_eq_1_511;
RR = sqrt(surf_x.^2 + surf_y.^2);
ZZ = surf_z;
flux(3) = calculateFlux(10722, [0 0 0 0 0 0], RR(1:50:end), ZZ(1:50:end));

load LCFS_QHS_Rstart_eq_1_512;
RR = sqrt(surf_x.^2 + surf_y.^2);
ZZ = surf_z;
flux(4) = calculateFlux(10722, [0 0 0 0 0 0], RR(1:50:end), ZZ(1:50:end));

load LCFS_QHS_Rstart_eq_1_513;
RR = sqrt(surf_x.^2 + surf_y.^2);
ZZ = surf_z;
flux(5) = calculateFlux(10722, [0 0 0 0 0 0], RR(1:50:end), ZZ(1:50:end));

load LCFS_QHS_Rstart_eq_1_514;
RR = sqrt(surf_x.^2 + surf_y.^2);
ZZ = surf_z;
flux(6) = calculateFlux(10722, [0 0 0 0 0 0], RR(1:50:end), ZZ(1:50:end));

load LCFS_QHS_Rstart_eq_1_515;
RR = sqrt(surf_x.^2 + surf_y.^2);
ZZ = surf_z;
flux(7) = calculateFlux(10722, [0 0 0 0 0 0], RR(1:50:end), ZZ(1:50:end));

load LCFS_QHS_Rstart_eq_1_516;
RR = sqrt(surf_x.^2 + surf_y.^2);
ZZ = surf_z;
flux(8) = calculateFlux(10722, [0 0 0 0 0 0], RR(1:50:end), ZZ(1:50:end));

save tor_flux
