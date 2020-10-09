%Code 11.1: S2RT/R Backscattering from Rayleigh Layer with Diffuse Upper
%Boundary

%Description: Code computes sigma_0_vv and sigma_0_hh for a weakly
%scattering Rayleigh layer with albedo a < 0.2. The layer is above a ground
%surface characterized by the PRISM model (code 10.5).

%Input Variables:
    %eps: complex dielectric constant of ground surface
    %f: frequency (GHz)
    %s: rms height of ground surface (m)
    %a: Single-scattering albedo of Rayleigh layer (unitless)
    %kappa_e: extinction coefficient of the Rayleigh layer (Np/m)
    %d: thickness of the Rayleigh layer (m)
    %theta: Incidence angle (degrees)

%Output Variables: 
    %sigma_0_vv (dB)
    %sigma_0_hh (dB)
    
%Book reference: Section 11-2.2 and eq. 11.23 with n = 2

%Matlab code:

function [sigma_0_vv, sigma_0_hh] = S2RTR_DiffuseUB1(eps_v_r, eps_v_i, h)
% [eps_v_r, eps_v_i] = RelDielConst_Vegetation(5.405,0.5);
% eps = eps_v_r - eps_v_i*j;
% eps,f,s,a,kappa_e,d, theta
eps = eps_v_r - eps_v_i*j;
f = 5.405;
s = 0.025;

d = h;
a = 0.015;
kappa_e = 1.695;

theta = 34.81;

theta_r = theta *pi/180; %transform to radian

kappa_s = a .* kappa_e ; %scattering coefficient

%-- call the PRISM-1 surface scattering model
[sig_s_vv, sig_s_hh, tmp] = PRISM1_ForwardModel(eps,theta,s,f);
sig_s_vv  = 10.^(sig_s_vv/10); 
sig_s_hh  = 10.^(sig_s_hh/10); 

%-- calculate transmissivity in layer
tau = kappa_e .*d .* sec(theta_r);
T = exp(-tau);

%-- calculate reflectivity of lower interface
[t1, t2, gammah, gammav, t3, t4, t5, t6] = ReflTransm_PlanarBoundary(1, eps, theta);


%-- calculate the total backscattering coefficient according to eq 11.23
% sigma_0_vv = T.^2 .*sig_s_vv + 0.75*a * cos(theta_r) .*(1- T.^2) ...
%     .*(1 + gammav.^2 .* T.^2) + 3 * 2 * kappa_s * d * gammav * T.^2;
sigma_0_vv = 0.75*a * cos(theta_r) .*(1- T.^2) ...
    .*(1 + gammav.^2 .* T.^2) + 3 * 2 * kappa_s * d * gammav * T.^2;
sigma_0_hh = T.^2 .*sig_s_hh + 0.75*a * cos(theta_r) .*(1- T.^2) ...
    .*(1 + gammah.^2 .* T.^2) + 3 * 2 * kappa_s * d * gammah * T.^2;


sigma_0_vv = 10* log10(sigma_0_vv);
sigma_0_hh = 10* log10(sigma_0_hh);


end