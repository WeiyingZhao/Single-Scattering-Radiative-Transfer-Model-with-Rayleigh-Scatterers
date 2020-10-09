
function [sig_s_vv, sig_s_hh, tmp] = PRISM1_ForwardModel(eps,theta_deg,s,f) 

%     var result = new Array();

%     var a=[];

    
    ks = s * (2*pi *f/0.3); %// calculate roughness parameter

    theta = theta_deg * (pi/180.0); %// incidence angle in radian

    gamma0 = Fresn_Refl0(eps); %//reflectivity (normal incidence)

    a = Fresn_Refl(eps, theta); 
    gammav = a(1);
    gammah = a(2);


    p =  1.- (2.*theta/pi)^(1./(3.*gamma0)) * exp(-ks); 
    p = p*p;

    q = 0.23 * sqrt(gamma0) *( 1. - exp(-ks));

    g = 0.70 *(1. - exp(-0.65 * ks^1.8));

    sigvv = g * cos(theta)^3 /sqrt(p) * (gammav + gammah);

    sig_s_vv = 10.*log(sigvv)/log(10);     %// vv
    sig_s_hh = 10.*log(sigvv * p)/log(10); %// hh
    tmp = 10.*log(sigvv * q)/log(10); %// vh

end


%Code 10.5: PRISM-1  - Forward Model

%Description: Code computes sigma_0 for all three polarization
%combinations, given the surface parameters.

%Input Variables:
    %eps = eps' - j eps'': Complex dielectric constant of the scattering
    %medium
    %theta: Incidence angle (degrees)
    %s: rms height (m)
    %f: Frequency (GHz)
    
%Output Products:
    %sigma_0_vv (dB)
    %sigma_0_hh (dB)
    %sigma_0_hv (dB)
    
%Book Reference: Section 10-5

%Matlab Code: 
% 
% function [sig_0_vv sig_0_hh sig_0_hv] = PRISM1_ForwardModel(eps,theta,s,f)
% 
% ks = s .* (2*pi *f/0.3); % calculate roughness parameter
% 
% theta = theta .* (pi/180.0); % incidence angle in radian
% 
% gamma0 = Fresn_Refl0(eps); %reflectivity (normal incidence)
% 
% [gammav, gammah] = Fresn_Refl(eps, theta); 
% 
% p = ( 1- (2 .*theta/pi).^(1./(3*gamma0)) .* exp(-ks) ).^2; 
% q = 0.23 .* sqrt(gamma0) .*( 1 - exp(-ks));
% 
% g = 0.70 .*(1 - exp(-0.65 .* ks.^1.8));
% 
% sigvv = g .* (cos(theta)).^3 ./sqrt(p) .* (gammav + gammah);
% 
% sig_0_vv = 10*log10(sigvv);
% sig_0_hh = 10*log10(sigvv .* p);
% sig_0_hv = 10*log10(sigvv .* q);
% 
% end
% 
% function gamma0 = Fresn_Refl0(eps)
% % calculates Fresnel reflectivity at normal incidence.
% gamma0 = (abs((1 - sqrt(eps)) ./ ( 1+ sqrt(eps)))).^2;
% end
% 
% %-----------------------------------------------------------------
% function [gammav, gammah] = Fresn_Refl(eps,theta)
% % calculates Fresnel reflectivities of v and h-polarizations at given set
% % of incidence angles.
% 
% [rho_v, rho_h] = refl_coef(theta, 1, eps);
% gammav = (abs(rho_v)).^2;
% gammah = (abs(rho_h)).^2;
% 
% end
% 
% %----------------------------------------------------------------
% function [rho_v, rho_h] = refl_coef(the1, eps1, eps2)
% 
% % calculates the v and h-polarized reflection coefficients of a plane
% % dielectric surface
% 
% n1 = sqrt(eps1);
% n2 = sqrt(eps2);
% costh2 = sqrt(1 - (n1.*sin(the1)./n2).^2);
% 
% 
% rho_v = -(n2.*cos(the1) - n1.*costh2)./(n2.*cos(the1)+n1.*costh2);
% rho_h = (n1.*cos(the1) - n2.*costh2)./(n1.*cos(the1)+n2.*costh2);
% end