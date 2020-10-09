function [eps_v_r, eps_v_i,kappa_s, kappa_a, sigma_0_vv, sigma_0_vv_canopy, sigma_0_vv_gcg, sigma_0_vv_gc]=function_TS_backscattering_estimation(S2_LAI_TS, S2_CW_TS)

%% 
C_water = 10.*S2_LAI_TS.*S2_CW_TS;
mv = C_water * 0.3988;%C_water * leaf thickness * water density
rho = 0.3;
mg = mv./(mv + (1 - mv)*rho);
%
h =[10,15,16,17,26,30,35.0326,58.6863,71.4251,82.6688,82.1204,81.1737,16,16];% height of the wheat
h = h/100;
%
f = 5.405; %Frequency (GHz)
s = 0.025; %rms height (m)
% a = 0.015;
% kappa_e = 1.695;
theta = 39.5925;

% compute real and imaginery part
for i =1:1:size(mg,2)
    [eps_v_r(i), eps_v_i(i)] = RelDielConst_Vegetation(f, mg(i));%4-9.2 Dielectric Model
    [kappa_e(i), kappa_s(i), kappa_a(i)] = function_extinction_scattering_Rayleigh(eps_v_r(i), eps_v_i(i));
    [sigma_0_vv(i), sigma_0_vv_canopy(i), sigma_0_vv_gcg(i), sigma_0_vv_gc(i)] = ...
        S2RTR_DiffuseUB_Kas(eps_v_r(i), eps_v_i(i), h(i), f, s, kappa_a(i), kappa_s(i), theta);
end




end