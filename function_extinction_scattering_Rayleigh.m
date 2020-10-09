function [k_e, k_s, k_a] = function_extinction_scattering_Rayleigh(eps_v_r, eps_v_i)


v = 0.005;%the volume fraction of scatterers
% v = (4/3*pi*r^3) * Nv;
lambda= 0.056;%lambda is the wavelength (in metres)
k = 2*pi/lambda;%
r =  0.00186;% 0.0001;%radius is the disc radius (metres) 1cm for partile radius
epsilon = eps_v_r + eps_v_i*1i;
K1 = (epsilon - 1)/(epsilon + 2);

k_s = 2*v*k^4*r^3*(abs(K1))^2;
k_a = v*k*eps_v_i.*(abs(3./(epsilon+2)))^2;

k_e = k_a + k_s;

end