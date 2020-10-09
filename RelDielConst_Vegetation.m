%Code 4.8: Relative Dielectric Constant of Vegetation

%Description: Code computes the real and imaginary parts of the relative
    %dielectric constant of vegetation material, such as corn leaves, in 
    %the microwave region.   
    
%Input Variables:
    %f: frequency in GHz
    %mg: Gravimetric moisture content 0< mg< 1
    
%Output Products:
    %eps_v_r: real part of dielectric constant
    %eps_v_i: imaginary part of dielectric constant
%Book Reference: Section 4-9

%Example call: [A B] = RelDielConst_Vegetation(f,mg)
    %Computes the real and imaginary components of the permitivity of
    %vegetation vased on frequency and gravimetric moisture content.
    
%MATLAB Code

function [eps_v_r, eps_v_i] = RelDielConst_Vegetation(f,mg)

S = 15; % salinity

%-- free water in leaves

sigma_i = 0.17.*S - 0.0013 .* S.^2; 

eps_w_r = 4.9 + 74.4 ./( 1 + (f/18).^2);

eps_w_i =  74.4 .*(f/18) ./( 1 + (f/18).^2) + 18*sigma_i ./f ;

% bound water in leaves

eps_b_r = 2.9 + 55*(1+ sqrt(f/0.36))./( (1+ sqrt(f/0.36)).^2 + (f/0.36));

eps_b_i = 55*sqrt(f/0.36) ./ ( (1+ sqrt(f/0.36)).^2 + (f/0.36));


% empirical fits
v_fw = mg .*( 0.55 * mg - 0.076);
v_bw = 4.64 .*mg.^2 ./(1 + 7.36 * mg.^2); 

eps_r = 1.7 - 0.74 *mg + 6.16 .* mg.^2;


eps_v_r = eps_r + v_fw .* eps_w_r + v_bw .*eps_b_r;
eps_v_i =  v_fw .* eps_w_i + v_bw .*eps_b_i;

  end