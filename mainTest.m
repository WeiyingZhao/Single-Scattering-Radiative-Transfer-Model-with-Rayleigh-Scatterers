% (1)Import LAI and CW images
% (2)Interpolate
% (3)Time series sigma_vv computation over HS01
% (4)Extinction and scattering in a Rayleigh medium
% (5)Compute canopy backscattering with k_a and k_s of HS01
% (6)Time series canopy backscattering computation with empirical a and k_e
% (7)bacscattering value estimation with FC01 time series images
% (8)bacscattering value estimation with FC01 time series data
% (9)Soil permittivity estimation based on NCP insitu data
% (10)Soil backscattering estimationd
% (11)Gravimetric water content computation
% (12)Time series transmitivities computation



%% (9)Soil permittivity estimation based on NCP insitu data
%Input Variables:
    %f: frequency in GHz
    %t: temperature in C
    %rho_b: bulk density in g/cm3 (typical value is 1.7 g/cm3)
    %S: Sand Fraction ( 0< S< 1)
    %C: Clay Fraction ( 0< C< 1)
    %mv: Volumetric Water content 0<mv<1
%% Preparing mv and temperature

%% Soil_permittivity_HS01 
%HS01_mv: [0.1103, 0.2639, 0.1109]
% [0.1103, 0.1103, 0.2461, 0.2667, 0.2133, 0.1393]             
%HS01_time: [3/29/2017, 5/4/2017, 5/31/2017]   
%HS01_Temperature: [20.85, 21.69, 23.71, 25.93, 31.97, 30.75]
%S2_time [3/9/2017, 3/29/2017, 4/18/2017, 4/28/2017, 5/18/2017, 5/28/2017]
%        [736763      736783      736803      736813      736833      736843]
% S2_LAI = readTimeValueNew('S2Henshui_2017_lai_HS01.csv'); % S2_LAI(2,:)
% 
% x = [736783, 736819, 736846];
% v = [0.1103, 0.2639, 0.1109];
% xq = 736783:1:736846;
% vq = interp1(x,v,xq,'spline');
% vq(find(xq==736843))


%%
clear
addpath '.\data\'
addpath C:\ADataPreparation
S2_LAI = readTimeValueNew('S2Henshui_2017_lai_HS01.csv');
%%
theta = 34.81;
theta_r = theta *pi/180; %transform to radian
f = 5.405; %Frequency (GHz)
rho_b = 1.7;
S = 0.25;
rms = 0.01;%0.025;%rms
C = 0.25;
mv = [0.1103, 0.1103, 0.2461, 0.2667, 0.2133, 0.1393];
mv = mv*1.7;
T = [20.85, 21.69, 23.71, 25.93, 31.97, 30.75];
[M,N] = size(mv);
for i = 1:1:N
    [epsr(i), epsi(i)] = RelDielConst_Soil(f,T(i),rho_b, mv(i),S,C);
    eps = epsr(i) - epsi(i)*j;
    [sig_s_vv(i), sig_s_hh, tmp] = PRISM1_ForwardModel(eps,theta_r,rms,f);
end

figure;
subplot(1,4,1);plot(S2_LAI(1,:), mv);datetick('x','mmm');xlabel('Time');ylabel('Volumetric moisture');
subplot(1,4,2);plot(S2_LAI(1,:), T);datetick('x','mmm');xlabel('Time');ylabel('Temperature (degree)');
subplot(1,4,3);
plot(S2_LAI(1,:), epsr);hold on;
plot(S2_LAI(1,:), epsi);
legend('epsilon-r', 'epsilon-i')
datetick('x','mmm');xlabel('Time');ylabel('Permittivity of soil');
subplot(1,4,4);plot(S2_LAI(1,:), sig_s_vv);
datetick('x','mmm');xlabel('Time');ylabel('Bacscattering (dB)');



%% Soil_permittivity_FC01
%FC01_mv: [0.2769, 0.4525, 0.0751]
%FC01_mv = [0.2769,0.2769,0.2769,0.2769,0.2769, 0.3477, 0.4617, 0.4805, 0.3386, 0.2324, 0.1779, 0.0751,0.0751,0.0751]
%FC01_time: [2017/3/31, 5/6/2017, 2017/06/02]
%FC01_Temperature: [14.23, 17.61, 16.21, 21.03, 16.17, 18.77, 22.97, 22.41, 25.85, 31.33, 25.75, 28.19, 30.05, 30.49]

% 
% x = [736785, 736821, 736848];
% v = [0.2769, 0.4525, 0.0751];
% xq = 736785:1:736848;
% vq = interp1(x,v,xq,'spline');
% vq(find(xq==736843))
%%
clear
addpath '.\data\'
addpath C:\ADataPreparation
S2_LAI = readTimeValueNew('S2Henshui_2017_lai_FC01.csv');
%%
theta = 39.5925;
theta_r = theta *pi/180; %transform to radian
rms = 0.025;%rms
f = 5.405; %Frequency (GHz)
rho_b = 1.7;
S = 0.25;
C = 0.25;
mv = [0.2769,0.2769,0.2769,0.2769,0.2769, 0.3477, 0.4617, 0.4805, 0.3386, 0.2324, 0.1779, 0.0751,0.0751,0.0751];
mv = mv*1.7;
T = [14.23, 17.61, 16.21, 21.03, 16.17, 18.77, 22.97, 22.41, 25.85, 31.33, 25.75, 28.19, 30.05, 30.49];
[M,N] = size(mv);
for i = 1:1:N
    [epsr(i), epsi(i)] = RelDielConst_Soil(f,T(i),rho_b, mv(i),S,C);
    eps = epsr(i) - epsi(i)*j;
    [sig_s_vv(i), sig_s_hh, tmp] = PRISM1_ForwardModel(eps,theta_r,rms,f);
end

figure;
subplot(1,4,1);plot(S2_LAI(1,:), mv);datetick('x','mmm');xlabel('Time');ylabel('Volumetric moisture');
subplot(1,4,2);plot(S2_LAI(1,:), T);datetick('x','mmm');xlabel('Time');ylabel('Temperature (degree)');
subplot(1,4,3);
plot(S2_LAI(1,:), epsr);hold on;
plot(S2_LAI(1,:), epsi);
legend('epsilon-r', 'epsilon-i')
datetick('x','mmm');xlabel('Time');ylabel('Permittivity of soil');
subplot(1,4,4);plot(S2_LAI(1,:), sig_s_vv);
datetick('x','mmm');xlabel('Time');ylabel('Bacscattering (dB)');



%% (8)bacscattering value estimation with FC01 time series data
clear
addpath '.\data\'
addpath C:\ADataPreparation
S1VV_FC01_TS = readS1TimeValue('S1VV_asc_2017_142_FC01_meanImgs.csv');
S2_LAI_TS = readTimeValueNew('S2Henshui_2017_lai_FC01.csv');
S2_CW_TS = readTimeValueNew('S2Henshui_2017_cw_FC01.csv');
S2_CW_TS(2,:) = S2_CW_TS(2,:)/10.0;
%% function version of bacscattering value estimation
[eps_v_r, eps_v_i,kappa_s, kappa_a, sigma_0_vv, sigma_0_vv_canopy, sigma_0_vv_gcg, sigma_0_vv_gc]=function_TS_backscattering_estimation(S2_LAI_TS(2,:), S2_CW_TS(2,:));
kappa_e = kappa_s + kappa_a;


%% 
C_water = 10.*S2_LAI_TS(2,:).*S2_CW_TS(2,:);
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
r =  0.00097;% 0.0001;%radius is the disc radius (metres) 1cm for partile radius
% compute real and imaginery part
for i =1:1:size(mg,2)
    [eps_v_r(i), eps_v_i(i)] = RelDielConst_Vegetation(f, mg(i));%4-9.2 Dielectric Model
    [kappa_e(i), kappa_s(i), kappa_a(i)] = function_extinction_scattering_Rayleigh(eps_v_r(i), eps_v_i(i));
    [sigma_0_vv(i), sigma_0_vv_canopy(i), sigma_0_vv_gcg(i), sigma_0_vv_gc(i)] = ...
        S2RTR_DiffuseUB_Kas(eps_v_r(i), eps_v_i(i), h(i), f, s, kappa_a(i), kappa_s(i), theta);
end


%%
figure;
plot(S2_CW_TS(1,:),sigma_0_vv(:), 'm');hold on;
plot(S2_CW_TS(1,:),sigma_0_vv_canopy(:), 'g');hold on;
plot(S2_CW_TS(1,:),sigma_0_vv_gcg(:), 'c');hold on;
plot(S2_CW_TS(1,:),sigma_0_vv_gc(:), 'r');hold on;
plot(S1VV_FC01_TS(1,:),S1VV_FC01_TS(2,:), 'b');hold on;
legend('Estimated VV','Canopy VV', 'Ground canopy ground VV', 'Ground canopy VV', 'S1-VV')
datetick('x','mmm');title('Estimated vv VS S1-VV');
xlabel('Time');ylabel('Backscattering value (dB)')

%%
figure;
subplot(1,2,1);
plot(S2_CW_TS(1,:),kappa_a(:), 'g');hold on;
plot(S2_CW_TS(1,:),kappa_s(:), 'c');hold on;
plot(S2_CW_TS(1,:),kappa_e(:), 'r');hold on;
legend('absorption coefficient', 'scattering coefficient', 'extinction coefficient')
datetick('x','mmm');%title('Estimated vv VS S1-VV');
xlabel('Time');ylabel('Coefficient values')

subplot(1,2,2);
plot(S2_CW_TS(1,:),sigma_0_vv(:), 'm');hold on;
plot(S2_CW_TS(1,:),sigma_0_vv_canopy(:), 'g');hold on;
plot(S2_CW_TS(1,:),sigma_0_vv_gcg(:), 'c');hold on;
plot(S2_CW_TS(1,:),sigma_0_vv_gc(:), 'r');hold on;
plot(S1VV_FC01_TS(1,:),S1VV_FC01_TS(2,:), 'b');hold on;
legend('Estimated VV','Canopy VV', 'Ground canopy ground VV', 'Ground canopy VV', 'S1-VV')
datetick('x','mmm');title('Estimated vv VS S1-VV');
xlabel('Time');ylabel('Backscattering value (dB)')


%% (12)Time series transmitivities computation
clear
addpath '.\data\'
S1VHVV01a = readS1TimeValue('S1_VV_timeSeries_HS01a.csv');
S1VHVV01 = readS1TimeValue('S1_VV_timeSeries_HS01.csv');

S2_LAI = readTimeValueNew('S2Henshui_2017_lai_HS01.csv'); % S2_LAI(2,:)
S2_CW = readTimeValueNew('S2Henshui_2017_cw_HS01.csv');
S2_CW(2,:) = S2_CW(2,:)/10.0;

C_water = 10.*S2_LAI(2,:).*S2_CW(2,:);
mv = C_water * 0.3988;%C_water * leaf thickness * water density
rho = 0.3;
mg = mv./(mv + (1 - mv)*rho);
%
h = [15,  20.4, 49.73, 68, 75, 76.4];% height of the wheat
d = h/100;
%
f = 5.405;
s = 0.025;
a = 0.015;
kappa_e = 1.695;
theta = 34.81;

%%
theta_r = theta *pi/180; %transform to radian
kappa_s = a .* kappa_e ; %scattering coefficient
%-- calculate transmissivity in layer
tau = kappa_e .*d .* sec(theta_r);
T = exp(-tau);
r = T.^2;%[0.5383, 0.4307, 0.1283, 0.0603, 0.0452, 0.0427]




%% (6)Time series canopy backscattering computation with empirical a and k_e
clear
addpath '.\data\'
S1VHVV01a = readS1TimeValue('S1_VV_timeSeries_HS01a.csv');
S1VHVV01 = readS1TimeValue('S1_VV_timeSeries_HS01.csv');

S2_LAI = readTimeValueNew('S2Henshui_2017_lai_HS01.csv'); % S2_LAI(2,:)
S2_CW = readTimeValueNew('S2Henshui_2017_cw_HS01.csv');
S2_CW(2,:) = S2_CW(2,:)/10.0;

C_water = 10.*S2_LAI(2,:).*S2_CW(2,:);
mv = C_water * 0.3988;%C_water * leaf thickness * water density
rho = 0.3;
mg = mv./(mv + (1 - mv)*rho);
%
h = [15,  20.4, 49.73, 68, 75, 76.4];% height of the wheat
h = h/100;
%
f = 5.405;
s = 0.025;
a = 0.015;
kappa_e = 1.695;
theta = 34.81;

% compute real and imaginery part
for i =1:1:size(mg,2)
    [eps_v_r(i), eps_v_i(i)] = RelDielConst_Vegetation(f, mg(i));%4-9.2 Dielectric Model
    [sigma_0_vv(i), sigma_0_vv_canopy(i), sigma_0_vv_gcg(i), sigma_0_vv_gc(i)] = S2RTR_DiffuseUB2_aKe(eps_v_r(i), eps_v_i(i), a, h(i), f, s,  kappa_e,  theta);
%     [sigma_0_vv(i), sigma_0_hh(i), sigma_0_vv_canopy(i)] = S2RTR_DiffuseUB2(eps_v_r(i), eps_v_i(i), a, h(i), f, s, kappa_e, theta);
end

figure;
subplot(2,3,1);plot(S2_LAI(1,:),S2_LAI(2,:));datetick('x','mmm');title('S2-LAI');
subplot(2,3,2);plot(S2_CW(1,:),S2_CW(2,:));datetick('x','mmm');title('S2-CW');
subplot(2,3,3);plot(S2_CW(1,:),h(:));datetick('x','mmm');title('Height of Wheat (m)');

subplot(2,3,4);
plot(S2_CW(1,:),mv(:),'r');hold on;
plot(S2_CW(1,:),mg(:),'b');hold on;datetick('x','mmm');title('Estimated mv (r) and mg (b)');

subplot(2,3,5);
plot(S2_CW(1,:),eps_v_r(:));hold on;
plot(S2_CW(1,:),eps_v_i(:));hold on;datetick('x','mmm');title('Estimated real (r) and imaginary (b) part');

subplot(2,3,6);
plot(S2_CW(1,:),sigma_0_vv_canopy(:), 'g');hold on;
plot(S2_CW(1,:),sigma_0_vv_gcg(:), 'c');hold on;
plot(S2_CW(1,:),sigma_0_vv_gc(:), 'r');hold on;
plot(S1VHVV01(1,:),S1VHVV01(2,:), 'b');hold on;
legend('Canopy VV', 'Ground canopy ground VV', 'Ground canopy VV', 'S1-VV')
datetick('x','mmm');title('Estimated vv VS S1-VV');
xlabel('Time');ylabel('Backscattering value (dB)')


% subplot(2,3,6);plot(S2_CW(1,:),sigma_0_vv(:), 'r');hold on;
% plot(S2_CW(1,:),sigma_0_vv_canopy(:), 'g');hold on;
% plot(S1VHVV01(1,:),S1VHVV01(2,:), 'b');hold on;
% legend('Estimated vv', 'canopy vv', 'S1-VV')
% datetick('x','mmm');title('Estimated vv VS S1-VV');




%% (10)Soil backscattering estimation
clear

theta_r = theta *pi/180; %transform to radian

%-- call the PRISM-1 surface scattering model Section 10-5
[sig_s_vv, sig_s_hh, tmp] = PRISM1_ForwardModel(eps,theta,s,f);
sig_s_vv  = 10.^(sig_s_vv/10); 
sig_s_hh  = 10.^(sig_s_hh/10); 


%% (3)Time series sigma_vv computation over HS01
clear
addpath '.\data\'
addpath C:\ADataPreparation
S1VHVV01a = readS1TimeValue('S1_VV_timeSeries_HS01a.csv');
S1VHVV01 = readS1TimeValue('S1_VV_timeSeries_HS01.csv');

S2_LAI = readTimeValueNew('S2Henshui_2017_lai_HS01.csv'); % S2_LAI(2,:)
S2_CW = readTimeValueNew('S2Henshui_2017_cw_HS01.csv');
S2_CW(2,:) = S2_CW(2,:)/10.0;

C_water = 10.*S2_LAI(2,:).*S2_CW(2,:);
mv = C_water * 0.3988;%C_water * leaf thickness * water density
rho = 0.3;
mg = mv./(mv + (1 - mv)*rho);
%
h = [15,  20.4, 49.73, 68, 75, 76.4];% height of the wheat
h = h/100;
%
f = 5.405;
s = 0.025;
a = 0.015;
kappa_e = 1.695;
theta = 34.81;

% compute real and imaginery part
for i =1:1:size(mg,2)
    [eps_v_r(i), eps_v_i(i)] = RelDielConst_Vegetation(f, mg(i));%4-9.2 Dielectric Model
    [sigma_0_vv(i), sigma_0_hh(i), sigma_0_vv_canopy(i)] = S2RTR_DiffuseUB2(eps_v_r(i), eps_v_i(i), a, h(i), f, s, kappa_e, theta);
end

figure;
subplot(2,3,1);plot(S2_LAI(1,:),S2_LAI(2,:));datetick('x','mmm');title('S2-LAI');
subplot(2,3,2);plot(S2_CW(1,:),S2_CW(2,:));datetick('x','mmm');title('S2-CW');
subplot(2,3,3);plot(S2_CW(1,:),h(:));datetick('x','mmm');title('Height of Wheat (m)');

subplot(2,3,4);
plot(S2_CW(1,:),mv(:),'r');hold on;
plot(S2_CW(1,:),mg(:),'b');hold on;datetick('x','mmm');title('Estimated mv (r) and mg (b)');

subplot(2,3,5);
plot(S2_CW(1,:),eps_v_r(:));hold on;
plot(S2_CW(1,:),eps_v_i(:));hold on;datetick('x','mmm');title('Estimated real (r) and imaginary (b) part');

subplot(2,3,6);plot(S2_CW(1,:),sigma_0_vv(:), 'r');hold on;
plot(S2_CW(1,:),sigma_0_vv_canopy(:), 'g');hold on;
plot(S1VHVV01(1,:),S1VHVV01(2,:), 'b');hold on;
legend('Estimated vv', 'canopy vv', 'S1-VV')
datetick('x','mmm');title('Estimated vv VS S1-VV');



%% (5)Compute canopy backscattering with k_a and k_s of HS01
% clear
% theta = 39.578;
clear
addpath C:\ADataPreparation
addpath '.\data\'
% S1VHVV01a = readS1TimeValue('S1_VV_timeSeries_HS01a.csv');
S1VHVV01 = readS1TimeValue('S1_VV_timeSeries_HS01.csv');

S2_LAI = readTimeValueNew('S2Henshui_2017_lai_HS01.csv'); % S2_LAI(2,:)
S2_CW = readTimeValueNew('S2Henshui_2017_cw_HS01.csv');
S2_CW(2,:) = S2_CW(2,:)/10.0;

C_water = 10.*S2_LAI(2,:).*S2_CW(2,:);
mv = C_water * 0.3988;%C_water * leaf thickness * water density
rho = 0.3;
mg = mv./(mv + (1 - mv)*rho);
%
h = [15,  20.4, 49.73, 68, 75, 76.4];% height of the wheat
h = h/100;
%
f = 5.405; %Frequency (GHz)
s = 0.025; %rms height (m)
% a = 0.015;
% kappa_e = 1.695;
theta = 34.81;

% compute real and imaginery part
for i =1:1:size(mg,2)
    [eps_v_r(i), eps_v_i(i)] = RelDielConst_Vegetation(f, mg(i));%4-9.2 Dielectric Model
    [kappa_e(i), kappa_s(i), kappa_a(i)] = function_extinction_scattering_Rayleigh(eps_v_r(i), eps_v_i(i));
    [sigma_0_vv(i), sigma_0_vv_canopy(i), sigma_0_vv_gcg(i), sigma_0_vv_gc(i)] = ...
        S2RTR_DiffuseUB_Kas(eps_v_r(i), eps_v_i(i), h(i), f, s, kappa_a(i), kappa_s(i), theta);
end
%%
figure;
subplot(1,2,1);
plot(S2_CW(1,:),kappa_a(:), 'g');hold on;
plot(S2_CW(1,:),kappa_s(:), 'c');hold on;
plot(S2_CW(1,:),kappa_e(:), 'r');hold on;
legend('absorption coefficient', 'scattering coefficient', 'extinction coefficient')
datetick('x','mmm');%title('Estimated vv VS S1-VV');
xlabel('Time');ylabel('Coefficient values')

subplot(1,2,2);
plot(S2_CW(1,:),sigma_0_vv_canopy(:), 'g');hold on;
plot(S2_CW(1,:),sigma_0_vv_gcg(:), 'c');hold on;
plot(S2_CW(1,:),sigma_0_vv_gc(:), 'r');hold on;
plot(S1VHVV01(1,:),S1VHVV01(2,:), 'b');hold on;
legend('Canopy VV', 'Ground canopy ground VV', 'Ground canopy VV', 'S1-VV')
datetick('x','mmm');title('Estimated vv VS S1-VV');
xlabel('Time');ylabel('Backscattering value (dB)')

%%
figure;
% plot(S2_CW(1,:),sigma_0_vv(:), 'r');hold on;
plot(S2_CW(1,:),sigma_0_vv_canopy(:), 'g');hold on;
plot(S2_CW(1,:),sigma_0_vv_gcg(:), 'c');hold on;
plot(S2_CW(1,:),sigma_0_vv_gc(:), 'r');hold on;
plot(S1VHVV01(1,:),S1VHVV01(2,:), 'b');hold on;
legend('Canopy VV', 'Ground canopy ground VV', 'Ground canopy VV', 'S1-VV')
datetick('x','mmm');title('Estimated vv VS S1-VV');
xlabel('Time');ylabel('Backscattering value (dB)')




%% (11)Gravimetric water content computation
clear
addpath '.\data\'
addpath C:\ADataPreparation
%% import S2-LAI/CW
% need remove the NaN images
S2_LAI = readTimeValueNew('S2Henshui_2017_lai_HS01.csv'); % S2_LAI(2,:)
S2_CW = readTimeValueNew('S2Henshui_2017_cw_HS01.csv');
S2_CW(2,:) = S2_CW(2,:)/10.0;
%% old method
C_water = 10.*S2_LAI(2,:).*S2_CW(2,:);
mv1 = C_water * 0.3988;%C_water * leaf thickness * water density
rho = 0.3;
mg1 = mv1./(mv1 + (1 - mv1)*rho);
%% new method
mv = S2_CW(2,:)./(0.02*S2_LAI(2,:));
mg = mv./(mv+(1-mv).*0.3);


mg_insitu = [0.53096, 0.76, 0.57517];
mg_time = [736783, 736819, 736846];
% [3/29/2017, 2017/5/4, 5/31/2017];

%%
figure;
subplot(1,3,1);plot(S2_LAI(1,:),S2_LAI(2,:));datetick('x','mmm');title('S2-LAI');
subplot(1,3,2);plot(S2_CW(1,:),S2_CW(2,:));datetick('x','mmm');title('S2-CW');
subplot(1,3,3);
plot(S2_CW(1,:),mv(:),'r');hold on;
plot(S2_CW(1,:),mg(:),'b');hold on;

plot(S2_CW(1,:),mv1(:),'r--');hold on;
plot(S2_CW(1,:),mg1(:),'b--');hold on;datetick('x','mmm');title('Estimated mv (r) and mg (b)');
scatter(mg_time, mg_insitu, '*');
legend('new m_v','new m_g','old m_v','old m_g','insitu m_g')
figure;









%% (1)Import LAI and CW images
clear
addpath '.\data\'
addpath C:\ADataPreparation
[S2_LAI,longrid1,latgrid1,Infor1] = geoimread('S2Henshui_2017_lai_HS01.tif');
[S2_CW,longrid1,latgrid1,Infor1] = geoimread('S2Henshui_2017_cw_HS01.tif');

S2_LAI_0428 = double(S2_LAI(:,:,10))./1000;
S2_CW_0428 = double(S2_CW(:,:,10))./10000;
% figure;
% subplot(1,2,1);imagesc(S2_LAI_0428);colorbar
% subplot(1,2,2);imagesc(S2_CW_0428);colorbar

C_water = 10.*S2_LAI_0428.*S2_CW_0428;

h = 0.7;% height of the wheat
mv = C_water * 0.3988;%C_water * leaf thickness * water density
rho = 0.3;
mg = mv./(mv + (1 - mv)*rho);
angle = 34.81;
% figure;
% subplot(2,2,1);imagesc(S2_LAI_0428);colorbar;title('S2-LAI')
% subplot(2,2,2);imagesc(S2_CW_0428);colorbar;title('S2-CW')
% subplot(2,2,3);imagesc(mv);colorbar;title('volume moisture')
% subplot(2,2,4);imagesc(mg);colorbar;title('gravimetric moisture')

% figure;
% subplot(1,2,1);imagesc(mv);colorbar;title('volume moisture')
% subplot(1,2,2);imagesc(mg);colorbar;title('gravimetric moisture')

% compute real and imaginery part
[M,N] = size(mg);
for i = 1:1:M
    for j = 1:1:N
        if mg(i,j) > 0
            [eps_v_r(i,j), eps_v_i(i,j)] = RelDielConst_Vegetation(5.405, mg(i,j));
            [sigma_0_vv(i,j), sigma_0_hh(i,j)] = S2RTR_DiffuseUB1(eps_v_r(i,j), eps_v_i(i,j), h);
        end 
%         [eps_v_r, eps_v_i] = RelDielConst_Vegetation(5.405,0.5);
    end
end
sigma_0_vv((sigma_0_vv==0)) = nan;
sigma_0_hh((sigma_0_hh==0)) = nan;
% figure;
% subplot(2,2,1);imagesc(eps_v_r);colorbar;title('Real Part')
% subplot(2,2,2);imagesc(eps_v_i);colorbar;title('Imaginary Part')
% subplot(2,2,3);imagesc(sigma_0_vv);colorbar;title('sigma-vv')
% subplot(2,2,4);imagesc(sigma_0_hh);colorbar;title('sigma-hh')

[S1_VV_0426_HS01,longrid1,latgrid1,Infor1] = geoimread('collectionVV_noisy.tif');
[S1_VV_0426_HS01,longrid1,latgrid1,Infor1] = geoimread('S1_VV_0426_HS01.tif');

figure;
subplot(2,4,1);imagesc(S2_LAI_0428);colorbar;title('(1)S2-LAI');set(gca,'XTick',[], 'YTick', [])
subplot(2,4,2);imagesc(S2_CW_0428);colorbar;title('(2)S2-CW');set(gca,'XTick',[], 'YTick', [])
subplot(2,4,3);imagesc(mv);colorbar;title('(3)Volume moisture');set(gca,'XTick',[], 'YTick', [])
subplot(2,4,4);imagesc(mg);colorbar;title('(4)Gravimetric moisture');set(gca,'XTick',[], 'YTick', [])

subplot(2,4,5);imagesc(eps_v_r);colorbar;title('(5)Real Part');set(gca,'XTick',[], 'YTick', [])
subplot(2,4,6);imagesc(eps_v_i);colorbar;title('(6)Imaginary Part');set(gca,'XTick',[], 'YTick', [])
subplot(2,4,7);imagesc(sigma_0_vv);colorbar;title('(7)Canopy related sigma-vv');set(gca,'XTick',[], 'YTick', [])
subplot(2,4,8);imagesc(S1_VV_0426_HS01);colorbar;title('(8)S1-VV-0426-denoised');set(gca,'XTick',[], 'YTick', [])




%% (6)Time series canopy backscattering computation with empirical a and k_e
clear
addpath '.\data\'
addpath C:\ADataPreparation
S1VHVV01 = readS1TimeValue('S1VV_asc_2017_142_FC01_meanImgs.csv');
%% import S2-LAI/CW
% need remove the NaN images
S2_LAI = readTimeValueNew('S2Henshui_2017_lai_FC01.csv');
S2_CW = readTimeValueNew('S2Henshui_2017_cw_FC01.csv');
S2_CW(2,:) = S2_CW(2,:)/10.0;

C_water = 10.*S2_LAI(2,:).*S2_CW(2,:);
mv = C_water * 0.3988;%C_water * leaf thickness * water density
rho = 0.3;
mg = mv./(mv + (1 - mv)*rho);
%
h =[10,15,16,17,26,30,35.0326,58.6863,71.4251,82.6688,82.1204,81.1737,16,16];% height of the wheat
h = h/100;
%
f = 5.405; %Frequency (GHz)
s = 0.025; %rms height (m)
a = 0.015;
kappa_e = 1.695;
theta = 39.5925;
%

% compute real and imaginery part
for i =1:1:size(mg,2)
    [eps_v_r(i), eps_v_i(i)] = RelDielConst_Vegetation(f, mg(i));%4-9.2 Dielectric Model
    [sigma_0_vv(i), sigma_0_vv_canopy(i), sigma_0_vv_gcg(i), sigma_0_vv_gc(i)] = S2RTR_DiffuseUB2_aKe(eps_v_r(i), eps_v_i(i), a, h(i), f, s,  kappa_e,  theta);
%     [sigma_0_vv(i), sigma_0_hh(i), sigma_0_vv_canopy(i)] = S2RTR_DiffuseUB2(eps_v_r(i), eps_v_i(i), a, h(i), f, s, kappa_e, theta);
end

figure;
subplot(2,3,1);plot(S2_LAI(1,:),S2_LAI(2,:));datetick('x','mmm');title('S2-LAI');
subplot(2,3,2);plot(S2_CW(1,:),S2_CW(2,:));datetick('x','mmm');title('S2-CW');
subplot(2,3,3);plot(S2_CW(1,:),h(:));datetick('x','mmm');title('Height of Wheat (m)');

subplot(2,3,4);
plot(S2_CW(1,:),mv(:),'r');hold on;
plot(S2_CW(1,:),mg(:),'b');hold on;datetick('x','mmm');title('Estimated mv (r) and mg (b)');

subplot(2,3,5);
plot(S2_CW(1,:),eps_v_r(:));hold on;
plot(S2_CW(1,:),eps_v_i(:));hold on;datetick('x','mmm');title('Estimated real (r) and imaginary (b) part');

subplot(2,3,6);
plot(S2_CW(1,:),sigma_0_vv_canopy(:), 'g');hold on;
plot(S2_CW(1,:),sigma_0_vv_gcg(:), 'c');hold on;
plot(S2_CW(1,:),sigma_0_vv_gc(:), 'r');hold on;
plot(S1VHVV01(1,:),S1VHVV01(2,:), 'b');hold on;
legend('Canopy VV', 'Ground canopy ground VV', 'Ground canopy VV', 'S1-VV')
datetick('x','mmm');title('Estimated vv VS S1-VV');
xlabel('Time');ylabel('Backscattering value (dB)')


%% (7)bacscattering value estimation with FC01 time series images
clear
addpath '.\data\'
addpath C:\ADataPreparation
%% read S1 denoised data with masks
[S1VV_FC01,longrid1,latgrid1,Infor1] = geoimread('S1VV_asc_2017_142_FC01_meanImgs.tif');
[S1VV_FC01_w,longrid1,latgrid1,Infor1] = geoimread('S1VV_asc_2017_142_FC01_whole.tif');
S1VV_FC01_TS = readS1TimeValue('S1VV_asc_2017_142_FC01_meanImgs.csv');
S1VV_FC01_img = S1VV_FC01(:,:,2:11);
S1VV_FC01_w_Img = S1VV_FC01_w(:,:,2:11);
%% import S2-LAI/CW
% need remove the NaN images
[S2LAI_FC01,longrid1,latgrid1,Infor1] = geoimread('S2Henshui_2017_lai_FC01.tif');
[S2CW_FC01,longrid1,latgrid1,Infor1] = geoimread('S2Henshui_2017_cw_FC01.tif');
S2_LAI_TS = readTimeValueNew('S2Henshui_2017_lai_FC01.csv');
S2_CW_TS = readTimeValueNew('S2Henshui_2017_cw_FC01.csv');

[S2cbrownTime S2cbrownTSerie] = csvimport('S2Henshui_2017_lai_FC01.csv','columns', [1, 2],'noHeader', true);%mean VV test area
[M,N]=size(S2cbrownTime);
k=1;
for i=1:1:M    
    if abs(S2cbrownTSerie{i}) > 0
        data(1,k)= datenum(S2cbrownTime{i},'mm/dd/yyyy');
        data(2,k)= str2num(S2cbrownTSerie{i})/1000.0;
        S2LAI_FC01_img(:,:,k) = double(S2LAI_FC01(:,:,i))/1000.0;
        S2CW_FC01_img(:,:,k) = double(S2CW_FC01(:,:,i))/10000.0;
        k = k+1;        
    end
end

%%
[M,N,Num] = size(S2CW_FC01_img);
for i = 1:1:M
    for j = 1:1:N
        if sum(S2LAI_FC01_img(i,j,:)) > 0
            [eps_v_r, eps_v_i,kappa_s, kappa_a, sigma_0_vv, sigma_0_vv_canopy, sigma_0_vv_gcg, sigma_0_vv_gc]...
                =function_TS_backscattering_estimation(squeeze(S2LAI_FC01_img(i,j,:))', squeeze(S2CW_FC01_img(i,j,:))');
            eps_v_r_w(i,j,:) = eps_v_r;
            eps_v_i_w(i,j,:) = eps_v_i;
            kappa_s_w(i,j,:) = kappa_s;
            kappa_a_w(i,j,:) = kappa_a;
            sigma_0_vv_w(i,j,:) = sigma_0_vv;
            sigma_0_vv_canopy_w(i,j,:) = sigma_0_vv_canopy;
            sigma_0_vv_gcg_w(i,j,:) = sigma_0_vv_gcg;
            sigma_0_vv_gc_w(i,j,:) = sigma_0_vv_gc;
        end
    end
end
%%
sigma_0_vv_canopy_w((sigma_0_vv_canopy_w==0))=nan;
sigma_0_vv_gcg_w((sigma_0_vv_gcg_w==0))=nan;
sigma_0_vv_gc_w((sigma_0_vv_gc_w==0))=nan;
figure;
subplot(2,2,1);imagesc(sigma_0_vv_canopy_w(:,:,10));
subplot(2,2,2);imagesc(sigma_0_vv_gcg_w(:,:,10));
subplot(2,2,3);imagesc(sigma_0_vv_gc_w(:,:,10));
subplot(2,2,4);imagesc(S1VV_FC01_w(:,:,9));


%%
figure;
histogram(sigma_0_vv_canopy_w(:,:,10));hold on;
histogram(sigma_0_vv_gcg_w(:,:,10));hold on;
histogram(sigma_0_vv_gc_w(:,:,10));hold on;
histogram(S1VV_FC01_w(:,:,9));hold on;
axis([-25 -2 0 3500])
legend('canopy','ground canopy ground','ground canopy','S1-VV')
set(gca,'XTick',[], 'YTick', [])

%%
S1_mean = mean(S1VV_FC01_img,3);
gcg_w = sigma_0_vv_gcg_w(:,:,10);
gcg_w((S1_mean==nan))=nan;
figure;imagesc(gcg_w);









%% (2)Interpolate
clear
%%
x = [736785, 736821, 736846];
v = [23.8, 78.2, 79.8];
xq = 736785:1:736846;
vq = interp1(x,v,xq,'spline');

vq1 = interp1(x,v,xq);



%%
x = [30, 66, 93];
v = [20.4, 73.2, 76.4];
xq = 30:1:93;
vq = interp1(x,v,xq,'spline');
vq1 = interp1(x,v,xq);
figure;
plot(xq,vq);hold on;
plot(xq,vq1);hold on;



%% (4)Extinction and scattering in a Rayleigh medium
clear
%% function with input eps_v_r(i,j), eps_v_i(i,j)
eps_v_r = 21.7411;
eps_v_i = 7.6965;

[k_e, k_s, k_a] = function_extinction_scattering_Rayleigh(eps_v_r, eps_v_i);

%% single point
v = 0.005;%the volume fraction of scatterers
% v = (4/3*pi*r^3) * Nv;
lambda= 0.056;%lambda is the wavelength (in metres)
k = 2*pi/lambda;%
r = 0.01;%radius is the disc radius (metres) 1cm for partile radius
epsilon = 21.7411 + 7.6965j;
K1 = (epsilon - 1)/(epsilon + 2);

k_s = 2*v*k^4*r^3*(abs(K1))^2;
k_a = v*k*imag(epsilon)*(abs(3/(epsilon+2)))^2;

k_e = k_a + k_s;






