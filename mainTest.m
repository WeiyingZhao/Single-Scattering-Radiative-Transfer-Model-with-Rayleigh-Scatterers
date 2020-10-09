%% bacscattering value estimation with FC01 time series data
clear
addpath '.\data\'
addpath C:\ADataPreparation
S1VV_FC01_TS = readS1TimeValue('S1VV_asc_2017_142_FC01_meanImgs.csv');
S2_LAI_TS = readTimeValueNew('S2_2017_lai_FC01.csv');
S2_CW_TS = readTimeValueNew('S2_2017_cw_FC01.csv');
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

