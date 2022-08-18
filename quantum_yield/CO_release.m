%% import irradiance data for the lamp; this is given at 0.5 m and 500 W (units mW/m^2*nm)
irradiance = importdata('66142.dat');
% correct to 300 W by multiplying with 3/5
irradiance(:,2) = irradiance(:,2)*3/5;


%% import transmittance data for 3 filters round wavelengths to nearest whole number 
transmittance = importdata('transmittance_filters.csv');
T_data = transmittance.data(); 
T_data(:,2:2:6) = T_data(:,2:2:6)./100;
T_data(:,1:2:5) = round(T_data(:,1:2:5));
T_280_cutoff = flipud(T_data(1:2:end-1,5:6));

%% truncate irradiance data 200 - 800 nm and correct irradiance for a given filter
filter_low = 200;
filter_high = 800;
irradiance = irradiance(irradiance(:,1)>=filter_low & irradiance(:,1)<=filter_high, :);
corrected_irradiance_280_cutoff(:,1) = irradiance(:,1);
corrected_irradiance_280_cutoff(:,2) = T_280_cutoff(:,2).* irradiance(:,2);

%% correct irradiance for distance from the reactor different from 0.5 m for initial data
distance = 0.3;
corrected_irradiance_280_cutoff(:,2) = corrected_irradiance_280_cutoff(:,2)*0.5^2/distance^2;


%% diameter of the round bottom flask in cm
diam = 4.1; 

% depth of the solution in cm
depth = 1.5;

% calculate the area irradiated by the lamp in m^2
area = pi*diam/2*(depth)*1e-4;

% calculate the power per wavelength
power_spectrum_280_cutoff = corrected_irradiance_280_cutoff(:,2) * area;


%% import extinction coefficients for the substrate
% wavelenght in nm versus extinction coefficient 
e_ox = importdata('oxetanone_extinction_2.txt');
e_ox = flipud(e_ox(:,[1 5]));
e_cb = importdata('cyclobutanone_extinction_2.txt');
e_cb = flipud(e_cb(:,[1 5]));
e_az = importdata('azetidinone_extinction_2.txt');
e_az = flipud(e_az(:,[1 5]));
e_th = importdata('thientanone_extinction_2.txt');
e_th = flipud(e_th(:,[1 5]));

%%
% Get wavelengths and extinctions coeffs
% The wavelengths upper limit is 326nm which has an index of 64.
upper_lim_idx = 64;
w_ox = e_ox(1:2:end,1);
w_ox = w_ox(1:upper_lim_idx);
f_ox = e_ox(1:2:end,2);
f_ox = f_ox(1:upper_lim_idx);

w_cb = e_cb(1:2:end,1);
w_cb = w_cb(1:upper_lim_idx);
f_cb = e_cb(1:2:end,2);
f_cb = f_cb(1:upper_lim_idx);

w_az = e_az(1:2:end,1);
w_az = w_az(1:upper_lim_idx);
f_az = e_az(1:2:end,2);
f_az = f_az(1:upper_lim_idx);

w_th = e_th(1:2:end,1);
w_th = w_th(1:upper_lim_idx);
f_th = e_th(1:2:end,2);
f_th = f_th(1:upper_lim_idx);
% change to W from mW
power_spectrum_280_cutoff = power_spectrum_280_cutoff/1000;
power_spectrum_280_cutoff = power_spectrum_280_cutoff(1:upper_lim_idx);

%% enter concentrations of solutions
concs = [1.6e-1 8e-2 4e-2 2e-2 1e-2];


%% get the integral between power spectrum and extinction coefficients

theta_max = acos(1-2*depth/diam);
irradiated_solution_area = integral2(@(phi, theta) (diam/2)^2.*sin(theta),0,pi,0,theta_max);
integ_ox = [];
integ_cb = [];
integ_az = [];
for i = 1:size(concs,2)
    integ_ox_path = [];
    integ_cb_path = [];
    integ_az_path = [];
    for j = 1:size(f_ox,1)
        integ_ox_path = cat(1,integ_ox_path,integral2(@(phi, theta) (1-10.^(-f_ox(j)*concs(i)*diam*sin(theta).*sin(phi))).*sin(theta)*(diam/2)^2,0,pi,0,theta_max));
        integ_cb_path = cat(1,integ_cb_path,integral2(@(phi, theta) (1-10.^(-f_cb(j)*concs(i)*diam*sin(theta).*sin(phi))).*sin(theta)*(diam/2)^2,0,pi,0,theta_max));
        integ_az_path = cat(1,integ_az_path,integral2(@(phi, theta) (1-10.^(-f_az(j)*concs(i)*diam*sin(theta).*sin(phi))).*sin(theta)*(diam/2)^2,0,pi,0,theta_max));
    end
    integ_ox = cat(1,integ_ox, trapz(w_ox*1e-9,power_spectrum_280_cutoff.*w_ox.*integ_ox_path));
    integ_cb = cat(1,integ_cb, trapz(w_cb*1e-9,power_spectrum_280_cutoff.*w_cb.*integ_cb_path));
    integ_az = cat(1,integ_az, trapz(w_az*1e-9,power_spectrum_280_cutoff.*w_az.*integ_az_path));
end

integ_ox = (1/irradiated_solution_area)*integ_ox;
integ_cb = (1/irradiated_solution_area)*integ_cb;
integ_az = (1/irradiated_solution_area)*integ_az;

%% Calculate rate of absorbed photons, I_a from three experiments
% V = 7e-6 m^3 volume of the reaction cell
% h = 6.63e-34 Js Planck's constant
% c = 3e8 m/s speed of light
% N_a = 6.022e23 Avogadro's constant

V = 7e-6;
h = 6.63e-34;
c = 3e8;
N_a = 6.022e23;
I_a_ox = (1./(1e3*N_a*V*h*c)).*integ_ox;
I_a_cb = (1./(1e3*N_a*V*h*c)).*integ_cb;
I_a_az = (1./(1e3*N_a*V*h*c)).*integ_az;


%% initial rate for oxetanone
ox_data = importdata('CO_release_oxetanone.csv');
ox_data = ox_data.data;
t = ox_data(:,3);
CO_ppm_10 = ox_data(:,4)*1e-6*28e-3/(22.4);
CO_ppm_20 = ox_data(~isnan(ox_data(:,5)),5)*1e-6*28e-3/(22.4);
CO_ppm_40 = ox_data(~isnan(ox_data(:,6)),6)*1e-6*28e-3/(22.4);
CO_ppm_80 = ox_data(~isnan(ox_data(:,7)),7)*1e-6*28e-3/(22.4);
CO_ppm_160 = ox_data(~isnan(ox_data(:,8)),8)*1e-6*28e-3/(22.4);
rate_10_ox = [1 0]*lsqr([t ones(size(CO_ppm_10,1),1)], CO_ppm_10);
rate_20_ox = [1 0]*lsqr([t(1:size(CO_ppm_20,1)) ones(size(CO_ppm_20,1),1)], CO_ppm_20);
rate_40_ox = [1 0]*lsqr([t(1:size(CO_ppm_40,1)) ones(size(CO_ppm_40,1),1)], CO_ppm_40);
rate_80_ox = [1 0]*lsqr([t(1:size(CO_ppm_80,1)) ones(size(CO_ppm_80,1),1)], CO_ppm_80);
rate_160_ox = [1 0]*lsqr([t(1:size(CO_ppm_160,1)) ones(size(CO_ppm_160,1),1)], CO_ppm_160);

%% initial rates for cyclobutanone
cb_data = importdata('CO_release_cyclobutanone.csv');
cb_data = cb_data.data;
t = cb_data(:,3);
CO_ppm_10 = cb_data(:,4)*1e-6*28e-3/(22.4);
CO_ppm_20 = cb_data(~isnan(cb_data(:,5)),5)*1e-6*28e-3/(22.4);
CO_ppm_40 = cb_data(~isnan(cb_data(:,6)),6)*1e-6*28e-3/(22.4);
CO_ppm_80 = cb_data(~isnan(cb_data(:,7)),7)*1e-6*28e-3/(22.4);
CO_ppm_160 = cb_data(~isnan(cb_data(:,8)),8)*1e-6*28e-3/(22.4);
rate_10_cb = [1 0]*lsqr([t ones(size(CO_ppm_10,1),1)], CO_ppm_10);
rate_20_cb = [1 0]*lsqr([t(1:size(CO_ppm_20,1)) ones(size(CO_ppm_20,1),1)], CO_ppm_20);
rate_40_cb = [1 0]*lsqr([t(1:size(CO_ppm_40,1)) ones(size(CO_ppm_40,1),1)], CO_ppm_40);
rate_80_cb = [1 0]*lsqr([t(1:size(CO_ppm_80,1)) ones(size(CO_ppm_80,1),1)], CO_ppm_80);
rate_160_cb = [1 0]*lsqr([t(1:size(CO_ppm_160,1)) ones(size(CO_ppm_160,1),1)], CO_ppm_160);

%% initial rates for azetidinone
az_data = importdata('CO_release_azetidinone.csv');
az_data = az_data.data;
t = az_data(:,3);
CO_ppm_10 = az_data(:,4)*1e-6*28e-3/(22.4);
CO_ppm_20 = az_data(~isnan(az_data(:,5)),5)*1e-6*28e-3/(22.4);
CO_ppm_40 = az_data(~isnan(az_data(:,6)),6)*1e-6*28e-3/(22.4);
CO_ppm_80 = az_data(~isnan(az_data(:,7)),7)*1e-6*28e-3/(22.4);
CO_ppm_160 = az_data(~isnan(az_data(:,8)),8)*1e-6*28e-3/(22.4);
rate_10_az = [1 0]*lsqr([t ones(size(CO_ppm_10,1),1)], CO_ppm_10);
rate_20_az = [1 0]*lsqr([t(1:size(CO_ppm_20,1)) ones(size(CO_ppm_20,1),1)], CO_ppm_20);
rate_40_az = [1 0]*lsqr([t(1:size(CO_ppm_40,1)) ones(size(CO_ppm_40,1),1)], CO_ppm_40);
rate_80_az = [1 0]*lsqr([t(1:size(CO_ppm_80,1)) ones(size(CO_ppm_80,1),1)], CO_ppm_80);
rate_160_az = [1 0]*lsqr([t(1:size(CO_ppm_160,1)) ones(size(CO_ppm_160,1),1)], CO_ppm_160);


%% Initial rates from three experiments
rates_oxetanone = flipud(1/(V*1e3)*[rate_10_ox rate_20_ox rate_40_ox rate_80_ox rate_160_ox].');
rates_azetidinone = flipud(1/(V*1e3)*[rate_10_az rate_20_az rate_40_az rate_80_az rate_160_az].');
rates_cyclobutanone = flipud(1/(V*1e3)*[rate_10_cb rate_20_cb rate_40_cb rate_80_cb rate_160_cb].');



%%
subplot(3,1,1)
plot(w_ox, f_ox, 'b-', 'LineWidth',2)
hold
plot(w_az, f_az, 'r-', 'LineWidth',2)
plot(w_cb, f_cb, 'k-', 'LineWidth',2)
plot(w_th, f_th, 'g-', 'LineWidth', 2)
xlim([200 400])
xlabel('wavelength [nm]')
ylabel('extinction coefficient [M^{-1}cm^{-1}]')
legend('3-oxetanone', '3-azetidinone', 'cyclobutanone', '3-thietanone')

subplot(3,1,2)
plot(concs, rates_oxetanone, 'bo')
hold
plot(concs, rates_azetidinone, 'ro')
plot(concs, rates_cyclobutanone, 'ko')
xlabel('concentration [M]')
ylabel('initial rate [M s^{-1}]')
xlim([0 0.17])

subplot(3,1,3)
c = polyfit(I_a_ox, rates_oxetanone,1);
y_est = polyval(c,I_a_ox)
plot(I_a_ox, rates_oxetanone, 'b*')
hold
plot(I_a_ox,y_est,'b-','LineWidth',2)
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])

plot(I_a_az, rates_azetidinone, 'r*')
c = polyfit(I_a_az, rates_azetidinone,1);
y_est = polyval(c,I_a_az)
plot(I_a_az,y_est,'r-','LineWidth',2)
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])


plot(I_a_cb, rates_cyclobutanone, 'k*')
c = polyfit(I_a_cb, rates_cyclobutanone,1);
y_est = polyval(c,I_a_cb)
plot(I_a_cb,y_est,'k-','LineWidth',2)
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])

xlabel('photon absorption rate [M s^{-1}]')
ylabel('initial rate [M s^{-1}]')

