clear;clc;close all

addpath("Matlab functions")

%% Dati iniziali
r0_vec = [-7368.038574853538, -7231.584293256432, -148.523707822187];            %[km]
v0_vec = [4.126512186315761, -3.956371322777358, -0.490613661500991];            %[km/s]

mu_terra = 398600;                                                               %[km^3/s^2]

start_time   = datetime(2001, 9, 11, 12, 00, 00);                                % YYYY-MM-DD-HH-min-sec
stop_time    = datetime(2001, 9, 12, 15, 00, 00);                                % YYYY-MM-DD-HH-min-sec
delta_t_sat_sample = 3*3600;                                                     %[s]

%% Calcolo parametri orbitali
[sat_param] = get_parametri_orbitali(r0_vec,v0_vec,mu_terra);


%% Calcolo orbita 
% Custom
dt = 1;                                                                           %[s]
[sat_orbit] = propagatore(sat_param,dt,delta_t_sat_sample,start_time,stop_time);  %% Da restituire valori t_sample
% PLOT
t_utc = start_time:seconds(dt):stop_time;
lla = eci2lla(sat_orbit.r_eci*1000, datevec(t_utc));
lat = lla(:,1);
lon = lla(:,2);
alt = lla(:,3);
fig = uifigure;
globe = geoglobe(fig);
geoplot3(globe,lat, lon, alt, 'LineWidth', 2, 'Color', 'r')

% Satellite Communications Toolbox
deltat_sample = 1;        %[s]
sc = satelliteScenario(start_time,stop_time,deltat_sample);
sat_orbit_tb = satellite(sc,sat_param.a*1000,sat_param.e,sat_param.i,sat_param.raan,sat_param.omega,0);
% PLOT
v = satelliteScenarioViewer(sc);














