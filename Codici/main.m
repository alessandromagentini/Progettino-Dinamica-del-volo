clear;clc;close all

% Aggiungo il path della cartella delle funciton
addpath("Matlab functions")

% Flags
groundtrack3_flag        = 0;      % per geoplot 3D della ground track
grundtrack2_flag         = 0;      % per geoplot 2D della ground track
plot_eci_flag            = 1;      % per plot (non globe) in ECI
toolbox_flag             = 0;      % per utilizzo satellite communication toolbox

%% Dati iniziali
r0_vec = [-7368.038574853538, -7231.584293256432, -148.523707822187];             %[Km]
v0_vec = [4.126512186315761, -3.956371322777358, -0.490613661500991];             %[Km/s]

mu_terra = 398600;                                                                %[km^3/s^2]

start_time   = datetime(2001, 9, 11, 12, 00, 00);                                 % YYYY-MM-DD-HH-min-sec
stop_time    = datetime(2001, 9, 12, 15, 00, 00);                                 % YYYY-MM-DD-HH-min-sec
delta_t_sat_sample = 3*3600;                                                      %[s]

%% Calcolo parametri orbitali
[sat_param] = get_parametri_orbitali(r0_vec,v0_vec,mu_terra);


%% Calcolo orbita 
%1) Custom
dt = 1;                                                                           %[s]
[sat_orbit] = propagatore(sat_param,dt,delta_t_sat_sample,start_time,stop_time);  %% Da restituire valori t_sample
% PLOT
plotter(sat_param,sat_orbit,grundtrack2_flag,groundtrack3_flag,plot_eci_flag)


%2) Satellite Communications Toolbox
if toolbox_flag == 1
    deltat_sample = 1;        %[s]
    sc = satelliteScenario(start_time,stop_time,deltat_sample);
    sat_orbit_tb = satellite(sc,sat_param.a,sat_param.e,sat_param.i,sat_param.raan,sat_param.omega,0);

    % PLOT
    v = satelliteScenarioViewer(sc);
end




