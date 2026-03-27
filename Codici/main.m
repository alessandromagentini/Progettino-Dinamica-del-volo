clear;clc;close all

% Aggiungo il path della cartella delle funciton
addpath("Matlab functions")

% Flags
groundtrack3_flag        = 1;      % per geoplot 3D della ground track
groundtrack2_flag        = 0;      % per geoplot 2D della ground track
plot_eci_flag            = 0;      % per plot (non globe) in ECI
satellite_tb_flag        = 0;      % per utilizzo satellite communication toolbox
simulink_flag            = 1;      % per utilizzo del modello simulink
confronta_risultati_flag = 1;      % per verifica dei risulati

%% Dati iniziali
r0_vec = [-7368.038574853538, -7231.584293256432, -148.523707822187];             %[Km]
v0_vec = [4.126512186315761, -3.956371322777358, -0.490613661500991];             %[Km/s]

mu_terra = 398600;                                                                %[km^3/s^2]

start_time   = datetime(2001, 9, 11, 12, 00, 00);                                 % YYYY-MM-DD-HH-min-sec
stop_time    = datetime(2001, 9, 11, 15, 00, 00);                                 % YYYY-MM-DD-HH-min-sec
delta_t_sat_sample = 3*3600;                                                      %[s]

%% Calcolo parametri orbitali
[sat_param] = get_parametri_orbitali(r0_vec,v0_vec,mu_terra);

%% Calcolo orbita 
%1) Custom
if groundtrack3_flag == 1 || groundtrack2_flag ==1
    dt = 1;                                                                           %[s]
    [sat_orbit] = propagatore(sat_param,dt,delta_t_sat_sample,start_time,stop_time);  %% Da restituire valori t_sample

    % PLOT
    plotter(sat_param,sat_orbit,groundtrack3_flag,groundtrack2_flag,plot_eci_flag)
    res.sat_orbit = sat_orbit;
else 
    sat_orbit.TA = 0;
end

%2) Satellite Communications Toolbox
if satellite_tb_flag == 1
    deltat_sample = 1;        %[s]
    sc_sat_tb = satelliteScenario(start_time,stop_time,deltat_sample);
    sat_orbit_sat_tb = satellite(sc_sat_tb,sat_param.a,sat_param.e,sat_param.i,sat_param.raan,sat_param.omega,0);

    % PLOT
    v = satelliteScenarioViewer(sc_sat_tb);
end

%3) Simulink model with aerospace toolbox
if simulink_flag ==1
    mission_duration = seconds(stop_time - start_time);
    
    % Simulazione
    simOut = sim("modello_simulink");

    % Salvo i dati
    posData  = simOut.yout{1}.Values;
    velData  = simOut.yout{2}.Values;
    timeData = simOut.yout{3}.Values;

    sc_aero_tb = satelliteScenario(start_time, stop_time, 1);
    sat_orbit_aero_tb = satellite(sc_aero_tb, posData, velData, "CoordinateFrame", "ecef");
    
    % PLOT
    groundTrack(sat_orbit_aero_tb);
    play(sc_aero_tb);

    res.aerotb_res = struct("pos",posData,"vel",velData,"time",timeData);
end

if confronta_risultati_flag == 1
    confronta_risultati(res)
end
