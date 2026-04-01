clear;clc;close all

% Aggiungo il path della cartella delle funciton
addpath("Matlab functions")

% Flags
groundtrack3_flag        = 0;      % per geoplot 3D della ground track
groundtrack2_flag        = 0;      % per geoplot 2D della ground track
plot_eci_flag            = 0;      % per plot (non globe) in ECI
satellite_tb_flag        = 0;      % per utilizzo satellite communication toolbox
simulink_flag            = 0;      % per plot del modello simulink
analisi_risultati_flag   = 1;      % per verifica dei risulati

%% Dati iniziali
r0_vec = [-7368.038574853538, -7231.584293256432, -148.523707822187];             %[Km]
v0_vec = [4.126512186315761, -3.956371322777358, -0.490613661500991];             %[Km/s]

mu_terra = 398600.4418;                                                           %[km^3/s^2]

start_time   = datetime(2026, 3, 26, 18, 00, 00);                                 % YYYY-MM-DD-HH-min-sec
stop_time    = datetime(2026, 3, 27,  6, 00, 00);                                 % YYYY-MM-DD-HH-min-sec
mission_duration = seconds(stop_time - start_time);                               %[s]
delta_t_sat_sample = 30 * 60;                                                     %[s]

%% Calcolo parametri orbitali
fprintf("Calcolo dei parametri orbitali...")
[sat_param] = get_parametri_orbitali(r0_vec,v0_vec,mu_terra);
fprintf(" completato!\n")

%% Calcolo orbita 
%1) Custom
dt = 1;                                                                           %[s]
fprintf("Propagazione dell'orbita con function custom...")
[sat_orbit] = propagatore(sat_param,dt,delta_t_sat_sample,start_time,stop_time); 
fprintf(" completato!\n")
if groundtrack3_flag == 1 || groundtrack2_flag == 1 || plot_eci_flag == 1  % PLOT
    plotter(sat_param,sat_orbit,groundtrack3_flag,groundtrack2_flag,plot_eci_flag)
end
res.sat_orbit = sat_orbit;

%2) Satellite Communications Toolbox
if satellite_tb_flag == 1
    deltat_sample = 60;        %[s]
    sc_sat_tb = satelliteScenario(start_time,stop_time,deltat_sample);
    sat_orbit_sat_tb = satellite(sc_sat_tb,sat_param.a,sat_param.e,sat_param.i,sat_param.raan,sat_param.omega,sat_param.TA0);

    % PLOT
    v = satelliteScenarioViewer(sc_sat_tb);
end

%3) Simulink model with aerospace toolbox
% Simulazione
fprintf("Propagazione dell'orbita con simulink e aerospace toolbox...")
simOut = sim("modello_simulink");
fprintf(" completato!\n")

% Salvo i dati
posData_icrf  = simOut.yout{1}.Values;
velData_icrf  = simOut.yout{2}.Values;
timeData      = simOut.yout{3}.Values;

sc_aero_tb = satelliteScenario(start_time, stop_time, 60);
sat_orbit_aero_tb = satellite(sc_aero_tb, posData_icrf, velData_icrf);

%% Salvo i saples richiesti                                                  
samples_date = (start_time:seconds(delta_t_sat_sample):stop_time)';
n_samples = length(samples_date);
r_samples_tb = zeros(n_samples, 3);
t_start_sim = timeData.Time(1); 
for k = 1:n_samples
    tempo_cercato = t_start_sim + (k - 1) * delta_t_sat_sample;  
    [scostamento, idx] = min(abs(timeData.Time - tempo_cercato));
    if scostamento < 1e-3
        r_samples_tb(k,:) = posData_icrf.Data(idx, :);                                   
    else 
        r_samples_tb(k,:) = interp1(timeData.Time, posData_icrf.Data, tempo_cercato, 'linear');
    end
end

res.aerotb_res = struct("pos_icrf",posData_icrf, "vel_icrf",velData_icrf, "time",timeData);

if simulink_flag == 1 % PLOT
    % groundTrack(sat_orbit_aero_tb);
    play(sc_aero_tb);
end

if analisi_risultati_flag == 1
    analisi_risultati(sat_param,res)
end
