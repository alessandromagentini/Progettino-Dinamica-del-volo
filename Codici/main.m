clear;clc;close all

addpath("Matlab functions")

%% Dati iniziali
r0_vec = [-7368.038574853538, -7231.584293256432, -148.523707822187];      %[km]
v0_vec = [4.126512186315761, -3.956371322777358, -0.490613661500991];      %[km/s]

mu_terra = 398600;                                                         %[km^3/s^2]

start_time  = datetime(2026, 3, 20, 12, 00, 00);                           % YYYY-MM-DD-HH-min-sec
end_time    = datetime(2026, 3, 21, 12, 00, 00);                           % YYYY-MM-DD-HH-min-sec
t_sat_sample = (3*(0:8))*3600;                                             %[s]

%% Calcolo parametri orbitali
[sat_param] = get_parametri_orbitali(r0_vec,v0_vec,mu_terra);


%% Calcolo orbita 
% Custom
[sat_orbit] = propagatore(sat_param,1,t_sat_sample); %% Da restituire valori t_sample











% % Aerospace Toolbox
% deltat_sample = 3*360;        %[s]
% r0_vec_t = (r0_vec*100)';     %[m/s]
% v0_vec_t = (v0_vec*100)';     %[m/s]
% time = start_time:seconds(deltat_sample):end_time;
% numericalPropOpts = Aero.spacecraft.NumericalPropagatorOptions( ...
%       ODESet = odeset(RelTol=1e-8,AbsTol=1e-8,MaxStep=300), ...
%                IncludeThirdBodyGravity=false, ...
%                ThirdBodyGravitySource=["Sun" "Moon"], ...
%                GravitationalPotentialModel="point-mass");
% [r,v] = propagateOrbit( ...
%       time, ...
%       r0_vec_t, ...
%       v0_vec_t, ...
%       NumericalPropagatorOptions=numericalPropOpts);

% % PLOT AL VOLO
% globe_uif=uifigure;
% globe_uif.Pointer='crosshair';
% for i = 1:length(sat_orbit.r_tool)
%     time_step = start_time + seconds(sat_orbit.t(i));
%     lla(i,:) = eci2lla(sat_orbit.r_tool(i,:), datevec(time_step));
% end
% r_terra = 6371000;
% lat = lla(:,1);
% lon = lla(:,2);
% alt = r_terra + lla(:,3);
% gg = geoglobe(globe_uif);
% geoplot3(gg, lat, lon, alt, 'r', 'LineWidth', 2);












