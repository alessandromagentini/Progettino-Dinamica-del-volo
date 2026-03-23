clear;clc;close all

addpath("Matlab functions")

%% Dati iniziali
r_vec = [-7368.038574853538, -7231.584293256432, -148.523707822187];    %[km]
v_vec = [4.126512186315761, -3.956371322777358, -0.490613661500991];    %[km/s]

mu_terra = 398600;                                                      %[km^3/s^2]

t_sat_sample = 3*(0:8);                                                 %[s]

%% Calcolo parametri orbitali
[sat_param] = get_parametri_orbitali(r_vec,v_vec,mu_terra);

%% Calcolo orbita
[sat_orbit] = propagatore(sat_param,1,t_sat_sample); %% Da restituire valori t_sample













