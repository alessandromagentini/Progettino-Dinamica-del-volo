function [dydt,data] = integratore(t,y, data)
% Integratore per ode45, per propagazione orbita con perturbazione J2. 
% I calcoli vengono fatti già nel sistema ECI

% Setup costanti
J2 = data.J2;
RT = data.raggio_terra;      % raggio Terra [m]
mu = data.mu;                % [m^3/s^2]

J2_flag = data.flags.J2_flag;
moon_flag = data.flags.moon_flag;
sun_flag = data.flags.sun_flag;

% Organizzazione dati input
r = y(1:3);          % XYZ
drdt = y(4:6);       % XYZ

% Calcolo accelerazione
g = -mu / norm(r)^3 * r;

X = r(1); Y = r(2); Z = r(3);

if J2_flag == 1
    a_J2 = -3/2 * J2 * mu * RT^2 / norm(r)^5 * ...   % Perturbazione J2
        [X*(1 - 5*Z^2/norm(r)^2);
        Y*(1 - 5*Z^2/norm(r)^2);
        Z*(3 - 5*Z^2/norm(r)^2)];
else
    a_J2 = 0;
end

if moon_flag == 1
    a_tb_moon = third_body(r,t,data,"Moon");  % Perturbazione di terzo corpo (Luna)
else
    a_tb_moon = 0;
end

if sun_flag == 1
    a_tb_sun = third_body(r,t,data,"Sun");  % Perturbazione di terzo corpo (Sole)
else
    a_tb_sun = 0;
end
% Calcolo parametri orbitali
par =  get_parametri_orbitali(r*1e-3,drdt*1e-3,mu);
data.i     = par.i;
data.raan  = par.raan;
data.omega = par.omega;


dvdt = g + a_J2 + a_tb_moon + a_tb_sun;
dydt = [drdt; dvdt];

end



