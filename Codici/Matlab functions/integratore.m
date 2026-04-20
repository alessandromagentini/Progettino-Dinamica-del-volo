function dydt = integratore(t,y)
% Integratore per ode45, per propagazione orbita con perturbazione J2. 
% I calcoli vengono fatti già nel sistema ECI

% Setup costanti
J2 = 1.08263e-3;
RT = 6371e3;        % raggio Terra [m]
mu_terra = 3.986004418e14;  %[m^3/s^2]

% Organizzazione dati input
r = y(1:3);          % XYZ
drdt = y(4:6);       % XYZ

% Calcoli
g = -mu_terra / norm(r)^3 * r;
X = r(1); Y = r(2); Z = r(3); % conversione dalla formula del prof: Z/r = sin(latitudine)
a_J2 = -3/2 * J2 * mu_terra * RT^2 / norm(r)^5 * ...   % Perturbazione J2
       [X*(1 - 5*Z^2/norm(r)^2);
        Y*(1 - 5*Z^2/norm(r)^2);
        Z*(3 - 5*Z^2/norm(r)^2)];

dvdt = g + a_J2;
dydt = [drdt; dvdt];

end