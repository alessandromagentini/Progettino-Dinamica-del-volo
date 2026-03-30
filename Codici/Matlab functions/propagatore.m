function [orbit_polar] = propagatore(orbit_param, dt, t_sample, startTime, stopTime)
% Funzione che calcola l'orbita a partire dai parametri orbitali nel
% sistema perifocale e poi restituisce l'orbita in ECI

%% Preparazione variabili
a         = orbit_param.a;
T         = orbit_param.T;
e         = orbit_param.e;
h         = norm(orbit_param.h_vec);
mu        = orbit_param.mu;
i         = orbit_param.i;
raan      = orbit_param.raan;
omega     = orbit_param.omega;
TA0       = deg2rad(orbit_param.TA0);
xi0       = deg2rad(orbit_param.xi0);
M_e0      = deg2rad(orbit_param.M_e0);
t0        = orbit_param.t0;

%% Funzioni
% r
r_fun = @(TA) (norm(h)^2/mu)/(1 + e*cos(TA));
% True anomaly
TA_fun = @(xi) 2*atan(sqrt((1 + e)/(1 - e)) * tan(xi/2));

n = round(seconds(stopTime - startTime)/dt);
% n = round(T/dt) + 1;  % numero intervalli dt su un periodo T dell'orbita

%% Inizializzazione
xi        = zeros(n+1,1); xi(1) = xi0;
TA        = zeros(n+1,1); TA(1) = TA0;
r_kepl    = zeros(n+1,3);         
v_kepl    = zeros(n+1,3);
t         = zeros(n+1,1);
fpa       = zeros(n+1,1);
v_ort     = zeros(n+1,1);
v_rad     = zeros(n+1,1);
M_e       = zeros(n+1,1); M_e(1) = M_e0;

r_eci    = zeros(n+1,3);
v_eci    = zeros(n+1,3);

%% Risolutore:
% t -> M_e -> xi -> TA
for idx = 1:(n+1)                         %eventualmente si porebbe mettere un parfor
    t(idx) = t0 + (idx * dt - dt);                                         % tempo trascorso in [s]
    M_e(idx) = 2*pi/T * t(idx);                                            % Anomalia media
    xi(idx) = KeplerE(M_e(idx),e);                                         % Anomalia Eccentrica
    TA(idx) = TA_fun(xi(idx));                                             % Calcolo True Anomaly                                          % conversione in [deg]
    r_kepl(idx)  = r_fun(TA(idx));                                         % Calcolo del raggio in [m]
    fpa(idx) = atan((e*sin(TA(idx)))/(1 + e*cos(TA(idx))));                % Flight path angle
    
    v_ort(idx)  = mu/h * (1 + e*cos(TA(idx)));                             % Velocità ortogonale
    v_rad(idx)  = mu/h * e * sin(TA(idx));                                 % Velocità radiale
    v_kepl(idx) = norm([v_ort(idx) v_rad(idx)]);                           % Velocità risultante

    %% Conversione in ECI:
    %% DA FARE?) scrivere la nostra funzione
    if TA(idx) < 0
        TA(idx) = 2*pi + TA(idx);
    end
    [rt, vt] = keplerian2ijk(a, e, i, raan, omega, rad2deg(TA(idx))); %funzione dell'aerospace toolbox temporanea/per verificare
    r_eci(idx,:) = rt';
    v_eci(idx,:) = vt';
end
t_date = (startTime:seconds(dt):stopTime)';                                % date istanti temporali 

%% Salvo i dati richiesti                                                  
samples_date = (startTime:seconds(t_sample):stopTime)';
n_samples = length(samples_date);
r_samples = zeros(n_samples, 3);
t_start_sim = t(1); 
for k = 1:n_samples
    tempo_cercato = t_start_sim + (k - 1) * t_sample;  
    [scostamento, idx] = min(abs(t - tempo_cercato));
    if scostamento < 1e-3
        r_samples(k,:) = r_eci(idx, :);                                   
    else 
        r_samples(k,:) = interp1(t, r_eci, tempo_cercato, 'linear');
    end
end
%% Organizzazione output
samples = struct("date",samples_date, "position",r_samples);
orbit_polar = struct("r_kepl",r_kepl, "v_ort",v_ort, "v_rad",v_rad, "v_kepl",v_kepl, ...
                     "TA",rad2deg(TA), "t",t,"t_date",t_date,"dt",dt, "M_e",rad2deg(M_e), "xi",rad2deg(xi), "r_eci",r_eci, "v_eci",v_eci, "fpa",fpa, "samples",samples);

end %propagatore



% Function locali
function E = KeplerE(M_e,e)
    
    fun = @(x) x - e*sin(x) - M_e; 

    % find best x0
    if M_e <= pi
        x0 = M_e + e/2;
    else 
        x0 = M_e - e/2;
    end
    % opzionale: usare formula più complessa data a lezione

    E = fzero(fun, x0);          % fzero per risolvere in E
    % Opzionale: usare N-R scritto esplicitamente
end