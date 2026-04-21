function [orbit_polar] = propagatore(orbit_param, dt, t_sample, startTime, stopTime, flags)
%   Questa funzione propaga l'orbita di un satellite a partire dai parametri
%   orbitali iniziali. Supporta sia una propagazione analitica (Kepleriana) 
%   sia una propagazione numerica (ODE45) che include perturbazioni J2 
%   ed effetti di terzo corpo (Sole/Luna).
%
%   INPUT:
%       orbit_param : Struct contenente i parametri orbitali e costanti:
%                     .a, .e, .i, .raan, .omega, .TA0 [deg], .mu, .T, .h_vec,
%                     .M_e0 [deg], .xi0 [deg], .t0, .r0_vec, .v0_vec
%       dt          : Passo temporale fisso per il modello kepleriano [s]
%       t_sample    : Passo di campionamento per i dati di output filtrati [s]
%       startTime   : Datetime di inizio simulazione
%       stopTime    : Datetime di fine simulazione
%       flags       : Struct di controllo:
%                     flags utilizzate:
%                         .modello_propagatore : "keplerian" o "numerical"
%                         .moon_flag           : 1 attiva perturbazione Luna
%                         .sun_flag            : 1 attiva perturbazione Sole
%
%   OUTPUT:
%   orbit_polar : Struct contenente i risultati della propagazione:
%                .r_eci, .v_eci    : Posizione e velocità in ECI [m, m/s]
%                .t_date           : Vettore datetime degli istanti calcolati
%                .samples          : Struct con posizione campionata a t_sample
%                .i, .raan, .omega : Evoluzione dei parametri (se numerico)
%                .r_kepl, .v_kepl  : pos e vel in perifocale (solo kepleriano)
%                .TA, .M_e, .xi    : Anomalie (Vera, Media, Eccentrica) [deg]

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
r0        = orbit_param.r0_vec;
v0        = orbit_param.v0_vec;

J2 = 1.08263*1e-3;
RT = 6371*1e3; %[m]

%% Funzioni
% r
r_fun = @(TA) (norm(h)^2/mu)/(1 + e*cos(TA));
% True anomaly
TA_fun = @(xi) 2*atan(sqrt((1 + e)/(1 - e)) * tan(xi/2));

n = round(seconds(stopTime - startTime)/dt);
% n = round(T/dt) + 1;  % numero intervalli dt su un periodo T dell'orbita

%% Inizializzazione
xi          = zeros(n+1,1); xi(1) = xi0;
TA          = zeros(n+1,1); TA(1) = TA0;
r_kepl      = zeros(n+1,1);         
v_kepl      = zeros(n+1,1);
t           = zeros(n+1,1);
fpa         = zeros(n+1,1);
v_ort       = zeros(n+1,1);
v_rad       = zeros(n+1,1);
M_e         = zeros(n+1,1); M_e(1) = M_e0;

r_kepl_vec  = zeros(n+1,3);
v_kepl_vec  = zeros(n+1,3);
r_eci       = zeros(n+1,3);
v_eci       = zeros(n+1,3);

%% Risolutore:
% 1) Kepleriano 
type = flags.modello_propagatore;
if type == "keplerian"   % t -> M_e -> xi -> TA
    for idx = 1:(n+1)                         %eventualmente si porebbe mettere un parfor
        t(idx) = t0 + (idx * dt - dt);                                         % tempo trascorso in [s]
        M_e(idx) = 2*pi/T * t(idx);                                            % Anomalia media
        xi(idx) = KeplerE(M_e(idx),e);                                         % Anomalia Eccentrica
        TA(idx) = TA_fun(xi(idx));                                             % Calcolo True Anomaly                          
        r_kepl(idx)  = r_fun(TA(idx));                                         % Calcolo del raggio in [m]
        fpa(idx) = atan((e*sin(TA(idx)))/(1 + e*cos(TA(idx))));                % Flight path angle
        
        v_ort(idx)  = mu/h * (1 + e*cos(TA(idx)));                             % Velocità ortogonale
        v_rad(idx)  = mu/h * e * sin(TA(idx));                                 % Velocità radiale
        v_kepl(idx) = norm([v_ort(idx) v_rad(idx)]);                           % Velocità risultante
    
        if TA(idx) < 0
            TA(idx) = 2*pi + TA(idx);
        end
    
        %% Conversione in ECI:
        [rt, vt] = keplerian2ijk(a, e, i, raan, omega, rad2deg(TA(idx))); %funzione dell'aerospace toolbox temporanea/per verificare
        r_eci(idx,:) = rt';
        v_eci(idx,:) = vt';
    
        % % custom
        % r_kepl_vec(idx,:) = [r_kepl(idx)*cos(TA(idx)); r_kepl(idx)*sin(TA(idx)); 0];
        % v_kepl_vec(idx,:) = (mu/h) * [-sin(TA(idx)); e + cos(TA(idx)); 0];
        % [r_eci(idx,:), v_eci(idx,:)] = perifocale2ECI((r_kepl_vec(idx,:)*1000)',(v_kepl_vec(idx,:)*1000)',i,raan,omega);
    
        %[r_eci(idx,:), v_eci(idx,:)] = orb2eci(mu,[a,e,deg2rad(i),deg2rad(raan),deg2rad(omega),TA(idx)]);
    end
    t_date = (startTime:seconds(dt):stopTime)';  % date istanti temporali 
    
% 2) Numerico
elseif type == "numerical"
    r0 = r0'*1e3; v0 = v0'*1e3; %converto in [m]
    tol = 1e-13;
    max_step = 600;  %[s]
    t_span = [0, seconds(stopTime - startTime)];
    
    % Anticipo creando una tabella per gli effetti di terzo corpo (poi andrà interpolata)
    if flags.moon_flag == 1
        TB.moon_table = get_thirdBody_table(startTime,t_span(2), 5400, "Moon"); % 16 valori/giorno
    else
        TB.moon_table = NaN;
    end

    if flags.sun_flag == 1
        TB.sun_table = get_thirdBody_table(startTime,t_span(2), 10800, "Sun"); % 8 valore/giorno
    else
        TB.sun_table = NaN;
    end
    data = struct("orbit_param",orbit_param,"r0",r0,"v0",v0,"startTime",startTime,"t_span",t_span, ...
        "mu",mu*1e9,"J2",J2,"raggio_terra",RT, "TB",TB, "flags",flags);

    % INTEGRAZIONE
    int_opt = get_integration_opt(tol,max_step); % aggiungere function per manovra
    [t, yout] = ode45(@(t,y) integratore(t,y,data), t_span, [r0; v0], int_opt);
    r_eci = yout(:, 1:3);
    v_eci = yout(:, 4:6);

    t_date = startTime + seconds(t); % date istanti temporali 

    % Ricostruisco i parametri orbitali per vederne l'andamento
    i      = zeros(length(t),1);
    raan   = zeros(length(t),1);
    omega  = zeros(length(t),1);
    for idx = 1:length(t)
        par = get_parametri_orbitali(r_eci(idx,:)*1e-3,v_eci(idx,:)*1e-3, mu);
        i(idx)     = par.i;
        raan(idx)  = par.raan;
        omega(idx) = par.omega;
    end

else
    error("type must be numerical or keplerian")
end

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
                     "TA",rad2deg(TA), "t",t,"t_date",t_date,"dt",dt, "M_e",rad2deg(M_e), "xi",rad2deg(xi), ...
                     "r_eci",r_eci, "v_eci",v_eci, "fpa",fpa, ...
                     "i",i, "raan",raan, "omega", omega, ...
                     "samples",samples);

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

