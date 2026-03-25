function [orbit_polar] = propagatore(orbit_param, dt, t_sample,startTime,stopTime)
% Funzione che calcola l'orbita a partire dai parametri orbitali nel
% sistema perifocale e poi restituisce l'orbita in ECI

%% Da restituire valori t_sample %%

%% Preparazione variabili
a         = orbit_param.a;
T         = orbit_param.T;
e         = orbit_param.e;
h         = norm(orbit_param.h_vec);
mu        = orbit_param.mu;
i         = orbit_param.i;
raan      = orbit_param.raan;
omega     = orbit_param.omega;

%% Funzioni
% r
r_fun = @(TA) (norm(h)^2/mu)/(1 + e*cos(TA));
% True anomaly
TA_fun = @(xi) 2*atan(sqrt((1 + e)/(1 - e)) * tan(xi/2));

n = round(seconds(stopTime - startTime)/dt);
% n = round(T/dt) + 1;  % numero intervalli dt


%% Calcolo M_e - anomalia media per ogni dt
M_e = ((2*pi/T) .* dt*(0:n))'; 

%% Inizializzazione
xi        = zeros(n+1,1);
TA        = zeros(n+1,1);
r_kepl    = zeros(n+1,3);
v_kepl    = zeros(n+1,3);
t         = zeros(n+1,1);

r_eci    = zeros(n+1,3);
v_eci    = zeros(n+1,3);

%% Risolutore:
% t -> M_e -> xi -> TA
for idx = 1:(n+1)
    xi(idx) = KeplerE(M_e(idx),e);           % Risoluzione equazione di Keplero
    TA(idx) = TA_fun(xi(idx));               % Calcolo True Anomaly
    TA(idx) = rad2deg(TA(idx));              % conversione in [deg]
    r_kepl(idx)  = r_fun(TA(idx));           % Calcolo del raggio
    t(idx)  = idx * dt - dt;                 % tempo trascorso

    %% da aggiungere: calcolo della velocità

    %% Conversione in ECI:
    %% DA FARE?) scrivere la nostra funzione
    if TA(idx) < 0
        TA(idx) = 360 + TA(idx);
    end
    [rt, vt] = keplerian2ijk(a, e, i, raan, omega, TA(idx)); %funzione dell'aerospace toolbox temporanea/per verificare
    r_eci(idx,:) = rt';
    v_eci(idx,:) = vt';
end

%% Organizzazione output
orbit_polar = struct("r_kepl",r_kepl, "v_kepl",v_kepl, "TA",TA, "t",t, "M_e",M_e, "xi",xi, "r_eci",r_eci, "v_eci",v_eci);

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