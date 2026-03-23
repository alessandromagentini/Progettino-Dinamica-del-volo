function [orbit_polar] = propagatore(orbit_param, dt, t_sample)

%% Da restituire valori t_sample %%

%% Preparazione variabili
T   = orbit_param.T;
e   = orbit_param.e;
h   = norm(orbit_param.h_vec);
mu  = orbit_param.mu;

%% Funzioni
% r
r_fun = @(TA) (norm(h)^2/mu)/(1 + e*cos(TA));
% True anomaly
TA_fun = @(xi) 2*atan(sqrt((1 + e)/(1 - e)) * tan(xi/2));

n = round(T/dt) + 1;  % numero intervalli dt

%% Calcolo M_e - anomalia media per ogni dt
M_e = (2*pi/T) .* dt*(0:n); 

%% Inizializzazione
xi   = zeros(n+1,1);
TA   = zeros(n+1,1);
r    = zeros(n+1,1);
t    = zeros(n+1,1);

%% Risolutore:
% t -> M_e -> xi -> TA
for i = 1:(n+1)
    xi(i) = KeplerE(M_e(i),e);
    TA(i) = TA_fun(xi(i));          % Calculate true anomaly
    r(i) = r_fun(TA(i));            % Calculate radius using the defined function
    t(i) = i * dt - dt;                  % Calculate time at each sample
end

%% Organizzazione output
orbit_polar = struct("r",r,"TA",TA,"t",t,"M_e",M_e,"xi",xi);

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