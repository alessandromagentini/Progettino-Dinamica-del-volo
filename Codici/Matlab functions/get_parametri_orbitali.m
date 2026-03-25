function [obj_param] = get_parametri_orbitali(r,v, mu)

% k - versore
k =[0 0 1];

% h - momento angolare specifico (vettore)
h = cross(r,v);

% E - Energia specifica
E = (norm(v)^2)/2 - mu/norm(r);

% Identifica forma dell'orbita tramite E
if E < 0
    type = "elliptical";

    % a - semiasse maggiore
    a = -mu/(2*E);

    % e - eccentricità
    e_vec = (cross(v,h)./mu - r./norm(r))';
    e = norm(e_vec);                                                       % e = sqrt(1+(2*norm(h)^2*E)/mu^2);

    % T - periodo orbitale
    T = (2*pi/sqrt(mu))*a^(3/2);                                           % [s]

    % N - Vettore della linea dei nodi
    N = cross(k,h)';

    % i - inclinazione
    if h(3)>0
    i = acos(h(3)/norm(h));                                                % [rad]
    else 
    i = pi - acos(h(3)/norm(h));  %controllare formula appunti             % [rad]
    end
    i = rad2deg(i);                                                        % [deg]
    % opzionale: aggiungere identificazione orbita prograda/retrograda


    % raan - right ascension of ascending node
    if N(2) > 0
        raan = acos(N(1)/norm(N));                                         % [rad]
    else
        raan = 2*pi - acos(N(1)/norm(N));                                  % [rad]
    end
    raan = rad2deg(raan);                                                  % [deg]    

    % omega - argomento di pericentro
    cosomega = dot(N,e_vec)/(norm(N)*e); 
    cosomega = max(-1, min(1, cosomega));                                            % forza in [-1, 1] per evitare errori numerici
    if e_vec(3) > 0
        omega = acos(cosomega);                            % [rad]
    else
        omega = 2*pi - acos(cosomega);                     % [rad]
    end
    omega = rad2deg(omega);                                                % [deg] 
            
    % Ricompatta tutto in uno struct
    obj_param = struct("r0_vec",r, "v0_vec",v, "h_vec",h, "E",E, ...
                       "a",a, "e_vec",e_vec, "e",e, "T",T, ...
                       "i",i, "raan",raan, "omega",omega, ...
                       "mu",mu, "type",type);

elseif E == 0
    type = "parabolic";

    % da completare

else 
    type = "hyperbolic";

    % da completare

end








