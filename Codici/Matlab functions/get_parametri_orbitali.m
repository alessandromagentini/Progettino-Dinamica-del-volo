function [obj_param] = get_parametri_orbitali(r,v, mu)

% Function che calcola i parametri orbitali classici da stato r, v.
%
%   Questa funzione determina la tipologia di orbita (circolare, ellittica, 
%   parabolica o iperbolica) e calcola i principali parametri geometrici 
%   e orientativi dell'orbita.
%
%   INPUT:
%       r        - Vettore posizione nel sistema ECI [1x3] o [3x1]  [km]
%       v        - Vettore velocità nel sistema ECI [1x3] o [3x1]   [km/s]
%       mu       - Parametro gravitazionale standard del corpo      [km^3/s^2]
%
%   OUTPUT:
%       obj_param - Struct contenente i parametri calcolati:
%           .r0_vec - raggio iniziale                               [km]
%           .v0_vec - velocità iniziale                             [km/s]
%           .h_vec  - Vettore momento angolare specifico            [km^2/s]
%           .E      - Energia specifica                             [km^2/s^2]
%           .a      - Semiasse maggiore                             [m]
%           .rp, .ra- Raggio al pericentro e apocentro              [m]
%           .e_vec  - vettore eccentricità                          [-]
%           .e      - Eccentricità scalare                          [-]
%           .T      - Periodo orbitale (solo per orbite chiuse)     [s]
%           .i      - Inclinazione                                  [deg]
%           .raan   - Ascensione retta del nodo ascendente (Omega)  [deg]
%           .omega  - Argomento di pericentro                       [deg]
%           .mu     - Parametro gravitazionale standard del corpo   [km^3/s^2]   
%           .type   - Stringa descrittiva del tipo di orbita        [-]


% k - versore
k =[0 0 1];

% h - momento angolare specifico (vettore)
h = cross(r,v);                                                            %[km^2/s]

% E - Energia specifica
E = (norm(v)^2)/2 - mu/norm(r);                                            %[km^2/s^2]

% Identifica forma dell'orbita tramite E e calcola i parametri orbitali
if E < 0    % orbita ellittica/circolare
    % a - semiasse maggiore
    a = -mu/(2*E);                                                         %[km]

    % e - eccentricità
    e_vec = (cross(v,h)./mu - r./norm(r))';
    e = norm(e_vec);                                                       % e = sqrt(1+(2*norm(h)^2*E)/mu^2);

    if e == 0
        type = "circular";
    else
        type = "elliptical";
    end

    % rp - raggio al pericentro  
    rp = a*(1 - e)*1000;                                                   %[m]
    % ra - raggio all'apogeo
    ra = a*(1 + e)*1000;                                                   %[m]

    % T - periodo orbitale
    T = (2*pi/sqrt(mu))*a^(3/2);                                           %[s]

    % N - Vettore della linea dei nodi
    N = cross(k,h)';

    % i - inclinazione
    i = acos(h(3)/norm(h));                                                %[rad]
    i = rad2deg(i);                                                        %[deg]
    % opzionale: aggiungere identificazione orbita prograda/retrograda

    % raan - right ascension of ascending node
    if N(2) > 0
        raan = acos(N(1)/norm(N));                                         %[rad]
    else
        raan = 2*pi - acos(N(1)/norm(N));                                  %[rad]
    end
    raan = rad2deg(raan);                                                  %[deg]    

    % omega - argomento di pericentro
    cosomega = dot(N,e_vec)/(norm(N)*e); 
    cosomega = max(-1, min(1, cosomega));                                  % forza in [-1, 1] per evitare errori numerici
    if e_vec(3) >= 0
        omega = acos(cosomega);                                            %[rad]
    else
        omega = 2*pi - acos(cosomega);                                     %[rad]
    end
    omega = rad2deg(omega);                                                %[deg] 
            
    % Ricompatta tutto in uno struct
    obj_param = struct("r0_vec",r, "v0_vec",v, "h_vec",h, "E",E, ...
                       "a",a*1000, "rp",rp, "ra",ra, "e_vec",e_vec, "e",e, ...
                       "T",T, "i",i, "raan",raan, "omega",omega, ...
                       "mu",mu, "type",type);

elseif E == 0   % orbita parabolica
    type = "parabolic";

    % e - eccentricità
    e_vec = (cross(v,h)./mu - r./norm(r))';
    e = 1;

    % rp - raggio di pericentro
    rp = norm(h)^2/(2*mu);
    % da completare

    obj_param = struct("r0_vec",r, "v0_vec",v, "h_vec",h, "E",E, ...
        "rp",rp, "e_vec",e_vec, "e",e,"mu",mu, "type",type);

else  % orbita iperbolica
    type = "hyperbolic";

    % da completare
    obj_param = NaN;
end








