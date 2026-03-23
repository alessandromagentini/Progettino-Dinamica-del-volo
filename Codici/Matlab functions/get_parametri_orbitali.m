function [obj_param] = get_parametri_orbitali(r,v, mu)

% h - momento angolare specifico
h = cross(r,v);

% E - Energia specifica
E = (norm(v)^2)/2 - mu/norm(r);

% Identifica forma dell'orbita tramite E
if E < 0
    type = "elliptical";

    % a - semiasse maggiore
    a = -mu/(2*E);

    % e - eccentricità
    e_vec = cross(v,h)./mu - r./norm(r);
    e = norm(e_vec);                            % e = sqrt(1+(2*norm(h)^2*E)/mu^2);

    % T - periodo orbitale
    T = (2*pi/sqrt(mu))*a^(3/2);                % [s]

    % Ricompatta tutto in uno struct
    obj_param = struct("r_vec",r,"v_vec",v,"h_vec",h,"E",E,"a",a,"e_vec",e_vec,"e",e,"T",T,"mu",mu,"type",type);

elseif E == 0
    type = "parabolic";

    % da completare

else 
    type = "hyperbolic";

    % da completare

end








