function [r_eci, v_eci] = perifocale2ECI(r_peri, v_peri, i, raan, omega)

i     = deg2rad(i);
raan  = deg2rad(raan);
omega = deg2rad(omega);

% Per conversione in ECI
    R_raan = [cos(raan)  -sin(raan)  0;
              sin(raan)   cos(raan)  0;
              0            0           1];
    R_i     = [1   0       0;
               0   cos(i)  -sin(i);
               0   sin(i)   cos(i)];
    R_omega = [cos(omega)  -sin(omega)  0;
               sin(omega)   cos(omega)  0;
               0            0           1];
    R =  R_omega * R_i * R_raan;

    r_eci = R * r_peri;
    v_eci = R * v_peri;
end

% Mancano la matrice di precessione dell'asse terrestre ecc. (da implementare se
% c'è tempo, ma l'errore nei nostri intervalli temporali dovrebbe essere
% invisibile)