function a_tb = third_body(r,t, data, third_body)
% Function per calcolare l'accelereazione dovuta alla perturbazione di
% effetti di terzo corpo con l'uso di function dell'aerospace toolbox

tb = third_body;
dateTime = data.startTime + seconds(t);
jd = juliandate(dateTime);

if tb == "Moon"
    mu_tb = 4.9048695e12;   % [m³/s²]
elseif tb == "Sun" 
    mu_tb  = 1.32712440e20;  % [m³/s²]
else
    error("thiird body non supportato")
end

r_tb = planetEphemeris(jd, 'Earth', tb) * 1e3;   % [m]
r_tb = r_tb';

d = r_tb - r;
a_tb = mu_tb * (d/norm(d)^3 - r_tb/(norm(r_tb)^3));
