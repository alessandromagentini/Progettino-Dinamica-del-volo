function [] = analisi_risultati(sat_param,res)
% Function per l'analisi e plot dei risultati

mu_terra = 398600.4418;   

% Estrazione dati
custom_data = res.sat_orbit;
aero_toolbox = res.aerotb_res;

start_time = custom_data.t_date(1);

cr = custom_data.r_eci; 
cv = custom_data.v_eci;
ct = custom_data.t - custom_data.t(1);     

tbr = aero_toolbox.pos_icrf.Data(3:end,:); % i primi step dell'integratore hanno passo troppo vicino a 0
tbv = aero_toolbox.vel_icrf.Data(3:end,:); 
tbt = seconds(datetime(aero_toolbox.time.Data(3:end,:)) - start_time);

tbr_interp = interp1(tbt, tbr, ct, "spline");
tbv_interp = interp1(tbt, tbv, ct, "spline");

%% Estrazione dati GMAT
addpath("..\Codici\GMAT")
if exist("Kep_param.csv","file")
    gmat_data = readtable("Kep_param.csv");%#ok
    t_gmat    = gmat_data.SAT_ElapsedSecs;
    e_gmat    = gmat_data.SAT_Earth_ECC;
    RA_gmat   = gmat_data.SAT_EarthICRF_RA;
    TA_gmat   = gmat_data.SAT_Earth_TA;
    raan_gmat = gmat_data.SAT_EarthICRF_RAAN;
    i_gmat    = gmat_data.SAT_EarthICRF_INC;
    aop_gmat  = gmat_data.SAT_EarthICRF_AOP;

    % Limito l'intervallo gmat allo stesso dei calcoli fatti in matlab
    [~, idx_gmat_end] = min(abs(t_gmat - ct(end)));
    t_gmat = t_gmat(1:idx_gmat_end);
    i_gmat = i_gmat(1:idx_gmat_end);
    raan_gmat = raan_gmat(1:idx_gmat_end);
    aop_gmat = aop_gmat(1:idx_gmat_end);
    % Associo i valori trami

    i_gmat_interp        = interp1(t_gmat,i_gmat,ct,"spline");
    raan_gmat_interp     = interp1(t_gmat,raan_gmat,ct,"spline");
    aop_gmat_interp      = interp1(t_gmat,aop_gmat,ct,"spline");
end

% Calcolo degli scarti tra toolbox e custom
delta_pos = sqrt(sum((cr - tbr_interp).^2, 2));
delta_vel = sqrt(sum((cv - tbv_interp).^2, 2));

scartor_percentuale = (delta_pos ./ norm(tbr_interp)) * 100;
scartov_percentuale = (delta_vel ./ norm(tbv_interp)) * 100;

% Andamento dei parametri orbitali toolbox
i = zeros(length(tbr_interp),1); raan  = zeros(length(tbr_interp),1);
T = zeros(length(tbr_interp),1); omega = zeros(length(tbr_interp),1);
for idx = 1:length(tbr_interp)
    param = get_parametri_orbitali(tbr_interp(idx,:)/1000,tbv_interp(idx,:)/1000, mu_terra);
    i(idx) = param.i;
    raan(idx) = param.raan;
    omega(idx) = param.omega;
    T(idx) = param.T;
end

%% --- PLOT ---

% Plot scostamenti tra i modelli
figure('Name', 'Confronto Propagatori', 'NumberTitle', 'off')
% Subplot 1: Scarto Assoluto posizione
subplot(2,2,1)
plot(ct/3600, delta_pos/1000, 'LineWidth', 1.5, 'Color', "r")
hold on
pr = polyfit(ct, delta_pos/1000, 1); 
r_reg = polyval(pr, ct);
plot(ct/3600, r_reg, 'g--', 'LineWidth', 2)
grid on
ylabel('Scostamento Posizione [km]')
xlabel('Tempo [ore]')
title('Differenza Assoluta')

% Subplot 2: Scarto Percentuale posizione
subplot(2,2,2)
plot(ct/3600, scartor_percentuale, 'LineWidth', 1.5, 'Color', "b")
hold on
percr = polyfit(ct, scartor_percentuale, 1); 
r_perc_reg = polyval(percr, ct);
plot(ct/3600, r_perc_reg, 'g--', 'LineWidth', 2)
grid on
ylabel('Scostamento [%]')
xlabel('Tempo [ore]')
title('Scostamento r Percentuale')

% Subplot 3: Scarto Assoluto velocità
subplot(2,2,3)
plot(ct/3600, delta_vel/1000, 'LineWidth', 1.5, 'Color', "r")
hold on
pv = polyfit(ct, delta_vel/1000, 1); 
v_reg = polyval(pv, ct);
plot(ct/3600, v_reg, 'g--', 'LineWidth', 2)
grid on
ylabel('Scostamento Velocità [km/s]')
xlabel('Tempo [ore]')
title('Differenza Assoluta')

% Subplot 4: Scarto Percentuale velocità
subplot(2,2,4)
plot(ct/3600, scartov_percentuale, 'LineWidth', 1.5, 'Color', "b")
hold on
percv = polyfit(ct, scartov_percentuale, 1); 
v_perc_reg = polyval(percv, ct);
plot(ct/3600, v_perc_reg, 'g--', 'LineWidth', 2)
grid on
ylabel('Scostamento [%]')
xlabel('Tempo [ore]')
title('Scostamento v Percentuale')

% Plot parametri orbitali
figure('Name','Andamento parametri orbitali','NumberTitle', 'off')
% i
subplot(3,1,1)
plot(ct/3600,i,'LineWidth', 1.5, 'Color', "b") %toolbox
hold on
plot(ct/3600, ones(length(i),1)*sat_param.i,'LineWidth', 1.5, 'Color',"g") %Custom
hold on
plot(ct/3600, i_gmat_interp,'LineWidth', 1.5, 'Color',[1, 0.5, 0])  %GMAT
grid on
ylabel('i [deg]')
xlabel('Tempo [ore]')
legend("Aerospace toolbox", "Propagatore custom","GMAT")
title('andamento inclinazione')

% raan
subplot(3,1,2)
plot(ct/3600,raan,'LineWidth', 1.5, 'Color', "b") %toolbox
hold on
plot(ct/3600,ones(length(i),1)*sat_param.raan,'LineWidth', 1.5, 'Color', "g") %Custom
hold on
plot(ct/3600, raan_gmat_interp,'LineWidth', 1.5, 'Color',[1, 0.5, 0])  %GMAT
grid on
ylabel('raan [deg]')
xlabel('Tempo [ore]')
legend("Aerospace toolbox", "Propagatore custom","GMAT")
title('andamento raan')

% Argomento di perigeo
subplot(3,1,3)
plot(ct/3600,omega,'LineWidth', 1.5, 'Color', "b") %toolbox
hold on
plot(ct/3600,ones(length(i),1)*sat_param.omega,'LineWidth', 1.5, 'Color', "g") %Custom
hold on
plot(ct/3600, aop_gmat_interp,'LineWidth', 1.5, 'Color',[1, 0.5, 0])  %GMAT
grid on
ylabel('AoP [deg]')
xlabel('Tempo [ore]')
legend("Aerospace toolbox", "Propagatore custom","GMAT")
title('andamento AoP')

