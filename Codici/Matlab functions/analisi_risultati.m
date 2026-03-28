function [] = analisi_risultati(sat_param,res)
% Function per l'analisi e plot dei risultati

mu_terra = 398600.4418;   

% Estrazione dati
custom_data = res.sat_orbit;
aero_toolbox = res.aerotb_res;

start_time = custom_data.t_date(1);

cr = custom_data.r_eci; 
cv = custom_data.v_eci;
ct = custom_data.t;     

tbr = aero_toolbox.pos_eci(3:end,:);
tbv = aero_toolbox.vel_eci(3:end,:);
tbt = seconds(datetime(aero_toolbox.time.Data(3:end,:)) - start_time);

tbr_interp = interp1(tbt, tbr, ct, 'spline');
tbv_interp = interp1(tbt, tbv, ct, 'spline');

delta_pos = sqrt(sum((cr + tbr_interp).^2, 2));
delta_vel = sqrt(sum((cv + tbv_interp).^2, 2));

scartor_percentuale = (delta_pos ./ norm(tbr_interp)) * 100;
scartov_percentuale = (delta_vel ./ norm(tbv_interp)) * 100;

% Andamento dei parametri orbitali
i = zeros(length(tbr_interp),1); raan  = zeros(length(tbr_interp),1); T = zeros(length(tbr_interp),1);
for idx = 1:length(tbr_interp)
    param = get_parametri_orbitali(tbr_interp(idx,:)/1000,tbv_interp(idx,:)/1000, mu_terra);
    i(idx) = param.i;
    raan(idx) = param.raan;
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
subplot(2,1,1)
plot(ct/3600,i,'LineWidth', 1.5, 'Color', "r")
hold on
plot(ct/3600, ones(length(i),1)*sat_param.i,'LineWidth', 1.5, 'Color',"m")
grid on
ylabel('i [deg]')
xlabel('Tempo [ore]')
legend("Aerospace toolbox", "Propagatore custom")
title('andamento inclinazione')

% raan
subplot(2,1,2)
plot(ct/3600,raan,'LineWidth', 1.5, 'Color', "b")
hold on
plot(ct/3600,ones(length(i),1)*sat_param.raan,'LineWidth', 1.5, 'Color', "g")
grid on
ylabel('raan [deg]')
xlabel('Tempo [ore]')
legend("Aerospace toolbox", "Propagatore custom")
title('andamento raan')

