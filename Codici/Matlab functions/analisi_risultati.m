function [] = analisi_risultati(res)
% Function per l'analisi e plot dei risultati

mu_terra = 398600.4418;   

% Estrazione dati
custom_data = res.sat_orbit;
aero_toolbox = res.aerotb_res;

start_time = custom_data.t_date(1);

cr = custom_data.r_eci; 
cv = custom_data.v_eci;
ct = custom_data.t;     

tbr = aero_toolbox.pos.Data(3:end,:);
tbv = aero_toolbox.vel.Data(3:end,:);
tbt = seconds(datetime(aero_toolbox.time.Data(3:end,:)) - start_time);

tbr_interp = interp1(tbt, tbr, ct, 'spline');
tbv_interp = interp1(tbt, tbv, ct, 'spline');

delta_pos = sqrt(sum((cr + tbr_interp).^2, 2));
% delta_vel

r_toolbox = norm(tbr_interp);
errore_percentuale = (delta_pos ./ r_toolbox) * 100;

% Andamento dei parametri orbitali
i = zeros(length(tbr_interp),1); raan  = zeros(length(tbr_interp),1); 
for idx = 1:length(tbr_interp)
    param = get_parametri_orbitali(tbr_interp(idx,:)/1000,tbv_interp(idx,:)/1000, mu_terra);
    i(idx) = param.i;
    raan(idx) = param.raan;
end

%% --- PLOT ---

% Plot scostamenti tra i modelli
figure('Name', 'Confronto Propagatori', 'NumberTitle', 'off')
% Subplot 1: Errore Assoluto
subplot(2,1,1)
plot(ct/3600, delta_pos/1000, 'LineWidth', 1.5, 'Color', "r")
grid on
ylabel('Scostamento Posizione [km]')
xlabel('Tempo [ore]')
title('Differenza Assoluta (Norma della distanza)')

% Subplot 2: Errore Percentuale
subplot(2,1,2)
plot(ct/3600, errore_percentuale, 'LineWidth', 1.5, 'Color', "b")
grid on
ylabel('Scostamento [%]')
xlabel('Tempo [ore]')
title('Scostamento Percentuale')

% Plot parametri orbitali
figure('Name','Andamento parametri orbitali','NumberTitle', 'off')
% i
subplot(2,1,1)
plot(ct/3600,i,'LineWidth', 1.5, 'Color', "r")
grid on
ylabel('i [deg]')
xlabel('Tempo [ore]')
title('andamento inclinazione')

% raan
subplot(2,1,2)
plot(ct/3600,raan,'LineWidth', 1.5, 'Color', "b")
grid on
ylabel('raan [deg]')
xlabel('Tempo [ore]')
title('andamento raan')

