function [] = confronta_risultati(res)
% Estrazione dati
custom_data = res.sat_orbit;
aero_toolbox = res.aerotb_res;

start_time = custom_data.t_date(1);

cr = custom_data.r_eci; % [N x 3]
ct = custom_data.t;     % [N x 1]

tbr = aero_toolbox.pos.Data(3:end);
tbt = seconds(datetime(aero_toolbox.time.Data(3:end,:)) - start_time);

tbr_interp = interp1(tbt, tbr, ct, 'spline');

delta_pos = sqrt(sum((cr - tbr_interp).^2, 2));

r_toolbox = sqrt(sum(tbr_interp.^2, 2));
errore_percentuale = (delta_pos ./ r_toolbox) * 100;

% --- PLOT ---
figure('Name', 'Confronto Propagatori', 'NumberTitle', 'off')

% Subplot 1: Errore Assoluto
subplot(2,1,1)
plot(ct/3600, delta_pos, 'LineWidth', 1.5, 'Color', "r")
grid on
ylabel('Scostamento Posizione [m o km]')
title('Differenza Assoluta (Norma della distanza)')

% Subplot 2: Errore Percentuale
subplot(2,1,2)
plot(ct/3600, errore_percentuale, 'LineWidth', 1.5, 'Color', "b")
grid on
ylabel('Scostamento [%]')
xlabel('Tempo [ore]')
title('Errore Percentuale Relativo alla Distanza dal Centro')
end
