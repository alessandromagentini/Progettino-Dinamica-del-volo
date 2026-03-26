function [] = plotter(sat_param,sat_orbit,groundtrack3_flag, grundtrack2_flag,plot_eci)

R_terra = 6378e3;                                                          %[m]
start_time = sat_orbit.t_date(1); stop_time = sat_orbit.t_date(end);
dt = sat_orbit.dt;

% PLOT
t_utc = start_time:seconds(dt):stop_time;
lla = eci2lla(sat_orbit.r_eci, datevec(t_utc));
lat = lla(:,1);
lon = lla(:,2);
alt = lla(:,3);
if groundtrack3_flag == 1
    fig = uifigure;
    globe = geoglobe(fig);
    geoplot3(globe,lat, lon, alt, 'LineWidth', 2, 'Color', 'r')            % alt = 0 per avere la ground track "effettiva"
end
if grundtrack2_flag == 1
    figure();
    geoplot(lat, lon, 'r-', 'LineWidth', 1.5);
    geobasemap('satellite');
    title('Ground Track');
end
if plot_eci == 1
    figure('Name', 'Visualizzazione Orbita ECI', 'Color', 'w');
    hold on; grid on; axis equal;
    
    % Rappresentazione della Terra (Sfera)
    [x_s, y_s, z_s] = sphere(50);
    surf(x_s*R_terra, y_s*R_terra, z_s*R_terra, ...
         'FaceColor', [0.3 0.5 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    % Plot della traiettoria
    plot3(sat_orbit.r_eci(:,1), sat_orbit.r_eci(:,2), sat_orbit.r_eci(:,3), ...
          'Color', [0 0.4470 0.7410], 'LineWidth', 2);

    % Plot parametri orbitali
    e_dir = sat_param.e_vec / norm(sat_param.e_vec);
    r_peri_vec = e_dir * (sat_param.rp); 
    r_apo_vec  = -e_dir * (sat_param.ra);
    % Marker per Pericentro (Magenta) e Apocentro (Arancione)
    plot3(r_peri_vec(1), r_peri_vec(2), r_peri_vec(3), 'p', 'MarkerSize', 10, ...
        'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm');
    text(r_peri_vec(1), r_peri_vec(2), r_peri_vec(3), '  R_p', ...
    'Color', 'm', 'FontWeight', 'bold', 'FontSize', 10);
    plot3(r_apo_vec(1), r_apo_vec(2), r_apo_vec(3), 'p', 'MarkerSize', 10, ...
        'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.5 0]);
    text(r_apo_vec(1), r_apo_vec(2), r_apo_vec(3), '  R_a', ...
    'Color', [1 0.5 0], 'FontWeight', 'bold', 'FontSize', 10);

    % Plot samples
    plot3(sat_orbit.samples.position(:,1), sat_orbit.samples.position(:,2), sat_orbit.samples.position(:,3), ...
        'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', 'Samples (3h)');
    r_samp = sat_orbit.samples.position;
    t_samp_date = sat_orbit.samples.date; 
    for k = 1:size(r_samp, 1) % marker temporale dei samples
        label_text = sprintf('  #%d (%s)', k, datestr(t_samp_date(k), 'HH:MM'));
        text(r_samp(k,1) + 200, r_samp(k,2) + 200, r_samp(k,3), label_text, ...
            'FontSize', 8, ...
            'Color', 'k', ...
            'FontWeight', 'bold');
    end

    % --- AGGIUNTA VERSORI ASSI ECI ---
    L = R_terra * 1.5; % Lunghezza delle frecce
    % Asse X (Rosso) - Punto d'Ariete
    quiver3(0,0,0, L, 0, 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5); 
    text(L, 0, 0, '  X (Vernal)', 'Color', 'r', 'FontWeight', 'bold');
    
    % Asse Y (Verde)
    quiver3(0,0,0, 0, L, 0, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    text(0, L, 0, '  Y', 'Color', 'g', 'FontWeight', 'bold');
    
    % Asse Z (Blu) - Polo Nord Celeste
    quiver3(0,0,0, 0, 0, L, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    text(0, 0, L, '  Z (Polo)', 'Color', 'b', 'FontWeight', 'bold');
    % ---------------------------------

    xlabel('X_{ECI} [m]'); ylabel('Y_{ECI} [m]'); zlabel('Z_{ECI} [m]');
    title(['Orbita ECI dal ', datestr(start_time), ' al ', datestr(stop_time)]);
    view(135, 30); % Vista angolata 
    hold off
end