function a_tb = third_body(r,t, tb_table)
% Function per calcolare l'accelereazione dovuta alla perturbazione di
% effetti di terzo corpo. Fa interpolazione dei dati della tabella ottenuta
% con la function: get_thirdBody_table, richiamata prima dell'integrazione
% per ottimizzazione del codice

r_tb = interp1(tb_table.time, tb_table.pos', t, "linear")';
    
d = r_tb - r;
a_tb = tb_table.mu * (d/norm(d)^3 - r_tb/(norm(r_tb)^3));

end


% %% -------------------VERSIONE VECCHIA----------------------------%%
% % (sostituita con get_thirdBody_table prima dell'integrazione)

% function a_tb = third_body(r,t, data, third_body)
% tb = third_body;
% dateTime = data.startTime + seconds(t);
% jd = juliandate(dateTime);
% 
% if tb == "Moon"
%     mu_tb = 4.9048695e12;   % [m³/s²]
% elseif tb == "Sun" 
%     mu_tb  = 1.32712440e20;  % [m³/s²]
% else
%     error("third body non supportato")
% end
% 
% r_tb = planetEphemeris(jd, 'Earth', tb) * 1e3;   % [m]
% r_tb = r_tb';
% 
% d = r_tb - r;
% a_tb = mu_tb * (d/norm(d)^3 - r_tb/(norm(r_tb)^3));
% end
