function tb_table = get_thirdBody_table(startTime, duration, stepSize, thirdBody)
% Genera una tabella di posizioni del terzo corpo rispetto alla Terra
% usando le function del toolbox

% INPUT:
% startTime: datetime di inizio
% duration: durata totale in secondi
% stepSize: ogni quanti secondi calcolare la posizione (es. 300s o 600s)
% bodyName: 'Sun' o 'Moon'

t_vec = (0:stepSize:duration)';
dateTime_vec = startTime + seconds(t_vec);
jd_vec = juliandate(dateTime_vec);

pos_km = planetEphemeris(jd_vec, 'Earth', thirdBody);

tb_table.time = t_vec;
tb_table.pos = pos_km' * 1e3; % Convertiamo in metri e trasponiamo [3 x N]

if strcmpi(thirdBody, 'Moon')
    tb_table.mu = 4.9048695e12;
elseif strcmpi(thirdBody, 'Sun')
    tb_table.mu = 1.32712440e20;
else
    error(" %s Non supportato", thirdBody)
end

end