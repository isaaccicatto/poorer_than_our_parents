% Parameterization
target = 'RoW';            % RoW = sum, MM = sum, A = P75, young, B = P25, old, C = P50
suffix = target;             % rename per region if you like (e.g. 'A','B')
file   = 'WPP_calibration.xlsx';
iso    = ["NOR" "ISL" "DNK" "GBR" "USA" "CAN" "GRC" "DEU" "BEL" "IRL" ...
          "FRA" "PRT" "ESP" "LUX" "NLD" "AUT" "CHE" "SWE" "FIN" "JPN" ...
          "ITA" "NZL" "AUS" "KOR" "POL" "HUN" "CZE" "EST" "ISR" "SVK" ...
          "SVN" "CHL" "CRI" "LTU" "LVA"];        % 35 countries (OECD ∩ high-income; no TUR, MEX, COL)
isoCol = "ISO3 Alpha-code";

% Aggregate cohorts
young  = "20-24";
middle = ["25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64"];
old    = ["65-69","70-74","75-79","80-84","85-89","90-94","95-99","100+"];

% Read sheets, filter region & year (Estimates < 2024, Zero-migration >= 2024)
readWPP = @(sheet) local_read(file, sheet);
Est = readWPP('Estimates');
Prj = readWPP('Zero-migration');
Est = Est(ismember(string(Est.(isoCol)), iso) & Est.Year< 2024, :);
Prj = Prj(ismember(string(Prj.(isoCol)), iso) & Prj.Year>=2024, :);
Wall = sortrows([Est; Prj], 'Year');

% Canonical year grid
Year = unique(Wall.Year);
nT = numel(Year);  nC = numel(iso);

% Default level series (only meaningful for RoW/M; NaN otherwise)
N_y = nan(nT,1);  N_m = nan(nT,1);  N_o = nan(nT,1);

switch target
    case {'MM', 'RoW'}
        % ---- population-weighted aggregate: sum levels, then ratios/rate ----
        Wg = groupsummary(Wall, 'Year', 'sum', [young middle old]);
        Wg.Properties.VariableNames = erase(Wg.Properties.VariableNames, 'sum_');
        Age = Wg{:, [young middle old]};
        N_y = Age(:,1);
        N_m = sum(Age(:, 2:1+numel(middle)), 2);
        N_o = sum(Age(:, 2+numel(middle):end), 2);
        psi_y = N_y ./ N_m;
        psi_o = N_o ./ N_m;
        n_m   = [NaN; N_m(2:end) ./ N_m(1:end-1) - 1];

    otherwise
        % Build [year x country] matrices with each country's TIME SERIES intact,
        % then take the percentile ACROSS countries, year by year, per variable.
        PSI_Y = nan(nT,nC);  PSI_O = nan(nT,nC);  N_M = nan(nT,nC);
        for c = 1:nC
            m  = string(Wall.(isoCol)) == iso(c);
            Wc = sortrows(Wall(m,:), 'Year');
            assert(numel(Wc.Year)==nT, 'Region %s does not cover all years.', iso(c));
            [tf,pos] = ismember(Wc.Year, Year);  assert(all(tf));
            Agec = Wc{:, [young middle old]};
            Nyc  = Agec(:,1);
            Nmc  = sum(Agec(:, 2:1+numel(middle)), 2);
            Noc  = sum(Agec(:, 2+numel(middle):end), 2);
            PSI_Y(pos,c) = Nyc ./ Nmc;
            PSI_O(pos,c) = Noc ./ Nmc;
            N_M(pos,c)   = [NaN; Nmc(2:end)./Nmc(1:end-1) - 1];  % growth WITHIN country
        end

        % pick the tail matching the label (independent percentile per variable)
        switch target
            case 'A'   % fast growth, few old, many young
                pn = 75;  po = 25;  py = 75;
            case 'B'     % slow growth, many old, few young
                pn = 25;  po = 75;  py = 25;
            case 'C'  % typical country
                pn = 50;  po = 50;  py = 50;
            otherwise
                error('invalid target: use RoW | MM | A | B | C');
        end
        n_m   = prctile(N_M,   pn, 2);   % NaN in row 1 (all-NaN) -> stays NaN
        psi_o = prctile(PSI_O, po, 2);
        psi_y = prctile(PSI_Y, py, 2);
end

% Smoothing through cubic spline
p_smooth = 0.01;
n_m_s   = csaps(Year, n_m,   p_smooth, Year);
psi_o_s = csaps(Year, psi_o, p_smooth, Year);
psi_y_s = csaps(Year, psi_y, p_smooth, Year);

% Diagnostic plot: raw vs smoothed
figure('Color','w');
subplot(3,1,1);
plot(Year, n_m, '.', 'MarkerSize', 8); hold on;
plot(Year, n_m_s, '-', 'LineWidth', 1.5);
title(sprintf('n\\_m (%s): raw (dots) vs smoothed (line)', target)); grid on;
subplot(3,1,2);
plot(Year, psi_o, '.', 'MarkerSize', 8); hold on;
plot(Year, psi_o_s, '-', 'LineWidth', 1.5);
title('\psi_o: raw vs smoothed'); grid on;
subplot(3,1,3);
plot(Year, psi_y, '.', 'MarkerSize', 8); hold on;
plot(Year, psi_y_s, '-', 'LineWidth', 1.5);
title('\psi_y: raw vs smoothed'); grid on;
xlabel('Year');

% Invert gamma and phi with smoothed trends
omega = 1 - 1/40;  nu = 1 - 1/5;
gamma_path = ((1 + n_m_s) .* psi_o_s - (1-omega)) ./ ...
             [psi_o_s(1); psi_o_s(1:end-1)];
phi_path   =  (1 + n_m_s) .* psi_y_s - nu .* ...
             [psi_y_s(1); psi_y_s(1:end-1)];

% ---- validity check: percentile-mixing can break the demographic identity ----
gok = gamma_path(~isnan(gamma_path));
pok = phi_path(~isnan(phi_path));
assert(all(gok > 0 & gok < 1), ...
    ['gamma outside of (0,1) in %s: inconsistent percentiles. ' ...
     'Smooth more (lower p_smooth) or go back to 25/75.'], target);
assert(all(pok > 0), ...
    'negative phi in %s: inconsistent percentiles.', target);

% Save (using the suffix)
S = struct();
S.Year = Year;
S.(sprintf('gamma_%s_path', suffix)) = gamma_path;
S.(sprintf('phi_%s_path',   suffix)) = phi_path;
S.(sprintf('n_m_%s',   suffix))      = n_m;
S.(sprintf('psi_y_%s', suffix))      = psi_y;
S.(sprintf('psi_o_%s', suffix))      = psi_o;
S.(sprintf('n_m_%s_s',   suffix))    = n_m_s;
S.(sprintf('psi_y_%s_s', suffix))    = psi_y_s;
S.(sprintf('psi_o_%s_s', suffix))    = psi_o_s;
matfile = fullfile(sprintf('demo_paths_%s.mat', suffix));
save(matfile, '-struct', 'S');

% CSV  (N_y/N_m/N_o are NaN for percentile targets — no coherent level exists)
T = table(Year, N_y, N_m, N_o, psi_y, psi_o, n_m, ...
          psi_y_s, psi_o_s, n_m_s, gamma_path, phi_path, ...
'VariableNames', {'Year', ...
        sprintf('N_y_%s',suffix), sprintf('N_m_%s',suffix), sprintf('N_o_%s',suffix), ...
        sprintf('psi_y_%s',suffix), sprintf('psi_o_%s',suffix), sprintf('n_m_%s',suffix), ...
        sprintf('psi_y_%s_s',suffix), sprintf('psi_o_%s_s',suffix), sprintf('n_m_%s_s',suffix), ...
        sprintf('gamma_%s_path',suffix), sprintf('phi_%s_path',suffix)});
csvfile = sprintf('demo_paths_%s.csv', suffix);
writetable(T, csvfile);

% print relevant variables (1970 -> 2070)
y0 = 1970;  y1 = 2070;
i0 = find(Year == y0);
i1 = find(Year == y1);
assert(~isempty(i0) && ~isempty(i1), 'Period outside the interval of Year.');
fprintf('\n=== Variation %d -> %d (region %s) ===\n', y0, y1, suffix);
fprintf('  d gamma = %+.10f   (%.10f -> %.10f)\n', ...
        gamma_path(i1)-gamma_path(i0), gamma_path(i0), gamma_path(i1));
fprintf('  d phi   = %+.10f   (%.10f -> %.10f)\n', ...
        phi_path(i1)-phi_path(i0),     phi_path(i0),   phi_path(i1));
fprintf('  d psi_o = %+.10f   (%.10f -> %.10f)\n', ...
        psi_o(i1)-psi_o(i0),           psi_o(i0),      psi_o(i1));
fprintf('  d psi_y = %+.10f   (%.10f -> %.10f)\n', ...
        psi_y(i1)-psi_y(i0),           psi_y(i0),      psi_y(i1));
fprintf('  d n_m   = %+.10f pp (%.10f -> %.10f)\n', ...
        (n_m(i1)-n_m(i0))*100,         n_m(i0),        n_m(i1));

% subfunction
function T = local_read(file, sheet)
    opts = detectImportOptions(file, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
    opts.VariableNamesRange = 'A17';
    opts.DataRange          = 'A18';
    T = readtable(file, opts);
end