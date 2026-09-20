% addpath '/Applications/Dynare/7.0-arm64/matlab' - Paste this in the terminal
% ============================================================================
% DEBT REFORM: unanticipated gradual rise in the debt rule b/y, 0.60 -> 0.90,
% raised-cosine over reform_span periods starting at t_reform (2025).
% b_y is EXOGENOUS (varexo); fy_target endogenous via three steady states
% (ini / dem with b_y=0.6 / ref with b_y=0.9), reform regime from t_reform on.
% Same architecture as the pension reforms: constant beliefs, scalar exo per
% iteration, exo paths recorded in the loop (export cannot drift).
% ============================================================================

%%%%%%% VARIABLES %%%%%%%
var

% Endogenous
n_MM 
psi_o_MM
psi_y_MM
n_y_MM
n_m_MM
n_o_MM
mu_o_MM
mu_m_MM
mu_y_MM
s_o_MM
s_m_MM
s_y_MM
h_m_MM
h_y_MM
a_MM
a_o_MM
a_m_MM
a_y_MM
z_o_MM
z_m_MM
z_y_MM
c_o_MM
c_m_MM
c_y_MM
y_MM
k_MM
w_MM
inv_MM
b_MM
gov_MM
e_MM
t_MM
f_MM
nx_MM
tilde_R_MM
g_MM
q_MM

% Auxiliary (endogenous)
tilde_s_m_MM
tilde_s_y_MM
tilde_s_o_MM
Omega_m_MM
Omega_y_MM
eta_o_MM
eta_m_MM
eta_y_MM
walras_MM
;

varexo gamma_MM, phi_MM, R_star, fy_target_MM, b_y;


%%%%%%% PARAMETERS %%%%%%%

parameters alpha, delta, beta, sigma, g_y, repl, x, omega, nu, Lambda_MM;

omega = 1 - 1 / 40;     % worklife = 40
nu = 1 - 1 / 5;     % university = 5
alpha = 0.67;
delta = 0.1;
beta = 1.025;
sigma = 0.5;
g_y = 0.25;
repl = 0.4;
x = 0.005;
Lambda_MM = 0.1;


%%%%%%% EQUATIONS %%%%%%%

model;

% 0. Normalization
g_MM = (1 + n_MM) * (1 + x);
q_MM = tilde_R_MM / g_MM;

% 1. Demographics
n_y_MM = (nu + phi_MM * 1 / psi_y_MM(-1)) - 1;
n_m_MM = (omega + (1 - nu) * psi_y_MM(-1)) - 1;
n_o_MM = ((1 - omega) * 1 / psi_o_MM(-1) + gamma_MM) - 1;
n_MM = phi_MM * tilde_s_m_MM(-1) - (1 - gamma_MM) * tilde_s_o_MM(-1);
psi_y_MM = (nu * psi_y_MM(-1) + phi_MM) * 1 / (1 + n_m_MM);
psi_o_MM = (1 - omega + gamma_MM * psi_o_MM(-1)) * 1 / (1 + n_m_MM);
tilde_s_y_MM = (nu * tilde_s_y_MM(-1) + phi_MM * tilde_s_m_MM(-1)) * 1 / (1 + n_MM);
tilde_s_m_MM = (omega * tilde_s_m_MM(-1) + (1 - nu) * tilde_s_y_MM(-1)) * 1 / (1 + n_MM);
tilde_s_o_MM = ((1 - omega) * tilde_s_m_MM(-1) + gamma_MM * tilde_s_o_MM(-1)) * 1 / (1 + n_MM); 

% 2. Household
%% MPCs
1 = (1 + gamma_MM(+1) * beta^sigma * tilde_R_MM(+1)^(sigma - 1) * (1 / mu_o_MM(+1))) * mu_o_MM; 
1 = (1 + beta^sigma * (Omega_m_MM(+1) * tilde_R_MM(+1))^(sigma - 1) * (1 / mu_m_MM(+1))) * mu_m_MM;
1 = (1 + beta^sigma * (Omega_y_MM(+1) * tilde_R_MM(+1))^(sigma - 1) * (1 / mu_y_MM(+1))) * mu_y_MM;
Omega_m_MM = omega + (1 - omega) * (mu_o_MM / mu_m_MM)^(1 / (1 - sigma));
Omega_y_MM = nu + (1 - nu) * (mu_m_MM / mu_y_MM)^(1 / (1 - sigma));

%% PDVs
s_o_MM = e_MM + gamma_MM(+1) * g_MM(+1) / (tilde_R_MM(+1) * (1 + n_o_MM(+1))) * s_o_MM(+1); 
s_m_MM = g_MM(+1) / (Omega_m_MM(+1) * tilde_R_MM(+1)) * ((Omega_m_MM(+1) - omega(+1)) / ((1 + n_o_MM(+1)) * psi_o_MM) * s_o_MM(+1) + omega(+1) / (1 + n_m_MM(+1)) * s_m_MM(+1)); 
s_y_MM = g_MM(+1) / (Omega_y_MM(+1) * tilde_R_MM(+1)) * ((Omega_y_MM(+1) - nu(+1)) * psi_y_MM / (1 + n_m_MM(+1)) * s_m_MM(+1) + nu(+1) / (1 + n_y_MM(+1)) * s_y_MM(+1)); 
h_m_MM = w_MM - t_MM + g_MM(+1) * omega(+1) / (Omega_m_MM(+1) * tilde_R_MM(+1) * (1 + n_m_MM(+1))) * h_m_MM(+1); 
h_y_MM = g_MM(+1) / (Omega_y_MM(+1) * tilde_R_MM(+1)) * ((Omega_y_MM(+1) - nu(+1)) * psi_y_MM / (1 + n_m_MM(+1)) * h_m_MM(+1) + nu(+1) / (1 + n_y_MM(+1)) * h_y_MM(+1));

% Assets
a_o_MM = q_MM * (1 - mu_o_MM) * (eta_o_MM(-1) + (1 - omega) * eta_m_MM(-1)) * a_MM(-1) + e_MM - mu_o_MM * s_o_MM;
a_m_MM = q_MM * (1 - mu_m_MM) * (omega * eta_m_MM(-1) + (1 - nu) * eta_y_MM(-1)) * a_MM(-1) + w_MM - t_MM - mu_m_MM * (h_m_MM + s_m_MM);
a_y_MM = q_MM * (1 - mu_y_MM) * nu * eta_y_MM(-1) * a_MM(-1) - mu_y_MM * (h_y_MM + s_y_MM);
eta_o_MM = a_o_MM / a_MM;
eta_m_MM = a_m_MM / a_MM;
eta_y_MM = a_y_MM / a_MM;
eta_y_MM + eta_m_MM + eta_o_MM = 1;

% Wealth
z_o_MM = q_MM * (eta_o_MM(-1) + (1 - omega) * eta_m_MM(-1)) * a_MM(-1) + s_o_MM;
z_m_MM = q_MM * (omega * eta_m_MM(-1) + (1 - nu) * eta_y_MM(-1)) * a_MM(-1) + h_m_MM + s_m_MM;
z_y_MM = q_MM * nu * eta_y_MM(-1) * a_MM(-1) + h_y_MM + s_y_MM;

% Consumption
c_o_MM = mu_o_MM * z_o_MM;
c_m_MM = mu_m_MM * z_m_MM;
c_y_MM = mu_y_MM * z_y_MM;


% 3. Firm
tilde_R_MM(+1) = (1 - alpha) * y_MM(+1) / k_MM * g_MM(+1) + (1 - delta);
y_MM = tilde_s_m_MM^alpha * (k_MM(-1) / g_MM)^(1 - alpha);
w_MM = alpha * y_MM;
inv_MM = k_MM - (1 - delta) / g_MM * k_MM(-1);


% 4. Government
b_MM = b_y * y_MM;
gov_MM = g_y * y_MM;
e_MM = repl * w_MM / tilde_s_m_MM * tilde_s_o_MM;
t_MM = e_MM + gov_MM + q_MM * b_MM(-1) - b_MM;


% 5. Balance of Payments
a_MM = k_MM + b_MM + f_MM;
nx_MM = f_MM - q_MM * f_MM(-1);


% 6. Imperfect Capital Mobility
tilde_R_MM = R_star + Lambda_MM * (exp(-(f_MM(-1)/y_MM(-1) - fy_target_MM)) - 1);


% Walras check
walras_MM = y_MM - c_o_MM - c_m_MM - c_y_MM - inv_MM - gov_MM - nx_MM;

end;


%%%%%%% CALIBRATION %%%%%%%

D = load('demo_paths_MM.mat');   
R_data = load('R_star_path.mat'); 
R_star_path = R_data.R_star;    

y_start = 1970;
y_end   = 2070;

% horizon defined HERE (before the reform block uses it)
steps        = 100;
horizon_plot = 200;

i0 = find(D.Year == y_start);
i1 = find(D.Year == y_end);
assert(~isempty(i0) && ~isempty(i1), 'Period outside the interval of Year');

Year_win        = D.Year(i0:i1);            
gamma_MM_data    = D.gamma_MM_path(i0:i1);    
phi_MM_data      = D.phi_MM_path(i0:i1);

gamma_MM_start = gamma_MM_data(1);            
gamma_MM_end   = gamma_MM_data(end);          
phi_MM_start   = phi_MM_data(1);
phi_MM_end     = phi_MM_data(end);
R_star_start = R_star_path(1);      
R_star_end   = R_star_path(end);

% ---- REFORM: public debt expansion (unanticipated gradual) -----------------
by_0        = 0.60;      % pre-reform debt rule
by_new      = 1.05;      % post-reform target
t_reform    = 55;        % onset (1970+55 = 2025)
reform_span = 45;        % phase-in of the rule AND target migration (cosine)

by_path = by_0 * ones(horizon_plot + 1, 1);   % index t+1 = model period t
for tt = t_reform:horizon_plot
    prog = min((tt - t_reform + 1)/reform_span, 1);
    by_path(tt+1) = by_0 + (by_new - by_0) * (1 - cos(pi * prog)) / 2;   % raised-cosine
end
by_start = by_path(1);
by_end   = by_path(end);

initval;                     
  gamma_MM = gamma_MM_start;   
  phi_MM   = phi_MM_start;
  R_star  = R_star_start; 
  b_y     = by_0;
end;

% Initial (exogeneous) BGP
oo_.exo_steady_state(strcmp(M_.exo_names,'gamma_MM')) = gamma_MM_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'phi_MM'))   = phi_MM_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'R_star'))  = R_star_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'b_y'))     = by_0;

% exo order assumed: [gamma; phi; R_star; fy_target; b_y]
% (confirm with disp(M_.exo_names) after first compile)
ys_ini = ext_path_MM_by_steadystate(zeros(M_.endo_nbr,1), ...
    [gamma_MM_start; phi_MM_start; R_star_start; 0; by_0], M_, options_);
ys_dem = ext_path_MM_by_steadystate(zeros(M_.endo_nbr,1), ...
    [gamma_MM_end;   phi_MM_end;   R_star_end;   0; by_0], M_, options_);
ys_ref = ext_path_MM_by_steadystate(zeros(M_.endo_nbr,1), ...
    [gamma_MM_end;   phi_MM_end;   R_star_end;   0; by_new], M_, options_);

iy  = strcmp(M_.endo_names,'y_MM');
iff = strcmp(M_.endo_names,'f_MM');

fy_target_MM_start = ys_ini(iff)/ys_ini(iy);
fy_target_MM_dem   = ys_dem(iff)/ys_dem(iy);   % demographic target (= baseline)
fy_target_MM_ref   = ys_ref(iff)/ys_ref(iy);   % post-reform target (b_y = 0.9)

fprintf('fy_target: start=%.4f  dem=%.4f  ref=%.4f  (salto reforma = %+.4f)\n', ...
        fy_target_MM_start, fy_target_MM_dem, fy_target_MM_ref, ...
        fy_target_MM_ref - fy_target_MM_dem);

oo_.exo_steady_state(strcmp(M_.exo_names,'fy_target_MM')) = fy_target_MM_start;

steady;
check;


%%%%%%% EXTENDED-PATH %%%%%%%

options_.noprint = 1;
perfect_foresight_setup(periods=400);

n_vars = M_.endo_nbr;

realized_path = zeros(n_vars, horizon_plot + 1);
id_gamma_MM     = find(strcmp(M_.exo_names, 'gamma_MM'));
id_phi_MM       = find(strcmp(M_.exo_names, 'phi_MM'));
id_R_star      = find(strcmp(M_.exo_names, 'R_star'));
id_fy_target_MM = find(strcmp(M_.exo_names, 'fy_target_MM'));    
id_by          = find(strcmp(M_.exo_names, 'b_y'));

% recorded exo paths (filled IN the loop; export uses these directly)
gamma_real     = zeros(horizon_plot+1,1);  gamma_real(1)     = gamma_MM_start;
phi_real       = zeros(horizon_plot+1,1);  phi_real(1)       = phi_MM_start;
R_star_real    = zeros(horizon_plot+1,1);  R_star_real(1)    = R_star_start;
fy_target_real = zeros(horizon_plot+1,1);  fy_target_real(1) = fy_target_MM_start;
by_real        = zeros(horizon_plot+1,1);  by_real(1)        = by_0;

current_state = oo_.steady_state;
realized_path(:, 1) = current_state;

disp('>>> STARTING ROLLING MIT SHOCKS LOOP... <<<');

for i = 1:horizon_plot

    % ---- demographics (unchanged) ------------------------------------------
    if i <= steps
        current_gamma_MM = gamma_MM_data(i + 1);  % i=1 -> gamma_MM_data(2) = 1971
        current_phi_MM   = phi_MM_data(i + 1);
        w_dem = (1 - cos(pi * i / steps)) / 2;
    else
        current_gamma_MM = gamma_MM_end;
        current_phi_MM   = phi_MM_end;
        w_dem = 1;
    end

    % ---- external target: demographic path + reform migration -------------
    fy_dem = fy_target_MM_start + (fy_target_MM_dem - fy_target_MM_start) * w_dem;
    prog   = min(max((i - t_reform + 1)/reform_span, 0), 1);
    w_ref  = (1 - cos(pi * prog)) / 2;
    current_fy_target_MM = fy_dem + (fy_target_MM_ref - fy_target_MM_dem) * w_ref;

    current_R_star = R_star_path(i + 1);
    current_by     = by_path(i + 1);

    % ---- record what is fed (export reads these; no drift possible) --------
    gamma_real(i+1)     = current_gamma_MM;
    phi_real(i+1)       = current_phi_MM;
    R_star_real(i+1)    = current_R_star;
    fy_target_real(i+1) = current_fy_target_MM;
    by_real(i+1)        = current_by;

    oo_.exo_simul(:, id_gamma_MM)     = current_gamma_MM;
    oo_.exo_simul(:, id_phi_MM)       = current_phi_MM;
    oo_.exo_simul(:, id_R_star)      = current_R_star;
    oo_.exo_simul(:, id_fy_target_MM) = current_fy_target_MM;
    oo_.exo_simul(:, id_by)          = current_by;
    
    oo_.endo_simul = repmat(current_state, 1, options_.periods + 2);
    
    perfect_foresight_solver;
    
    current_state = oo_.endo_simul(:, 2);
    realized_path(:, i+1) = current_state;
    
    fprintf('Year %d of %d (calendar %d) solved under constant beliefs.\n', ...
            i, horizon_plot, y_start + i);
end
disp('>>> THE LOOP IS FINISHED! <<<');


%% Check for PCM in BGP
fprintf('R_tilde_MM_end - R_star_end = %.2e\n', ...
    realized_path(strcmp(M_.endo_names,'tilde_R_MM'), end) - R_star_end);


%%%%%%% AUTOMATIC PLOTS %%%%%%%

get_series = @(name) realized_path(find(strcmp(M_.endo_names, name)), 1:horizon_plot);

vars_to_plot = {'y_MM', 'k_MM', 'f_MM', 'tilde_R_MM', ... 
                'c_y_MM', 'c_m_MM', 'c_o_MM', ...
                'eta_y_MM', 'eta_m_MM', 'eta_o_MM', ...
                'a_y_MM', 'a_m_MM', 'a_o_MM'};

% Figure 1: Macro
figure('Name', 'Macro Aggregates', 'Position', [100 100 1000 600]);
subplot(2,2,1); plot(0:horizon_plot-1, get_series('y_MM'), 'LineWidth', 1.5); title('Output (y\_MM)'); grid on;
subplot(2,2,2); plot(0:horizon_plot-1, get_series('k_MM'), 'LineWidth', 1.5); title('Capital (k\_MM)'); grid on;
subplot(2,2,3); plot(0:horizon_plot-1, get_series('f_MM'), 'LineWidth', 1.5); title('NFA (f\_MM)'); grid on;
subplot(2,2,4); plot(0:horizon_plot-1, get_series('tilde_R_MM'), 'LineWidth', 1.5); title('Interest Rate'); grid on;

% Figure 2: Consumption by cohort
figure('Name', 'Consumption by Cohort', 'Position', [100 100 1000 400]);
plot(0:horizon_plot-1, get_series('c_y_MM'), 'b', 'LineWidth', 1.5); hold on;
plot(0:horizon_plot-1, get_series('c_m_MM'), 'r', 'LineWidth', 1.5);
plot(0:horizon_plot-1, get_series('c_o_MM'), 'k', 'LineWidth', 1.5);
legend('Young', 'Middle', 'Old'); title('Consumption'); grid on;

% Figure 3: Wealth shares
figure('Name', 'Wealth Shares', 'Position', [100 100 1000 400]);
plot(0:horizon_plot-1, get_series('eta_y_MM'), 'b', 'LineWidth', 1.5); hold on;
plot(0:horizon_plot-1, get_series('eta_m_MM'), 'r', 'LineWidth', 1.5);
plot(0:horizon_plot-1, get_series('eta_o_MM'), 'k', 'LineWidth', 1.5);
legend('Young', 'Middle', 'Old'); title('\eta (wealth shares)'); grid on;

for f = 1:3
    figure(f);
    saveas(gcf, sprintf('output_MM/aging_MM_fig_%d.png', f));  
end


%%%%%%% EXPORTING DATA TO CSV %%%%%%%
disp('>>> EXPORTING FULL DUMP... <<<');

N = horizon_plot + 1;   % t = 0 (SS inicial) ate t = horizon_plot

% All endogeneous
T_full = array2table(realized_path(:, 1:N)', 'VariableNames', M_.endo_names');
T_full = addvars(T_full, (0:horizon_plot)', 'Before', 1, 'NewVariableNames', 'Time');

% All exogenous — straight from the recorded paths (no recompute, no drift)
T_full = addvars(T_full, gamma_real, phi_real, R_star_real, fy_target_real, by_real, ...
    'NewVariableNames', {'gamma_MM', 'phi_MM', 'R_star', 'fy_target_MM', 'b_y'});

% All (constant) parameters
T_full = [T_full, array2table(repmat(M_.params', N, 1), 'VariableNames', M_.param_names')];

% Variables in efficiency units
T_full.w_indiv_MM = T_full.w_MM  ./ T_full.tilde_s_m_MM;
T_full.t_indiv_MM = T_full.t_MM  ./ T_full.tilde_s_m_MM;
T_full.e_indiv_MM = T_full.e_MM  ./ T_full.tilde_s_o_MM;
T_full.c_y_indiv_MM = T_full.c_y_MM ./ T_full.tilde_s_y_MM;
T_full.c_m_indiv_MM = T_full.c_m_MM ./ T_full.tilde_s_m_MM;
T_full.c_o_indiv_MM = T_full.c_o_MM ./ T_full.tilde_s_o_MM;
T_full.a_y_indiv_MM = T_full.a_y_MM ./ T_full.tilde_s_y_MM;
T_full.a_m_indiv_MM = T_full.a_m_MM ./ T_full.tilde_s_m_MM;
T_full.a_o_indiv_MM = T_full.a_o_MM ./ T_full.tilde_s_o_MM;
T_full(:, startsWith(T_full.Properties.VariableNames,'AUX')) = [];

writetable(T_full, 'output_MM_by/transition_results_MM_by.csv');
fprintf('CSV: %d rows x %d columns\n', height(T_full), width(T_full));
