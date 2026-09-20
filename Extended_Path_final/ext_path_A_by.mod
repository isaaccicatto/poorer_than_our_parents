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
n_A 
psi_o_A
psi_y_A
n_y_A
n_m_A
n_o_A
mu_o_A
mu_m_A
mu_y_A
s_o_A
s_m_A
s_y_A
h_m_A
h_y_A
a_A
a_o_A
a_m_A
a_y_A
z_o_A
z_m_A
z_y_A
c_o_A
c_m_A
c_y_A
y_A
k_A
w_A
inv_A
b_A
gov_A
e_A
t_A
f_A
nx_A
tilde_R_A
g_A
q_A

% Auxiliary (endogenous)
tilde_s_m_A
tilde_s_y_A
tilde_s_o_A
Omega_m_A
Omega_y_A
eta_o_A
eta_m_A
eta_y_A
walras_A
;

varexo gamma_A, phi_A, R_star, fy_target_A, b_y;


%%%%%%% PARAMETERS %%%%%%%

parameters alpha, delta, beta, sigma, g_y, repl, x, omega, nu, Lambda_A;

omega = 1 - 1 / 40;     % worklife = 40
nu = 1 - 1 / 5;     % university = 5
alpha = 0.67;
delta = 0.1;
beta = 1.025;
sigma = 0.5;
g_y = 0.25;
repl = 0.4;
x = 0.005;
Lambda_A = 0.1;


%%%%%%% EQUATIONS %%%%%%%

model;

% 0. Normalization
g_A = (1 + n_A) * (1 + x);
q_A = tilde_R_A / g_A;

% 1. Demographics
n_y_A = (nu + phi_A * 1 / psi_y_A(-1)) - 1;
n_m_A = (omega + (1 - nu) * psi_y_A(-1)) - 1;
n_o_A = ((1 - omega) * 1 / psi_o_A(-1) + gamma_A) - 1;
n_A = phi_A * tilde_s_m_A(-1) - (1 - gamma_A) * tilde_s_o_A(-1);
psi_y_A = (nu * psi_y_A(-1) + phi_A) * 1 / (1 + n_m_A);
psi_o_A = (1 - omega + gamma_A * psi_o_A(-1)) * 1 / (1 + n_m_A);
tilde_s_y_A = (nu * tilde_s_y_A(-1) + phi_A * tilde_s_m_A(-1)) * 1 / (1 + n_A);
tilde_s_m_A = (omega * tilde_s_m_A(-1) + (1 - nu) * tilde_s_y_A(-1)) * 1 / (1 + n_A);
tilde_s_o_A = ((1 - omega) * tilde_s_m_A(-1) + gamma_A * tilde_s_o_A(-1)) * 1 / (1 + n_A); 

% 2. Household
%% MPCs
1 = (1 + gamma_A(+1) * beta^sigma * tilde_R_A(+1)^(sigma - 1) * (1 / mu_o_A(+1))) * mu_o_A; 
1 = (1 + beta^sigma * (Omega_m_A(+1) * tilde_R_A(+1))^(sigma - 1) * (1 / mu_m_A(+1))) * mu_m_A;
1 = (1 + beta^sigma * (Omega_y_A(+1) * tilde_R_A(+1))^(sigma - 1) * (1 / mu_y_A(+1))) * mu_y_A;
Omega_m_A = omega + (1 - omega) * (mu_o_A / mu_m_A)^(1 / (1 - sigma));
Omega_y_A = nu + (1 - nu) * (mu_m_A / mu_y_A)^(1 / (1 - sigma));

%% PDVs
s_o_A = e_A + gamma_A(+1) * g_A(+1) / (tilde_R_A(+1) * (1 + n_o_A(+1))) * s_o_A(+1); 
s_m_A = g_A(+1) / (Omega_m_A(+1) * tilde_R_A(+1)) * ((Omega_m_A(+1) - omega(+1)) / ((1 + n_o_A(+1)) * psi_o_A) * s_o_A(+1) + omega(+1) / (1 + n_m_A(+1)) * s_m_A(+1)); 
s_y_A = g_A(+1) / (Omega_y_A(+1) * tilde_R_A(+1)) * ((Omega_y_A(+1) - nu(+1)) * psi_y_A / (1 + n_m_A(+1)) * s_m_A(+1) + nu(+1) / (1 + n_y_A(+1)) * s_y_A(+1)); 
h_m_A = w_A - t_A + g_A(+1) * omega(+1) / (Omega_m_A(+1) * tilde_R_A(+1) * (1 + n_m_A(+1))) * h_m_A(+1); 
h_y_A = g_A(+1) / (Omega_y_A(+1) * tilde_R_A(+1)) * ((Omega_y_A(+1) - nu(+1)) * psi_y_A / (1 + n_m_A(+1)) * h_m_A(+1) + nu(+1) / (1 + n_y_A(+1)) * h_y_A(+1));

% Assets
a_o_A = q_A * (1 - mu_o_A) * (eta_o_A(-1) + (1 - omega) * eta_m_A(-1)) * a_A(-1) + e_A - mu_o_A * s_o_A;
a_m_A = q_A * (1 - mu_m_A) * (omega * eta_m_A(-1) + (1 - nu) * eta_y_A(-1)) * a_A(-1) + w_A - t_A - mu_m_A * (h_m_A + s_m_A);
a_y_A = q_A * (1 - mu_y_A) * nu * eta_y_A(-1) * a_A(-1) - mu_y_A * (h_y_A + s_y_A);
eta_o_A = a_o_A / a_A;
eta_m_A = a_m_A / a_A;
eta_y_A = a_y_A / a_A;
eta_y_A + eta_m_A + eta_o_A = 1;

% Wealth
z_o_A = q_A * (eta_o_A(-1) + (1 - omega) * eta_m_A(-1)) * a_A(-1) + s_o_A;
z_m_A = q_A * (omega * eta_m_A(-1) + (1 - nu) * eta_y_A(-1)) * a_A(-1) + h_m_A + s_m_A;
z_y_A = q_A * nu * eta_y_A(-1) * a_A(-1) + h_y_A + s_y_A;

% Consumption
c_o_A = mu_o_A * z_o_A;
c_m_A = mu_m_A * z_m_A;
c_y_A = mu_y_A * z_y_A;


% 3. Firm
tilde_R_A(+1) = (1 - alpha) * y_A(+1) / k_A * g_A(+1) + (1 - delta);
y_A = tilde_s_m_A^alpha * (k_A(-1) / g_A)^(1 - alpha);
w_A = alpha * y_A;
inv_A = k_A - (1 - delta) / g_A * k_A(-1);


% 4. Government
b_A = b_y * y_A;
gov_A = g_y * y_A;
e_A = repl * w_A / tilde_s_m_A * tilde_s_o_A;
t_A = e_A + gov_A + q_A * b_A(-1) - b_A;


% 5. Balance of Payments
a_A = k_A + b_A + f_A;
nx_A = f_A - q_A * f_A(-1);


% 6. Imperfect Capital Mobility
tilde_R_A = R_star + Lambda_A * (exp(-(f_A(-1)/y_A(-1) - fy_target_A)) - 1);


% Walras check
walras_A = y_A - c_o_A - c_m_A - c_y_A - inv_A - gov_A - nx_A;

end;


%%%%%%% CALIBRATION %%%%%%%

D = load('demo_paths_A.mat');   
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
gamma_A_data    = D.gamma_A_path(i0:i1);    
phi_A_data      = D.phi_A_path(i0:i1);

gamma_A_start = gamma_A_data(1);            
gamma_A_end   = gamma_A_data(end);          
phi_A_start   = phi_A_data(1);
phi_A_end     = phi_A_data(end);
R_star_start = R_star_path(1);      
R_star_end   = R_star_path(end);

% ---- REFORM: public debt expansion (unanticipated gradual) -----------------
by_0        = 0.60;      % pre-reform debt rule
by_new      = 0.90;      % post-reform target
t_reform    = 55;        % onset (1970+55 = 2025)
reform_span = 30;        % phase-in of the rule AND target migration (cosine)

by_path = by_0 * ones(horizon_plot + 1, 1);   % index t+1 = model period t
for tt = t_reform:horizon_plot
    prog = min((tt - t_reform + 1)/reform_span, 1);
    by_path(tt+1) = by_0 + (by_new - by_0) * (1 - cos(pi * prog)) / 2;   % raised-cosine
end
by_start = by_path(1);
by_end   = by_path(end);

initval;                     
  gamma_A = gamma_A_start;   
  phi_A   = phi_A_start;
  R_star  = R_star_start; 
  b_y     = by_0;
end;

% Initial (exogeneous) BGP
oo_.exo_steady_state(strcmp(M_.exo_names,'gamma_A')) = gamma_A_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'phi_A'))   = phi_A_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'R_star'))  = R_star_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'b_y'))     = by_0;

% exo order assumed: [gamma; phi; R_star; fy_target; b_y]
% (confirm with disp(M_.exo_names) after first compile)
ys_ini = ext_path_A_by_steadystate(zeros(M_.endo_nbr,1), ...
    [gamma_A_start; phi_A_start; R_star_start; 0; by_0], M_, options_);
ys_dem = ext_path_A_by_steadystate(zeros(M_.endo_nbr,1), ...
    [gamma_A_end;   phi_A_end;   R_star_end;   0; by_0], M_, options_);
ys_ref = ext_path_A_by_steadystate(zeros(M_.endo_nbr,1), ...
    [gamma_A_end;   phi_A_end;   R_star_end;   0; by_new], M_, options_);

iy  = strcmp(M_.endo_names,'y_A');
iff = strcmp(M_.endo_names,'f_A');

fy_target_A_start = ys_ini(iff)/ys_ini(iy);
fy_target_A_dem   = ys_dem(iff)/ys_dem(iy);   % demographic target (= baseline)
fy_target_A_ref   = ys_ref(iff)/ys_ref(iy);   % post-reform target (b_y = 0.9)

fprintf('fy_target: start=%.4f  dem=%.4f  ref=%.4f  (salto reforma = %+.4f)\n', ...
        fy_target_A_start, fy_target_A_dem, fy_target_A_ref, ...
        fy_target_A_ref - fy_target_A_dem);

oo_.exo_steady_state(strcmp(M_.exo_names,'fy_target_A')) = fy_target_A_start;

steady;
check;


%%%%%%% EXTENDED-PATH %%%%%%%

options_.noprint = 1;
perfect_foresight_setup(periods=400);

n_vars = M_.endo_nbr;

realized_path = zeros(n_vars, horizon_plot + 1);
id_gamma_A     = find(strcmp(M_.exo_names, 'gamma_A'));
id_phi_A       = find(strcmp(M_.exo_names, 'phi_A'));
id_R_star      = find(strcmp(M_.exo_names, 'R_star'));
id_fy_target_A = find(strcmp(M_.exo_names, 'fy_target_A'));    
id_by          = find(strcmp(M_.exo_names, 'b_y'));

% recorded exo paths (filled IN the loop; export uses these directly)
gamma_real     = zeros(horizon_plot+1,1);  gamma_real(1)     = gamma_A_start;
phi_real       = zeros(horizon_plot+1,1);  phi_real(1)       = phi_A_start;
R_star_real    = zeros(horizon_plot+1,1);  R_star_real(1)    = R_star_start;
fy_target_real = zeros(horizon_plot+1,1);  fy_target_real(1) = fy_target_A_start;
by_real        = zeros(horizon_plot+1,1);  by_real(1)        = by_0;

current_state = oo_.steady_state;
realized_path(:, 1) = current_state;

disp('>>> STARTING ROLLING MIT SHOCKS LOOP... <<<');

for i = 1:horizon_plot

    % ---- demographics (unchanged) ------------------------------------------
    if i <= steps
        current_gamma_A = gamma_A_data(i + 1);  % i=1 -> gamma_A_data(2) = 1971
        current_phi_A   = phi_A_data(i + 1);
        w_dem = (1 - cos(pi * i / steps)) / 2;
    else
        current_gamma_A = gamma_A_end;
        current_phi_A   = phi_A_end;
        w_dem = 1;
    end

    % ---- external target: demographic path + reform migration -------------
    fy_dem = fy_target_A_start + (fy_target_A_dem - fy_target_A_start) * w_dem;
    prog   = min(max((i - t_reform + 1)/reform_span, 0), 1);
    w_ref  = (1 - cos(pi * prog)) / 2;
    current_fy_target_A = fy_dem + (fy_target_A_ref - fy_target_A_dem) * w_ref;

    current_R_star = R_star_path(i + 1);
    current_by     = by_path(i + 1);

    % ---- record what is fed (export reads these; no drift possible) --------
    gamma_real(i+1)     = current_gamma_A;
    phi_real(i+1)       = current_phi_A;
    R_star_real(i+1)    = current_R_star;
    fy_target_real(i+1) = current_fy_target_A;
    by_real(i+1)        = current_by;

    oo_.exo_simul(:, id_gamma_A)     = current_gamma_A;
    oo_.exo_simul(:, id_phi_A)       = current_phi_A;
    oo_.exo_simul(:, id_R_star)      = current_R_star;
    oo_.exo_simul(:, id_fy_target_A) = current_fy_target_A;
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
fprintf('R_tilde_A_end - R_star_end = %.2e\n', ...
    realized_path(strcmp(M_.endo_names,'tilde_R_A'), end) - R_star_end);


%%%%%%% AUTOMATIC PLOTS %%%%%%%

get_series = @(name) realized_path(find(strcmp(M_.endo_names, name)), 1:horizon_plot);

vars_to_plot = {'y_A', 'k_A', 'f_A', 'tilde_R_A', ... 
                'c_y_A', 'c_m_A', 'c_o_A', ...
                'eta_y_A', 'eta_m_A', 'eta_o_A', ...
                'a_y_A', 'a_m_A', 'a_o_A'};

% Figure 1: Macro
figure('Name', 'Macro Aggregates', 'Position', [100 100 1000 600]);
subplot(2,2,1); plot(0:horizon_plot-1, get_series('y_A'), 'LineWidth', 1.5); title('Output (y\_A)'); grid on;
subplot(2,2,2); plot(0:horizon_plot-1, get_series('k_A'), 'LineWidth', 1.5); title('Capital (k\_A)'); grid on;
subplot(2,2,3); plot(0:horizon_plot-1, get_series('f_A'), 'LineWidth', 1.5); title('NFA (f\_A)'); grid on;
subplot(2,2,4); plot(0:horizon_plot-1, get_series('tilde_R_A'), 'LineWidth', 1.5); title('Interest Rate'); grid on;

% Figure 2: Consumption by cohort
figure('Name', 'Consumption by Cohort', 'Position', [100 100 1000 400]);
plot(0:horizon_plot-1, get_series('c_y_A'), 'b', 'LineWidth', 1.5); hold on;
plot(0:horizon_plot-1, get_series('c_m_A'), 'r', 'LineWidth', 1.5);
plot(0:horizon_plot-1, get_series('c_o_A'), 'k', 'LineWidth', 1.5);
legend('Young', 'Middle', 'Old'); title('Consumption'); grid on;

% Figure 3: Wealth shares
figure('Name', 'Wealth Shares', 'Position', [100 100 1000 400]);
plot(0:horizon_plot-1, get_series('eta_y_A'), 'b', 'LineWidth', 1.5); hold on;
plot(0:horizon_plot-1, get_series('eta_m_A'), 'r', 'LineWidth', 1.5);
plot(0:horizon_plot-1, get_series('eta_o_A'), 'k', 'LineWidth', 1.5);
legend('Young', 'Middle', 'Old'); title('\eta (wealth shares)'); grid on;

for f = 1:3
    figure(f);
    saveas(gcf, sprintf('output_A/aging_A_fig_%d.png', f));  
end


%%%%%%% EXPORTING DATA TO CSV %%%%%%%
disp('>>> EXPORTING FULL DUMP... <<<');

N = horizon_plot + 1;   % t = 0 (SS inicial) ate t = horizon_plot

% All endogeneous
T_full = array2table(realized_path(:, 1:N)', 'VariableNames', M_.endo_names');
T_full = addvars(T_full, (0:horizon_plot)', 'Before', 1, 'NewVariableNames', 'Time');

% All exogenous — straight from the recorded paths (no recompute, no drift)
T_full = addvars(T_full, gamma_real, phi_real, R_star_real, fy_target_real, by_real, ...
    'NewVariableNames', {'gamma_A', 'phi_A', 'R_star', 'fy_target_A', 'b_y'});

% All (constant) parameters
T_full = [T_full, array2table(repmat(M_.params', N, 1), 'VariableNames', M_.param_names')];

% Variables in efficiency units
T_full.w_indiv_A = T_full.w_A  ./ T_full.tilde_s_m_A;
T_full.t_indiv_A = T_full.t_A  ./ T_full.tilde_s_m_A;
T_full.e_indiv_A = T_full.e_A  ./ T_full.tilde_s_o_A;
T_full.c_y_indiv_A = T_full.c_y_A ./ T_full.tilde_s_y_A;
T_full.c_m_indiv_A = T_full.c_m_A ./ T_full.tilde_s_m_A;
T_full.c_o_indiv_A = T_full.c_o_A ./ T_full.tilde_s_o_A;
T_full.a_y_indiv_A = T_full.a_y_A ./ T_full.tilde_s_y_A;
T_full.a_m_indiv_A = T_full.a_m_A ./ T_full.tilde_s_m_A;
T_full.a_o_indiv_A = T_full.a_o_A ./ T_full.tilde_s_o_A;
T_full(:, startsWith(T_full.Properties.VariableNames,'AUX')) = [];

writetable(T_full, 'output_A_by/transition_results_A_by.csv');
fprintf('CSV: %d rows x %d columns\n', height(T_full), width(T_full));
