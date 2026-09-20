% addpath '/Applications/Dynare/7.0-arm64/matlab' - Paste this in the terminal
% ============================================================================
% PENSION FREEZE REFORM (wage de-indexation), unanticipated at t_reform.
%   e_B = (repl * w_B/tilde_s_m_B + pen_bar_B) * tilde_s_o_B
%   pre-reform : repl = repl_0, pen_bar_B = 0            (old rule, wage-indexed)
%   post-reform: repl = 0,      pen_bar_B = frozen level (computed ON THE FLY at
%                t_reform from the realized t_reform-1 pension; no second pass)
% freeze_real = true  -> pen_bar declines at rate x (price freeze). WARNING: its
% infinite-horizon BGP has zero pensions (degenerate); default false.
% Also fixed vs previous version: horizon_plot defined BEFORE use; exo paths
% recorded inside the loop (export can no longer drift from what was fed).
% ============================================================================

%%%%%%% VARIABLES %%%%%%%
var

% Endogenous
n_B 
psi_o_B
psi_y_B
n_y_B
n_m_B
n_o_B
mu_o_B
mu_m_B
mu_y_B
s_o_B
s_m_B
s_y_B
h_m_B
h_y_B
a_B
a_o_B
a_m_B
a_y_B
z_o_B
z_m_B
z_y_B
c_o_B
c_m_B
c_y_B
y_B
k_B
w_B
inv_B
b_B
gov_B
e_B
t_B
f_B
nx_B
tilde_R_B
g_B
q_B

% Auxiliary (endogenous)
tilde_s_m_B
tilde_s_y_B
tilde_s_o_B
Omega_m_B
Omega_y_B
eta_o_B
eta_m_B
eta_y_B
walras_B
;

varexo gamma_B, phi_B, R_star, fy_target_B, repl, pen_bar_B;


%%%%%%% PARAMETERS %%%%%%%

parameters alpha, delta, beta, sigma, b_y, g_y, x, omega, nu, Lambda_B;

omega = 1 - 1 / 40;     % worklife = 40
nu = 1 - 1 / 5;     % university = 5
alpha = 0.67;
delta = 0.1;
beta = 1.025;
sigma = 0.5;
b_y = 0.6;
g_y = 0.25;
x = 0.005;
Lambda_B = 0.1;


%%%%%%% EQUATIONS %%%%%%%

model;

% 0. Normalization
g_B = (1 + n_B) * (1 + x);
q_B = tilde_R_B / g_B;

% 1. Demographics
n_y_B = (nu + phi_B * 1 / psi_y_B(-1)) - 1;
n_m_B = (omega + (1 - nu) * psi_y_B(-1)) - 1;
n_o_B = ((1 - omega) * 1 / psi_o_B(-1) + gamma_B) - 1;
n_B = phi_B * tilde_s_m_B(-1) - (1 - gamma_B) * tilde_s_o_B(-1);
psi_y_B = (nu * psi_y_B(-1) + phi_B) * 1 / (1 + n_m_B);
psi_o_B = (1 - omega + gamma_B * psi_o_B(-1)) * 1 / (1 + n_m_B);
tilde_s_y_B = (nu * tilde_s_y_B(-1) + phi_B * tilde_s_m_B(-1)) * 1 / (1 + n_B);
tilde_s_m_B = (omega * tilde_s_m_B(-1) + (1 - nu) * tilde_s_y_B(-1)) * 1 / (1 + n_B);
tilde_s_o_B = ((1 - omega) * tilde_s_m_B(-1) + gamma_B * tilde_s_o_B(-1)) * 1 / (1 + n_B); 

% 2. Household
%% MPCs
1 = (1 + gamma_B(+1) * beta^sigma * tilde_R_B(+1)^(sigma - 1) * (1 / mu_o_B(+1))) * mu_o_B; 
1 = (1 + beta^sigma * (Omega_m_B(+1) * tilde_R_B(+1))^(sigma - 1) * (1 / mu_m_B(+1))) * mu_m_B;
1 = (1 + beta^sigma * (Omega_y_B(+1) * tilde_R_B(+1))^(sigma - 1) * (1 / mu_y_B(+1))) * mu_y_B;
Omega_m_B = omega + (1 - omega) * (mu_o_B / mu_m_B)^(1 / (1 - sigma));
Omega_y_B = nu + (1 - nu) * (mu_m_B / mu_y_B)^(1 / (1 - sigma));

%% PDVs
s_o_B = e_B + gamma_B(+1) * g_B(+1) / (tilde_R_B(+1) * (1 + n_o_B(+1))) * s_o_B(+1); 
s_m_B = g_B(+1) / (Omega_m_B(+1) * tilde_R_B(+1)) * ((Omega_m_B(+1) - omega(+1)) / ((1 + n_o_B(+1)) * psi_o_B) * s_o_B(+1) + omega(+1) / (1 + n_m_B(+1)) * s_m_B(+1)); 
s_y_B = g_B(+1) / (Omega_y_B(+1) * tilde_R_B(+1)) * ((Omega_y_B(+1) - nu(+1)) * psi_y_B / (1 + n_m_B(+1)) * s_m_B(+1) + nu(+1) / (1 + n_y_B(+1)) * s_y_B(+1)); 
h_m_B = w_B - t_B + g_B(+1) * omega(+1) / (Omega_m_B(+1) * tilde_R_B(+1) * (1 + n_m_B(+1))) * h_m_B(+1); 
h_y_B = g_B(+1) / (Omega_y_B(+1) * tilde_R_B(+1)) * ((Omega_y_B(+1) - nu(+1)) * psi_y_B / (1 + n_m_B(+1)) * h_m_B(+1) + nu(+1) / (1 + n_y_B(+1)) * h_y_B(+1));

% Assets
a_o_B = q_B * (1 - mu_o_B) * (eta_o_B(-1) + (1 - omega) * eta_m_B(-1)) * a_B(-1) + e_B - mu_o_B * s_o_B;
a_m_B = q_B * (1 - mu_m_B) * (omega * eta_m_B(-1) + (1 - nu) * eta_y_B(-1)) * a_B(-1) + w_B - t_B - mu_m_B * (h_m_B + s_m_B);
a_y_B = q_B * (1 - mu_y_B) * nu * eta_y_B(-1) * a_B(-1) - mu_y_B * (h_y_B + s_y_B);
eta_o_B = a_o_B / a_B;
eta_m_B = a_m_B / a_B;
eta_y_B = a_y_B / a_B;
eta_y_B + eta_m_B + eta_o_B = 1;

% Wealth
z_o_B = q_B * (eta_o_B(-1) + (1 - omega) * eta_m_B(-1)) * a_B(-1) + s_o_B;
z_m_B = q_B * (omega * eta_m_B(-1) + (1 - nu) * eta_y_B(-1)) * a_B(-1) + h_m_B + s_m_B;
z_y_B = q_B * nu * eta_y_B(-1) * a_B(-1) + h_y_B + s_y_B;

% Consumption
c_o_B = mu_o_B * z_o_B;
c_m_B = mu_m_B * z_m_B;
c_y_B = mu_y_B * z_y_B;


% 3. Firm
tilde_R_B(+1) = (1 - alpha) * y_B(+1) / k_B * g_B(+1) + (1 - delta);
y_B = tilde_s_m_B^alpha * (k_B(-1) / g_B)^(1 - alpha);
w_B = alpha * y_B;
inv_B = k_B - (1 - delta) / g_B * k_B(-1);


% 4. Government
b_B = b_y * y_B;
gov_B = g_y * y_B;
e_B = (repl * w_B / tilde_s_m_B + pen_bar_B) * tilde_s_o_B;
t_B = e_B + gov_B + q_B * b_B(-1) - b_B;


% 5. Balance of Payments
a_B = k_B + b_B + f_B;
nx_B = f_B - q_B * f_B(-1);


% 6. Imperfect Capital Mobility
tilde_R_B = R_star + Lambda_B * (exp(-(f_B(-1)/y_B(-1) - fy_target_B)) - 1);


% Walras check
walras_B = y_B - c_o_B - c_m_B - c_y_B - inv_B - gov_B - nx_B;

end;


%%%%%%% CALIBRATION %%%%%%%

D = load('demo_paths_B.mat');   
R_data = load('R_star_path.mat'); 
R_star_path = R_data.R_star;    

y_start = 1970;
y_end   = 2070;

% horizon defined HERE (was used-before-defined in the old file)
steps        = 100;
horizon_plot = 200;

i0 = find(D.Year == y_start);
i1 = find(D.Year == y_end);
assert(~isempty(i0) && ~isempty(i1), 'Period outside the interval of Year');

Year_win        = D.Year(i0:i1);            
gamma_B_data    = D.gamma_B_path(i0:i1);    
phi_B_data      = D.phi_B_path(i0:i1);

gamma_B_start = gamma_B_data(1);            
gamma_B_end   = gamma_B_data(end);          
phi_B_start   = phi_B_data(1);
phi_B_end     = phi_B_data(end);
R_star_start = R_star_path(1);      
R_star_end   = R_star_path(end);

% ---- REFORM: pension freeze (unanticipated wage de-indexation) -------------
repl_0      = 0.45;      % pre-reform replacement rate (wage-indexed rule)
t_reform    = 55;        % onset (1970+55 = 2025); pension frozen at 2024 level
reform_span = 40;        % migration of the EXTERNAL TARGET only (cosine)
freeze_real = true;     % true: pen_bar erodes at rate x (price freeze).
                         % WARNING: real-freeze BGP has zero pensions.

initval;                     
  gamma_B   = gamma_B_start;   
  phi_B     = phi_B_start;
  R_star    = R_star_start; 
  repl      = repl_0;
  pen_bar_B = 0;
end;

% Initial (exogeneous) BGP
oo_.exo_steady_state(strcmp(M_.exo_names,'gamma_B'))   = gamma_B_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'phi_B'))     = phi_B_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'R_star'))    = R_star_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'repl'))      = repl_0;
oo_.exo_steady_state(strcmp(M_.exo_names,'pen_bar_B')) = 0;

% exo order assumed: [gamma; phi; R_star; fy_target; repl; pen_bar]
% (confirm with disp(M_.exo_names) after first compile)
ys_ini = ext_path_B_repl_steadystate(zeros(M_.endo_nbr,1), ...
    [gamma_B_start; phi_B_start; R_star_start; 0; repl_0; 0], M_, options_);
ys_dem = ext_path_B_repl_steadystate(zeros(M_.endo_nbr,1), ...
    [gamma_B_end;   phi_B_end;   R_star_end;   0; repl_0; 0], M_, options_);

iy  = strcmp(M_.endo_names,'y_B');
iff = strcmp(M_.endo_names,'f_B');
ie  = strcmp(M_.endo_names,'e_B');
iso = strcmp(M_.endo_names,'tilde_s_o_B');

fy_target_B_start = ys_ini(iff)/ys_ini(iy);
fy_target_B_dem   = ys_dem(iff)/ys_dem(iy);   % demographic target (= baseline)
fy_target_B_ref   = fy_target_B_dem;          % placeholder; set at t_reform
pen_frozen        = NaN;                      % set on the fly at t_reform

oo_.exo_steady_state(strcmp(M_.exo_names,'fy_target_B')) = fy_target_B_start;

steady;
check;


%%%%%%% EXTENDED-PATH %%%%%%%

options_.noprint = 1;
perfect_foresight_setup(periods=400);

n_vars = M_.endo_nbr;

realized_path = zeros(n_vars, horizon_plot + 1);
id_gamma_B     = find(strcmp(M_.exo_names, 'gamma_B'));
id_phi_B       = find(strcmp(M_.exo_names, 'phi_B'));
id_R_star      = find(strcmp(M_.exo_names, 'R_star'));
id_fy_target_B = find(strcmp(M_.exo_names, 'fy_target_B')); 
id_repl        = find(strcmp(M_.exo_names, 'repl'));   
id_pen         = find(strcmp(M_.exo_names, 'pen_bar_B'));

% recorded exo paths (filled IN the loop; export uses these directly)
gamma_real     = zeros(horizon_plot+1,1);  gamma_real(1)     = gamma_B_start;
phi_real       = zeros(horizon_plot+1,1);  phi_real(1)       = phi_B_start;
R_star_real    = zeros(horizon_plot+1,1);  R_star_real(1)    = R_star_start;
fy_target_real = zeros(horizon_plot+1,1);  fy_target_real(1) = fy_target_B_start;
repl_real      = zeros(horizon_plot+1,1);  repl_real(1)      = repl_0;
pen_real       = zeros(horizon_plot+1,1);  pen_real(1)       = 0;

current_state = oo_.steady_state;
realized_path(:, 1) = current_state;

disp('>>> STARTING ROLLING MIT SHOCKS LOOP... <<<');

for i = 1:horizon_plot

    % ---- demographics (unchanged) ------------------------------------------
    if i <= steps
        current_gamma_B = gamma_B_data(i + 1);  % i=1 -> gamma_B_data(2) = 1971
        current_phi_B   = phi_B_data(i + 1);
        w_dem = (1 - cos(pi * i / steps)) / 2;
    else
        current_gamma_B = gamma_B_end;
        current_phi_B   = phi_B_end;
        w_dem = 1;
    end

    % ---- REFORM: freeze the pension at the realized t_reform-1 level -------
    if i == t_reform
        pen_frozen = current_state(ie) / current_state(iso);   % e_indiv at 2024
        ys_ref = ext_path_B_repl_steadystate(zeros(M_.endo_nbr,1), ...
            [gamma_B_end; phi_B_end; R_star_end; 0; 0; pen_frozen], M_, options_);
        fy_target_B_ref = ys_ref(iff)/ys_ref(iy);
        fprintf('>>> REFORM at t=%d: pen frozen at %.5f | fy_ref=%.4f (dem %.4f, var %.4f)\n', ...
            t_reform, pen_frozen, fy_target_B_ref, fy_target_B_dem, ...
            fy_target_B_ref - fy_target_B_dem);
    end

    if i < t_reform
        current_repl = repl_0;
        current_pen  = 0;
    else
        current_repl = 0;                    % wage link severed
        if freeze_real
            current_pen = pen_frozen / (1 + x)^(i - t_reform);  % price freeze
        else
            current_pen = pen_frozen;        % efficiency-unit freeze (de-indexation)
        end
    end

    % ---- external target: demographic path + reform migration -------------
    fy_dem = fy_target_B_start + (fy_target_B_dem - fy_target_B_start) * w_dem;
    prog   = min(max((i - t_reform + 1)/reform_span, 0), 1);
    w_ref  = (1 - cos(pi * prog)) / 2;
    current_fy_target_B = fy_dem + (fy_target_B_ref - fy_target_B_dem) * w_ref;

    current_R_star = R_star_path(i + 1);

    % ---- record what is fed (export reads these; no drift possible) --------
    gamma_real(i+1)     = current_gamma_B;
    phi_real(i+1)       = current_phi_B;
    R_star_real(i+1)    = current_R_star;
    fy_target_real(i+1) = current_fy_target_B;
    repl_real(i+1)      = current_repl;
    pen_real(i+1)       = current_pen;

    oo_.exo_simul(:, id_gamma_B)     = current_gamma_B;
    oo_.exo_simul(:, id_phi_B)       = current_phi_B;
    oo_.exo_simul(:, id_R_star)      = current_R_star;
    oo_.exo_simul(:, id_fy_target_B) = current_fy_target_B;
    oo_.exo_simul(:, id_repl)        = current_repl;
    oo_.exo_simul(:, id_pen)         = current_pen;
    
    oo_.endo_simul = repmat(current_state, 1, options_.periods + 2);
    
    perfect_foresight_solver;
    
    current_state = oo_.endo_simul(:, 2);
    realized_path(:, i+1) = current_state;
    
    fprintf('Year %d of %d (calendar %d) solved under constant beliefs.\n', ...
            i, horizon_plot, y_start + i);
end
disp('>>> THE LOOP IS FINISHED! <<<');


%% Check for PCM in BGP
fprintf('R_tilde_B_end - R_star_end = %.2e\n', ...
    realized_path(strcmp(M_.endo_names,'tilde_R_B'), end) - R_star_end);


%%%%%%% AUTOMATIC PLOTS %%%%%%%

get_series = @(name) realized_path(find(strcmp(M_.endo_names, name)), 1:horizon_plot);

vars_to_plot = {'y_B', 'k_B', 'f_B', 'tilde_R_B', ... 
                'c_y_B', 'c_m_B', 'c_o_B', ...
                'eta_y_B', 'eta_m_B', 'eta_o_B', ...
                'a_y_B', 'a_m_B', 'a_o_B'};

% Figure 1: Macro
figure('Name', 'Macro Aggregates', 'Position', [100 100 1000 600]);
subplot(2,2,1); plot(0:horizon_plot-1, get_series('y_B'), 'LineWidth', 1.5); title('Output (y\_B)'); grid on;
subplot(2,2,2); plot(0:horizon_plot-1, get_series('k_B'), 'LineWidth', 1.5); title('Capital (k\_B)'); grid on;
subplot(2,2,3); plot(0:horizon_plot-1, get_series('f_B'), 'LineWidth', 1.5); title('NFA (f\_B)'); grid on;
subplot(2,2,4); plot(0:horizon_plot-1, get_series('tilde_R_B'), 'LineWidth', 1.5); title('Interest Rate'); grid on;

% Figure 2: Consumption by cohort
figure('Name', 'Consumption by Cohort', 'Position', [100 100 1000 400]);
plot(0:horizon_plot-1, get_series('c_y_B'), 'b', 'LineWidth', 1.5); hold on;
plot(0:horizon_plot-1, get_series('c_m_B'), 'r', 'LineWidth', 1.5);
plot(0:horizon_plot-1, get_series('c_o_B'), 'k', 'LineWidth', 1.5);
legend('Young', 'Middle', 'Old'); title('Consumption'); grid on;

% Figure 3: Wealth shares
figure('Name', 'Wealth Shares', 'Position', [100 100 1000 400]);
plot(0:horizon_plot-1, get_series('eta_y_B'), 'b', 'LineWidth', 1.5); hold on;
plot(0:horizon_plot-1, get_series('eta_m_B'), 'r', 'LineWidth', 1.5);
plot(0:horizon_plot-1, get_series('eta_o_B'), 'k', 'LineWidth', 1.5);
legend('Young', 'Middle', 'Old'); title('\eta (wealth shares)'); grid on;

for f = 1:3
    figure(f);
    saveas(gcf, sprintf('output_B/aging_B_fig_%d.png', f));  
end


%%%%%%% EXPORTING DATA TO CSV %%%%%%%
disp('>>> EXPORTING FULL DUMP... <<<');

N = horizon_plot + 1;   % t = 0 (SS inicial) ate t = horizon_plot

% All endogeneous
T_full = array2table(realized_path(:, 1:N)', 'VariableNames', M_.endo_names');
T_full = addvars(T_full, (0:horizon_plot)', 'Before', 1, 'NewVariableNames', 'Time');

% All exogenous — straight from the recorded paths (no recompute, no drift)
T_full = addvars(T_full, gamma_real, phi_real, R_star_real, fy_target_real, ...
                 repl_real, pen_real, ...
    'NewVariableNames', {'gamma_B', 'phi_B', 'R_star', 'fy_target_B', ...
                         'repl', 'pen_bar_B'});

% All (constant) parameters
T_full = [T_full, array2table(repmat(M_.params', N, 1), 'VariableNames', M_.param_names')];

% Variables in efficiency units
T_full.w_indiv_B = T_full.w_B  ./ T_full.tilde_s_m_B;
T_full.t_indiv_B = T_full.t_B  ./ T_full.tilde_s_m_B;
T_full.e_indiv_B = T_full.e_B  ./ T_full.tilde_s_o_B;
T_full.c_y_indiv_B = T_full.c_y_B ./ T_full.tilde_s_y_B;
T_full.c_m_indiv_B = T_full.c_m_B ./ T_full.tilde_s_m_B;
T_full.c_o_indiv_B = T_full.c_o_B ./ T_full.tilde_s_o_B;
T_full.a_y_indiv_B = T_full.a_y_B ./ T_full.tilde_s_y_B;
T_full.a_m_indiv_B = T_full.a_m_B ./ T_full.tilde_s_m_B;
T_full.a_o_indiv_B = T_full.a_o_B ./ T_full.tilde_s_o_B;
T_full(:, startsWith(T_full.Properties.VariableNames,'AUX')) = [];

writetable(T_full, 'output_B_repl/transition_results_B_repl.csv');
fprintf('CSV: %d rows x %d columns\n', height(T_full), width(T_full));
