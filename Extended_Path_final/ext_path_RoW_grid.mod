% addpath '/Applications/Dynare/7.0-arm64/matlab' - Paste this in the terminal

%%%%%%% VARIABLES %%%%%%%
var

% Endogenous
n_RoW 
psi_o_RoW
psi_y_RoW
n_y_RoW
n_m_RoW
n_o_RoW
mu_o_RoW
mu_m_RoW
mu_y_RoW
s_o_RoW
s_m_RoW
s_y_RoW
h_m_RoW
h_y_RoW
a_RoW
a_o_RoW
a_m_RoW
a_y_RoW
z_o_RoW
z_m_RoW
z_y_RoW
c_o_RoW
c_m_RoW
c_y_RoW
y_RoW
k_RoW
w_RoW
inv_RoW
b_RoW
gov_RoW
e_RoW
t_RoW
tilde_R_RoW
g_RoW
q_RoW

% Auxiliary (endogenous)
tilde_s_m_RoW
tilde_s_y_RoW
tilde_s_o_RoW
Omega_m_RoW
Omega_y_RoW
eta_o_RoW
eta_m_RoW
eta_y_RoW
walras_RoW
;

varexo gamma_RoW, phi_RoW;


%%%%%%% PARAMETERS %%%%%%%

parameters alpha, delta, beta, sigma, b_y, g_y, repl, x, omega, nu;

omega = 1 - 1 / 40;     % worklife = 40
nu = 1 - 1 / 5;     % university = 5
alpha = 0.67;
delta = 0.1;
@#include "beta_setting.mod"
b_y = 0.6;
g_y = 0.25;
repl = 0.40;
x = 0.005;


%%%%%%% EQUATIONS %%%%%%%

model;

% 0. Normalization
g_RoW = (1 + n_RoW) * (1 + x);
q_RoW = tilde_R_RoW / g_RoW;

% 1. Demographics
n_y_RoW = (nu + phi_RoW * 1 / psi_y_RoW(-1)) - 1;
n_m_RoW = (omega + (1 - nu) * psi_y_RoW(-1)) - 1;
n_o_RoW = ((1 - omega) * 1 / psi_o_RoW(-1) + gamma_RoW) - 1;
n_RoW = phi_RoW * tilde_s_m_RoW(-1) - (1 - gamma_RoW) * tilde_s_o_RoW(-1);
psi_y_RoW = (nu * psi_y_RoW(-1) + phi_RoW) * 1 / (1 + n_m_RoW);
psi_o_RoW = (1 - omega + gamma_RoW * psi_o_RoW(-1)) * 1 / (1 + n_m_RoW);
tilde_s_y_RoW = (nu * tilde_s_y_RoW(-1) + phi_RoW * tilde_s_m_RoW(-1)) * 1 / (1 + n_RoW);
tilde_s_m_RoW = (omega * tilde_s_m_RoW(-1) + (1 - nu) * tilde_s_y_RoW(-1)) * 1 / (1 + n_RoW);
tilde_s_o_RoW = ((1 - omega) * tilde_s_m_RoW(-1) + gamma_RoW * tilde_s_o_RoW(-1)) * 1 / (1 + n_RoW); 

% 2. Household
%% MPCs
1 = (1 + gamma_RoW(+1) * beta^sigma * tilde_R_RoW(+1)^(sigma - 1) * (1 / mu_o_RoW(+1))) * mu_o_RoW; 
1 = (1 + beta^sigma * (Omega_m_RoW(+1) * tilde_R_RoW(+1))^(sigma - 1) * (1 / mu_m_RoW(+1))) * mu_m_RoW;
1 = (1 + beta^sigma * (Omega_y_RoW(+1) * tilde_R_RoW(+1))^(sigma - 1) * (1 / mu_y_RoW(+1))) * mu_y_RoW;
Omega_m_RoW = omega + (1 - omega) * (mu_o_RoW / mu_m_RoW)^(1 / (1 - sigma));
Omega_y_RoW = nu + (1 - nu) * (mu_m_RoW / mu_y_RoW)^(1 / (1 - sigma));

%% PDVs
s_o_RoW = e_RoW + gamma_RoW(+1) * g_RoW(+1) / (tilde_R_RoW(+1) * (1 + n_o_RoW(+1))) * s_o_RoW(+1); 
s_m_RoW = g_RoW(+1) / (Omega_m_RoW(+1) * tilde_R_RoW(+1)) * ((Omega_m_RoW(+1) - omega(+1)) / ((1 + n_o_RoW(+1)) * psi_o_RoW) * s_o_RoW(+1) + omega(+1) / (1 + n_m_RoW(+1)) * s_m_RoW(+1)); 
s_y_RoW = g_RoW(+1) / (Omega_y_RoW(+1) * tilde_R_RoW(+1)) * ((Omega_y_RoW(+1) - nu(+1)) * psi_y_RoW / (1 + n_m_RoW(+1)) * s_m_RoW(+1) + nu(+1) / (1 + n_y_RoW(+1)) * s_y_RoW(+1)); 
h_m_RoW = w_RoW - t_RoW + g_RoW(+1) * omega(+1) / (Omega_m_RoW(+1) * tilde_R_RoW(+1) * (1 + n_m_RoW(+1))) * h_m_RoW(+1); 
h_y_RoW = g_RoW(+1) / (Omega_y_RoW(+1) * tilde_R_RoW(+1)) * ((Omega_y_RoW(+1) - nu(+1)) * psi_y_RoW / (1 + n_m_RoW(+1)) * h_m_RoW(+1) + nu(+1) / (1 + n_y_RoW(+1)) * h_y_RoW(+1));

% Assets
a_o_RoW = q_RoW * (1 - mu_o_RoW) * (eta_o_RoW(-1) + (1 - omega) * eta_m_RoW(-1)) * a_RoW(-1) + e_RoW - mu_o_RoW * s_o_RoW;
a_m_RoW = q_RoW * (1 - mu_m_RoW) * (omega * eta_m_RoW(-1) + (1 - nu) * eta_y_RoW(-1)) * a_RoW(-1) + w_RoW - t_RoW - mu_m_RoW * (h_m_RoW + s_m_RoW);
a_y_RoW = q_RoW * (1 - mu_y_RoW) * nu * eta_y_RoW(-1) * a_RoW(-1) - mu_y_RoW * (h_y_RoW + s_y_RoW);
eta_o_RoW = a_o_RoW / a_RoW;
eta_m_RoW = a_m_RoW / a_RoW;
eta_y_RoW = a_y_RoW / a_RoW;
eta_y_RoW + eta_m_RoW + eta_o_RoW = 1;

% Wealth
z_o_RoW = q_RoW * (eta_o_RoW(-1) + (1 - omega) * eta_m_RoW(-1)) * a_RoW(-1) + s_o_RoW;
z_m_RoW = q_RoW * (omega * eta_m_RoW(-1) + (1 - nu) * eta_y_RoW(-1)) * a_RoW(-1) + h_m_RoW + s_m_RoW;
z_y_RoW = q_RoW * nu * eta_y_RoW(-1) * a_RoW(-1) + h_y_RoW + s_y_RoW;

% Consumption
c_o_RoW = mu_o_RoW * z_o_RoW;
c_m_RoW = mu_m_RoW * z_m_RoW;
c_y_RoW = mu_y_RoW * z_y_RoW;


% 3. Firm
tilde_R_RoW = (1 - alpha) * y_RoW / k_RoW(-1) * g_RoW + (1 - delta);  
y_RoW = tilde_s_m_RoW^alpha * (k_RoW(-1) / g_RoW)^(1 - alpha);
w_RoW = alpha * y_RoW;
inv_RoW = k_RoW - (1 - delta) / g_RoW * k_RoW(-1);


% 4. Government
b_RoW = b_y * y_RoW;
gov_RoW = g_y * y_RoW;
e_RoW = repl * w_RoW / tilde_s_m_RoW * tilde_s_o_RoW;
t_RoW = e_RoW + gov_RoW + q_RoW * b_RoW(-1) - b_RoW;


% 5. Assets Market Clearing
a_RoW = k_RoW + b_RoW;


% Walras check
walras_RoW = y_RoW - c_o_RoW - c_m_RoW - c_y_RoW - inv_RoW - gov_RoW;

end;


%%%%%%% CALIBRATION %%%%%%%

D = load('demo_paths_RoW.mat');    % Year, gamma_RoW_path, phi_RoW_path
y_start = 1970;
y_end   = 2070;

i0 = find(D.Year == y_start);
i1 = find(D.Year == y_end);
assert(~isempty(i0) && ~isempty(i1), 'Ano fora do intervalo do .mat');

Year_win          = D.Year(i0:i1);              % 1970..2070 (101 anos)
gamma_RoW_data    = D.gamma_RoW_path(i0:i1);    % vetor 101x1
phi_RoW_data      = D.phi_RoW_path(i0:i1);

gamma_RoW_start = gamma_RoW_data(1);            % 1970
gamma_RoW_end   = gamma_RoW_data(end);          % 2070
phi_RoW_start   = phi_RoW_data(1);
phi_RoW_end     = phi_RoW_data(end);

fprintf('WPP: gamma %.4f -> %.4f | phi %.4f -> %.4f\n', ...
        gamma_RoW_start, gamma_RoW_end, phi_RoW_start, phi_RoW_end);

initval;                     
  gamma_RoW = gamma_RoW_start; 
  phi_RoW   = phi_RoW_start;
end;

% Initial (exogenous) BGP
oo_.exo_steady_state(strcmp(M_.exo_names,'gamma_RoW')) = gamma_RoW_start;
oo_.exo_steady_state(strcmp(M_.exo_names,'phi_RoW'))   = phi_RoW_start;

steady;
check;

%%%%%%% EXTENDED-PATH %%%%%%%
options_.noprint = 1;
perfect_foresight_setup(periods=400);

steps = 100;              
horizon_plot = 100;       
n_vars = M_.endo_nbr;

realized_path = zeros(n_vars, horizon_plot + 1);

id_gamma_RoW = find(strcmp(M_.exo_names, 'gamma_RoW'));
id_phi_RoW   = find(strcmp(M_.exo_names, 'phi_RoW'));

current_state = oo_.steady_state;
realized_path(:, 1) = current_state;

disp('>>> STARTING ROLLING MIT SHOCKS LOOP... <<<');
for i = 1:horizon_plot
    
    if i <= steps
        % Atualiza apenas a demografia usando os dados da WPP
        current_gamma_RoW = gamma_RoW_data(i + 1);  
        current_phi_RoW   = phi_RoW_data(i + 1);
    else
        current_gamma_RoW = gamma_RoW_end;
        current_phi_RoW   = phi_RoW_end;
    end

    oo_.exo_simul(:, id_gamma_RoW) = current_gamma_RoW;
    oo_.exo_simul(:, id_phi_RoW)   = current_phi_RoW;
    
    oo_.endo_simul = repmat(current_state, 1, options_.periods + 2);
    
    perfect_foresight_solver;
    
    current_state = oo_.endo_simul(:, 2);
    realized_path(:, i+1) = current_state;
    
    fprintf('Year %d of %d (calendar %d) solved under constant beliefs.\n', ...
            i, horizon_plot, y_start + i);
end
disp('>>> THE LOOP IS FINISHED! <<<');