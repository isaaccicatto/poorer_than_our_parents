function [ys,params,check] = ext_path_A_steadystate(ys,exo,M_,options_)
check = 0;
alpha  = M_.params(strcmp(M_.param_names,'alpha'));
delta  = M_.params(strcmp(M_.param_names,'delta'));
beta   = M_.params(strcmp(M_.param_names,'beta'));
sigma  = M_.params(strcmp(M_.param_names,'sigma'));
g_y    = M_.params(strcmp(M_.param_names,'g_y'));
repl   = M_.params(strcmp(M_.param_names,'repl'));
x      = M_.params(strcmp(M_.param_names,'x'));
omega  = M_.params(strcmp(M_.param_names,'omega'));
nu     = M_.params(strcmp(M_.param_names,'nu'));
% exo order: [gamma; phi; R_star; fy_target; b_y]
gamma_A = exo(1); phi_A = exo(2); R_star = exo(3); b_y = exo(5);
% 1. Demographics
n_A = ((nu + omega) + sqrt((nu - omega)^2 + 4 * phi_A * (1 - nu))) / 2 - 1;
psi_y_A = phi_A / (1 + n_A - nu);
psi_o_A = (1 - omega) / (1 + n_A - gamma_A);
tilde_s_m_A = 1 / (1 + psi_y_A + psi_o_A);
tilde_s_y_A = psi_y_A * tilde_s_m_A;
tilde_s_o_A = psi_o_A * tilde_s_m_A;
n_y_A = n_A;
n_m_A = n_A;
n_o_A = n_A;
% 2. Normalization
tilde_R_A = R_star;
g_A = (1 + n_A) * (1 + x);
q_A = tilde_R_A / g_A;
% 3. Firm
kappa_A = (1 - alpha) * g_A / (tilde_R_A - (1 - delta));
y_A = (tilde_s_m_A^alpha * (kappa_A / g_A)^(1 - alpha))^(1 / alpha);
k_A = kappa_A * y_A;
w_A = alpha * y_A;
inv_A = (1 - (1 - delta) / g_A) * k_A;
% 4. Government
b_A = b_y * y_A;
gov_A = g_y * y_A;
e_A = repl * w_A / tilde_s_m_A * tilde_s_o_A;
t_A = e_A + gov_A - b_A * (1 - q_A);
% 5. Household
%% MPCs
mu_o_A = 1 - gamma_A * beta^sigma * tilde_R_A^(sigma - 1);
Om_m = @(m) omega + (1 - omega) * (mu_o_A / m)^(1 / (1 - sigma));
res_m = @(m) m - 1 + beta^sigma * (Om_m(m) * tilde_R_A)^(sigma - 1);
mu_m_A    = fzero(res_m, [1e-8, 1]);
Omega_m_A = Om_m(mu_m_A);
Om_y = @(m) nu + (1 - nu) * (mu_m_A / m)^(1 / (1 - sigma));
res_y = @(m) m - 1 + beta^sigma * (Om_y(m) * tilde_R_A)^(sigma - 1);
mu_y_A    = fzero(res_y, [1e-8, 1]);
Omega_y_A = Om_y(mu_y_A);
%% PDVs
s_o_A = tilde_R_A * e_A / (tilde_R_A - gamma_A * (1 + x));
s_m_A = (Omega_m_A - omega) / (Omega_m_A * tilde_R_A / (1 + x) - omega) * 1 / psi_o_A * s_o_A;
s_y_A = (Omega_y_A - nu) / (Omega_y_A * tilde_R_A / (1 + x) - nu) * psi_y_A * s_m_A;
h_m_A = Omega_m_A * tilde_R_A / (1 + x) / (Omega_m_A * tilde_R_A / (1 + x) - omega) * (w_A - t_A);
h_y_A = (Omega_y_A - nu) / (Omega_y_A * tilde_R_A / (1 + x) - nu) * psi_y_A * h_m_A;
% Assets
a_y_A = - mu_y_A * (h_y_A + s_y_A) / (1 - q_A * nu * (1 - mu_y_A));
a_m_A = (q_A * (1 - mu_m_A) * (1 - nu) * a_y_A + w_A - t_A - mu_m_A * (h_m_A + s_m_A)) / (1 - q_A * omega * (1 - mu_m_A));
a_o_A = (q_A * (1 - mu_o_A) * (1 - omega) * a_m_A + e_A - mu_o_A * s_o_A) / (1 - q_A * (1 - mu_o_A));
a_A   = a_y_A + a_m_A + a_o_A;
eta_y_A = a_y_A / a_A;
eta_m_A = a_m_A / a_A;
eta_o_A = a_o_A / a_A;
% Wealth
z_o_A = q_A * (eta_o_A + (1 - omega) * eta_m_A) * a_A + s_o_A;
z_m_A = q_A * (omega * eta_m_A + (1 - nu) * eta_y_A) * a_A + h_m_A + s_m_A;
z_y_A = q_A * nu * eta_y_A * a_A + h_y_A + s_y_A;
% Consumption
c_o_A = mu_o_A * z_o_A;
c_m_A = mu_m_A * z_m_A;
c_y_A = mu_y_A * z_y_A;
% 6. Balance of Payments
f_A = a_A - k_A - b_A;
nx_A = f_A * (1 - q_A);
fy_target_A = f_A / y_A;
% 7. Walras check
walras_A = y_A - c_o_A - c_m_A - c_y_A - inv_A - gov_A - nx_A;
for i=1:M_.endo_nbr
if strncmp(M_.endo_names{i},'AUX_EXO_LEAD',12)
        eval([M_.endo_names{i} ' = gamma_A;']);
end
end
params = M_.params;
for i=1:M_.endo_nbr, eval(['ys(' int2str(i) ') = ' M_.endo_names{i} ';']);
end