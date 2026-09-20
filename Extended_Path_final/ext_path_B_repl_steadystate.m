function [ys,params,check] = ext_path_B_repl_steadystate(ys,exo,M_,options_)
check = 0;
alpha  = M_.params(strcmp(M_.param_names,'alpha'));
delta  = M_.params(strcmp(M_.param_names,'delta'));
beta   = M_.params(strcmp(M_.param_names,'beta'));
sigma  = M_.params(strcmp(M_.param_names,'sigma'));
b_y    = M_.params(strcmp(M_.param_names,'b_y'));
g_y    = M_.params(strcmp(M_.param_names,'g_y'));
x      = M_.params(strcmp(M_.param_names,'x'));
omega  = M_.params(strcmp(M_.param_names,'omega'));
nu     = M_.params(strcmp(M_.param_names,'nu'));
% exo order: [gamma; phi; R_star; fy_target; repl; pen_bar]
gamma_B = exo(1); phi_B = exo(2); R_star = exo(3); repl = exo(5); pen_bar = exo(6);
% 1. Demographics
n_B = ((nu + omega) + sqrt((nu - omega)^2 + 4 * phi_B * (1 - nu))) / 2 - 1;
psi_y_B = phi_B / (1 + n_B - nu);
psi_o_B = (1 - omega) / (1 + n_B - gamma_B);
tilde_s_m_B = 1 / (1 + psi_y_B + psi_o_B);
tilde_s_y_B = psi_y_B * tilde_s_m_B;
tilde_s_o_B = psi_o_B * tilde_s_m_B;
n_y_B = n_B;
n_m_B = n_B;
n_o_B = n_B;
% 2. Normalization
tilde_R_B = R_star;
g_B = (1 + n_B) * (1 + x);
q_B = tilde_R_B / g_B;
% 3. Firm
kappa_B = (1 - alpha) * g_B / (tilde_R_B - (1 - delta));
y_B = (tilde_s_m_B^alpha * (kappa_B / g_B)^(1 - alpha))^(1 / alpha);
k_B = kappa_B * y_B;
w_B = alpha * y_B;
inv_B = (1 - (1 - delta) / g_B) * k_B;
% 4. Government
b_B = b_y * y_B;
gov_B = g_y * y_B;
e_B = (repl * w_B / tilde_s_m_B + pen_bar) * tilde_s_o_B;
t_B = e_B + gov_B - b_B * (1 - q_B);
% 5. Household
%% MPCs
mu_o_B = 1 - gamma_B * beta^sigma * tilde_R_B^(sigma - 1);
Om_m = @(m) omega + (1 - omega) * (mu_o_B / m)^(1 / (1 - sigma));
res_m = @(m) m - 1 + beta^sigma * (Om_m(m) * tilde_R_B)^(sigma - 1);
mu_m_B    = fzero(res_m, [1e-8, 1]);
Omega_m_B = Om_m(mu_m_B);
Om_y = @(m) nu + (1 - nu) * (mu_m_B / m)^(1 / (1 - sigma));
res_y = @(m) m - 1 + beta^sigma * (Om_y(m) * tilde_R_B)^(sigma - 1);
mu_y_B    = fzero(res_y, [1e-8, 1]);
Omega_y_B = Om_y(mu_y_B);
%% PDVs
s_o_B = tilde_R_B * e_B / (tilde_R_B - gamma_B * (1 + x));
s_m_B = (Omega_m_B - omega) / (Omega_m_B * tilde_R_B / (1 + x) - omega) * 1 / psi_o_B * s_o_B;
s_y_B = (Omega_y_B - nu) / (Omega_y_B * tilde_R_B / (1 + x) - nu) * psi_y_B * s_m_B;
h_m_B = Omega_m_B * tilde_R_B / (1 + x) / (Omega_m_B * tilde_R_B / (1 + x) - omega) * (w_B - t_B);
h_y_B = (Omega_y_B - nu) / (Omega_y_B * tilde_R_B / (1 + x) - nu) * psi_y_B * h_m_B;
% Assets
a_y_B = - mu_y_B * (h_y_B + s_y_B) / (1 - q_B * nu * (1 - mu_y_B));
a_m_B = (q_B * (1 - mu_m_B) * (1 - nu) * a_y_B + w_B - t_B - mu_m_B * (h_m_B + s_m_B)) / (1 - q_B * omega * (1 - mu_m_B));
a_o_B = (q_B * (1 - mu_o_B) * (1 - omega) * a_m_B + e_B - mu_o_B * s_o_B) / (1 - q_B * (1 - mu_o_B));
a_B   = a_y_B + a_m_B + a_o_B;
eta_y_B = a_y_B / a_B;
eta_m_B = a_m_B / a_B;
eta_o_B = a_o_B / a_B;
% Wealth
z_o_B = q_B * (eta_o_B + (1 - omega) * eta_m_B) * a_B + s_o_B;
z_m_B = q_B * (omega * eta_m_B + (1 - nu) * eta_y_B) * a_B + h_m_B + s_m_B;
z_y_B = q_B * nu * eta_y_B * a_B + h_y_B + s_y_B;
% Consumption
c_o_B = mu_o_B * z_o_B;
c_m_B = mu_m_B * z_m_B;
c_y_B = mu_y_B * z_y_B;
% 6. Balance of Payments
f_B = a_B - k_B - b_B;
nx_B = f_B * (1 - q_B);
fy_target_B = f_B / y_B;
% 7. Walras check
walras_B = y_B - c_o_B - c_m_B - c_y_B - inv_B - gov_B - nx_B;
for i=1:M_.endo_nbr
if strncmp(M_.endo_names{i},'AUX_EXO_LEAD',12)
        eval([M_.endo_names{i} ' = gamma_B;']);
end
end
params = M_.params;
for i=1:M_.endo_nbr, eval(['ys(' int2str(i) ') = ' M_.endo_names{i} ';']);
end
