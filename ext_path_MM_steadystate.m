function [ys,params,check] = ext_path_MM_steadystate(ys,exo,M_,options_)

check = 0;
alpha  = M_.params(strcmp(M_.param_names,'alpha'));
delta  = M_.params(strcmp(M_.param_names,'delta'));
beta   = M_.params(strcmp(M_.param_names,'beta'));
sigma  = M_.params(strcmp(M_.param_names,'sigma'));
b_y    = M_.params(strcmp(M_.param_names,'b_y'));
g_y    = M_.params(strcmp(M_.param_names,'g_y'));
repl   = M_.params(strcmp(M_.param_names,'repl'));
x      = M_.params(strcmp(M_.param_names,'x'));
omega  = M_.params(strcmp(M_.param_names,'omega'));
nu     = M_.params(strcmp(M_.param_names,'nu'));
gamma_MM = exo(1); phi_MM = exo(2); R_star = exo(3);


% 1. Demographics
n_MM = ((nu + omega) + sqrt((nu - omega)^2 + 4 * phi_MM * (1 - nu))) / 2 - 1;
psi_y_MM = phi_MM / (1 + n_MM - nu);
psi_o_MM = (1 - omega) / (1 + n_MM - gamma_MM);
tilde_s_m_MM = 1 / (1 + psi_y_MM + psi_o_MM);
tilde_s_y_MM = psi_y_MM * tilde_s_m_MM; 
tilde_s_o_MM = psi_o_MM * tilde_s_m_MM;
n_y_MM = n_MM; 
n_m_MM = n_MM; 
n_o_MM = n_MM;


% 2. Normalization
tilde_R_MM = R_star;
g_MM = (1 + n_MM) * (1 + x);
q_MM = tilde_R_MM / g_MM;


% 3. Firm
kappa_MM = (1 - alpha) * g_MM / (tilde_R_MM - (1 - delta));
y_MM = (tilde_s_m_MM^alpha * (kappa_MM / g_MM)^(1 - alpha))^(1 / alpha);
k_MM = kappa_MM * y_MM;
w_MM = alpha * y_MM;
inv_MM = (1 - (1 - delta) / g_MM) * k_MM;


% 4. Government
b_MM = b_y * y_MM;
gov_MM = g_y * y_MM;
e_MM = repl * w_MM / tilde_s_m_MM * tilde_s_o_MM;
t_MM = e_MM + gov_MM - b_MM * (1 - q_MM);


% 5. Household
%% MPCs
mu_o_MM = 1 - gamma_MM * beta^sigma * tilde_R_MM^(sigma - 1);

Om_m = @(m) omega + (1 - omega) * (mu_o_MM / m)^(1 / (1 - sigma));
res_m = @(m) m - 1 + beta^sigma * (Om_m(m) * tilde_R_MM)^(sigma - 1);
mu_m_MM    = fzero(res_m, [1e-8, 1]);
Omega_m_MM = Om_m(mu_m_MM);

Om_y = @(m) nu + (1 - nu) * (mu_m_MM / m)^(1 / (1 - sigma));
res_y = @(m) m - 1 + beta^sigma * (Om_y(m) * tilde_R_MM)^(sigma - 1);
mu_y_MM    = fzero(res_y, [1e-8, 1]);
Omega_y_MM = Om_y(mu_y_MM);

%% PDVs
s_o_MM = tilde_R_MM * e_MM / (tilde_R_MM - gamma_MM * (1 + x));
s_m_MM = (Omega_m_MM - omega) / (Omega_m_MM * tilde_R_MM / (1 + x) - omega) * 1 / psi_o_MM * s_o_MM;
s_y_MM = (Omega_y_MM - nu) / (Omega_y_MM * tilde_R_MM / (1 + x) - nu) * psi_y_MM * s_m_MM;
h_m_MM = Omega_m_MM * tilde_R_MM / (1 + x) / (Omega_m_MM * tilde_R_MM / (1 + x) - omega) * (w_MM - t_MM);
h_y_MM = (Omega_y_MM - nu) / (Omega_y_MM * tilde_R_MM / (1 + x) - nu) * psi_y_MM * h_m_MM;

% Assets
a_y_MM = - mu_y_MM * (h_y_MM + s_y_MM) / (1 - q_MM * nu * (1 - mu_y_MM));
a_m_MM = (q_MM * (1 - mu_m_MM) * (1 - nu) * a_y_MM + w_MM - t_MM - mu_m_MM * (h_m_MM + s_m_MM)) / (1 - q_MM * omega * (1 - mu_m_MM));
a_o_MM = (q_MM * (1 - mu_o_MM) * (1 - omega) * a_m_MM + e_MM - mu_o_MM * s_o_MM) / (1 - q_MM * (1 - mu_o_MM));
a_MM   = a_y_MM + a_m_MM + a_o_MM;
eta_y_MM = a_y_MM / a_MM;  
eta_m_MM = a_m_MM / a_MM;  
eta_o_MM = a_o_MM / a_MM;

% Wealth
z_o_MM = q_MM * (eta_o_MM + (1 - omega) * eta_m_MM) * a_MM + s_o_MM;
z_m_MM = q_MM * (omega * eta_m_MM + (1 - nu) * eta_y_MM) * a_MM + h_m_MM + s_m_MM;
z_y_MM = q_MM * nu * eta_y_MM * a_MM + h_y_MM + s_y_MM;

% Consumption
c_o_MM = mu_o_MM * z_o_MM;
c_m_MM = mu_m_MM * z_m_MM;
c_y_MM = mu_y_MM * z_y_MM;


% 6. Balance of Payments
f_MM = a_MM - k_MM - b_MM;
nx_MM = f_MM * (1 - q_MM);
fy_target_MM = f_MM / y_MM;


% 7. Walras check
walras_MM = y_MM - c_o_MM - c_m_MM - c_y_MM - inv_MM - gov_MM - nx_MM;

for i=1:M_.endo_nbr
    if strncmp(M_.endo_names{i},'AUX_EXO_LEAD',12)
        eval([M_.endo_names{i} ' = gamma_MM;']);
    end
end

params = M_.params;
for i=1:M_.endo_nbr, eval(['ys(' int2str(i) ') = ' M_.endo_names{i} ';']); 
end
