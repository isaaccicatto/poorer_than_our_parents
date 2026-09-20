function [ys,params,check] = ext_path_RoW_grid_steadystate(ys,exo,M_,options_)

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
gamma_RoW = exo(1); phi_RoW = exo(2);


% 1. Demographics
n_RoW = ((nu + omega) + sqrt((nu - omega)^2 + 4 * phi_RoW * (1 - nu))) / 2 - 1;
psi_y_RoW = phi_RoW / (1 + n_RoW - nu);
psi_o_RoW = (1 - omega) / (1 + n_RoW - gamma_RoW);
tilde_s_m_RoW = 1 / (1 + psi_y_RoW + psi_o_RoW);
tilde_s_y_RoW = psi_y_RoW * tilde_s_m_RoW; 
tilde_s_o_RoW = psi_o_RoW * tilde_s_m_RoW;
n_y_RoW = n_RoW; 
n_m_RoW = n_RoW; 
n_o_RoW = n_RoW;


% 2. Interest Rate
options_fzero = optimset('Display','off');
% initial guess
[tilde_R_RoW, ~, exitflag] = fzero(@get_excess_assets, 1.045, options_fzero);

if exitflag < 0
    check = 1; % SS fails
    return;
end


% 3. Normalization
g_RoW = (1 + n_RoW) * (1 + x);
q_RoW = tilde_R_RoW / g_RoW;


% 4. Firm
kappa_RoW = (1 - alpha) * g_RoW / (tilde_R_RoW - (1 - delta));
y_RoW = (tilde_s_m_RoW^alpha * (kappa_RoW / g_RoW)^(1 - alpha))^(1 / alpha);
k_RoW = kappa_RoW * y_RoW;
w_RoW = alpha * y_RoW;
inv_RoW = (1 - (1 - delta) / g_RoW) * k_RoW;


% 5. Government
b_RoW = b_y * y_RoW;
gov_RoW = g_y * y_RoW;
e_RoW = repl * w_RoW / tilde_s_m_RoW * tilde_s_o_RoW;
t_RoW = e_RoW + gov_RoW - b_RoW * (1 - q_RoW);


% 6. Household
% MPCs
mu_o_RoW = 1 - gamma_RoW * beta^sigma * tilde_R_RoW^(sigma - 1);

Om_m = @(m) omega + (1 - omega) * (mu_o_RoW / m)^(1 / (1 - sigma));
res_m = @(m) m - 1 + beta^sigma * (Om_m(m) * tilde_R_RoW)^(sigma - 1);
mu_m_RoW    = fzero(res_m, [1e-8, 1]);
Omega_m_RoW = Om_m(mu_m_RoW);

Om_y = @(m) nu + (1 - nu) * (mu_m_RoW / m)^(1 / (1 - sigma));
res_y = @(m) m - 1 + beta^sigma * (Om_y(m) * tilde_R_RoW)^(sigma - 1);
mu_y_RoW    = fzero(res_y, [1e-8, 1]);
Omega_y_RoW = Om_y(mu_y_RoW);

% PDVs
s_o_RoW = tilde_R_RoW * e_RoW / (tilde_R_RoW - gamma_RoW * (1 + x));
s_m_RoW = (Omega_m_RoW - omega) / (Omega_m_RoW * tilde_R_RoW / (1 + x) - omega) * 1 / psi_o_RoW * s_o_RoW;
s_y_RoW = (Omega_y_RoW - nu) / (Omega_y_RoW * tilde_R_RoW / (1 + x) - nu) * psi_y_RoW * s_m_RoW;
h_m_RoW = Omega_m_RoW * tilde_R_RoW / (1 + x) / (Omega_m_RoW * tilde_R_RoW / (1 + x) - omega) * (w_RoW - t_RoW);
h_y_RoW = (Omega_y_RoW - nu) / (Omega_y_RoW * tilde_R_RoW / (1 + x) - nu) * psi_y_RoW * h_m_RoW;

% Assets
a_y_RoW = - mu_y_RoW * (h_y_RoW + s_y_RoW) / (1 - q_RoW * nu * (1 - mu_y_RoW));
a_m_RoW = (q_RoW * (1 - mu_m_RoW) * (1 - nu) * a_y_RoW + w_RoW - t_RoW - mu_m_RoW * (h_m_RoW + s_m_RoW)) / (1 - q_RoW * omega * (1 - mu_m_RoW));
a_o_RoW = (q_RoW * (1 - mu_o_RoW) * (1 - omega) * a_m_RoW + e_RoW - mu_o_RoW * s_o_RoW) / (1 - q_RoW * (1 - mu_o_RoW));
a_RoW   = a_y_RoW + a_m_RoW + a_o_RoW;
eta_y_RoW = a_y_RoW / a_RoW;  
eta_m_RoW = a_m_RoW / a_RoW;  
eta_o_RoW = a_o_RoW / a_RoW;

% Wealth
z_o_RoW = q_RoW * (eta_o_RoW + (1 - omega) * eta_m_RoW) * a_RoW + s_o_RoW;
z_m_RoW = q_RoW * (omega * eta_m_RoW + (1 - nu) * eta_y_RoW) * a_RoW + h_m_RoW + s_m_RoW;
z_y_RoW = q_RoW * nu * eta_y_RoW * a_RoW + h_y_RoW + s_y_RoW;

% Consumption
c_o_RoW = mu_o_RoW * z_o_RoW;
c_m_RoW = mu_m_RoW * z_m_RoW;
c_y_RoW = mu_y_RoW * z_y_RoW;


% 6. Walras check
walras_RoW = y_RoW - c_o_RoW - c_m_RoW - c_y_RoW - inv_RoW - gov_RoW;


params = M_.params;
ys = zeros(M_.endo_nbr, 1); % guarantees the existence and feasibility of vector ys

for i = 1:M_.endo_nbr
    varname = M_.endo_names{i};
    
    if strncmp(varname, 'AUX_EXO_LEAD', 12)
        ys(i) = gamma_RoW;
    else
        eval(['ys(i) = ' varname ';']); 
    end
end


% Nested function: excess demand in capital market

    function err = get_excess_assets(R_guess)
        R_temp = R_guess;
        
        g_temp = (1 + n_RoW) * (1 + x);
        q_temp = R_temp / g_temp;

        kappa_temp = (1 - alpha) * g_temp / (R_temp - (1 - delta));
        y_temp = (tilde_s_m_RoW^alpha * (kappa_temp / g_temp)^(1 - alpha))^(1 / alpha);
        k_temp = kappa_temp * y_temp;
        w_temp = alpha * y_temp;

        b_temp = b_y * y_temp;
        gov_temp = g_y * y_temp;
        e_temp = repl * w_temp / tilde_s_m_RoW * tilde_s_o_RoW;
        t_temp = e_temp + gov_temp - b_temp * (1 - q_temp);

        mu_o_temp = 1 - gamma_RoW * beta^sigma * R_temp^(sigma - 1);

        Om_m = @(m) omega + (1 - omega) * (mu_o_temp / m)^(1 / (1 - sigma));
        res_m = @(m) m - 1 + beta^sigma * (Om_m(m) * R_temp)^(sigma - 1);
        mu_m_temp = fzero(res_m, [1e-8, 1], options_fzero);
        Omega_m_temp = Om_m(mu_m_temp);

        Om_y = @(m) nu + (1 - nu) * (mu_m_temp / m)^(1 / (1 - sigma));
        res_y = @(m) m - 1 + beta^sigma * (Om_y(m) * R_temp)^(sigma - 1);
        mu_y_temp = fzero(res_y, [1e-8, 1], options_fzero);
        Omega_y_temp = Om_y(mu_y_temp);

        s_o_temp = R_temp * e_temp / (R_temp - gamma_RoW * (1 + x));
        s_m_temp = (Omega_m_temp - omega) / (Omega_m_temp * R_temp / (1 + x) - omega) * 1 / psi_o_RoW * s_o_temp;
        s_y_temp = (Omega_y_temp - nu) / (Omega_y_temp * R_temp / (1 + x) - nu) * psi_y_RoW * s_m_temp;
        
        h_m_temp = Omega_m_temp * R_temp / (1 + x) / (Omega_m_temp * R_temp / (1 + x) - omega) * (w_temp - t_temp);
        h_y_temp = (Omega_y_temp - nu) / (Omega_y_temp * R_temp / (1 + x) - nu) * psi_y_RoW * h_m_temp;

        a_y_temp = - mu_y_temp * (h_y_temp + s_y_temp) / (1 - q_temp * nu * (1 - mu_y_temp));
        a_m_temp = (q_temp * (1 - mu_m_temp) * (1 - nu) * a_y_temp + w_temp - t_temp - mu_m_temp * (h_m_temp + s_m_temp)) / (1 - q_temp * omega * (1 - mu_m_temp));
        a_o_temp = (q_temp * (1 - mu_o_temp) * (1 - omega) * a_m_temp + e_temp - mu_o_temp * s_o_temp) / (1 - q_temp * (1 - mu_o_temp));
        a_temp   = a_y_temp + a_m_temp + a_o_temp;

        % excess supply (+) or excess demand (-) for assets
        err = a_temp - k_temp - b_temp;
    end

end