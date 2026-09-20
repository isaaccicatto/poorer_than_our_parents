% ============================================================================
% RUN_BETA_GRID — runs the transition for a vector of betas (sigma fixed) and
% stores each run's r* series in a column, side by side.
% Output: beta_grid_rstar.csv (Year + one column per beta) + .mat + figure.
%
% SETUP (do this ONCE, see instructions in the chat):
%   1) copy ext_path_XX.mod           -> ext_path_RoW_grid.mod
%      copy ext_path_XX_steadystate.m -> ext_path_RoW_grid_steadystate.m
%   2) in ext_path_grid.mod, DELETE the lines "beta = ...;" and "sigma = ...;"
%      and replace them with:   @#include "beta_setting.mod"
%   3) in ext_path_grid.mod: set horizon_plot = 60; and comment out the
%      PLOTS and CSV EXPORT blocks (they save time and are not needed here)
% ============================================================================
clear; close all; clc;
% ----------------------------- CONFIG ---------------------------------------
FNAME     = 'ext_path_RoW_grid';     % name of the copied .mod (without extension)
RVAR      = 'tilde_R_RoW';           % series to extract (for SOE: 'tilde_R_A')
beta_grid = 1.0200:0.00025:1.0300;   % <- edit the grid here
sigma_fix = 0.5;                     % fixed sigma (collinear with beta; see 9 runs)
T_keep    = 56;                      % 1970..2025
% ----------------------------------------------------------------------------
nb      = numel(beta_grid);
R_store = NaN(T_keep, nb);
ok      = false(1, nb);
for ib = 1:nb
    b = beta_grid(ib);
    fprintf('\n===== RUN %d/%d : beta = %.4f (sigma = %.2f) =====\n', ...
            ib, nb, b, sigma_fix);
% 1) write the include that the .mod reads during preprocessing
    fid = fopen('beta_setting.mod', 'w');
    fprintf(fid, 'beta = %.6f;\nsigma = %.6f;\n', b, sigma_fix);
    fclose(fid);
% 2) run dynare — NOCLEARALL is essential (otherwise it wipes the wrapper)
try
        eval(['dynare ' FNAME ' noclearall nolog']);
% 3) extract the r* series (columns of realized_path: 1 = t0 = 1970)
        r_idx = strcmp(M_.endo_names, RVAR);
        assert(any(r_idx), 'variable %s not found in M_.endo_names', RVAR);
        R_store(:, ib) = realized_path(r_idx, 1:T_keep)';
        ok(ib) = true;
catch ME
        warning('beta = %.4f FAILED: %s  (moving on to the next)', b, ME.message);
end
end
% ----------------------------- SAVE -----------------------------------------
Year = (1970:1970 + T_keep - 1)';
cn   = strrep(compose('b_%.4f', beta_grid), '.', 'p');
Tout = [table(Year) array2table(R_store, 'VariableNames', cn)];
writetable(Tout, 'beta_grid_rstar.csv');
save('beta_grid_rstar.mat', 'beta_grid', 'sigma_fix', 'R_store', 'Year', 'ok');
fprintf('\nSaved: beta_grid_rstar.csv / .mat  (%d of %d runs ok)\n', sum(ok), nb);
% ----------------------------- PLOT -----------------------------------------
figure('Position', [100 100 900 500]); hold on;
cmap = turbo(nb);
for ib = 1:nb
if ok(ib)
        plot(Year, R_store(:, ib), 'LineWidth', 1.4, 'Color', cmap(ib,:), ...
'DisplayName', sprintf('\\beta = %.4f', beta_grid(ib)));
end
end
xlabel('Year'); ylabel('gross r*'); grid on;
legend('Location', 'northeast'); title('Model r* path across \beta (\sigma fixed)');