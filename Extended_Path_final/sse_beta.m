load('beta_grid_rstar.mat')
HLW_year  = readmatrix('hlw_star.xlsx','Sheet','for_grid','Range','A1:A55');
HLW_rate  = readmatrix('hlw_star.xlsx','Sheet','for_grid','Range','B1:B55');
HLW_gross = 1 + HLW_rate(:)/100;      % IF rates are in net % (2.5)
HLW_year  = HLW_year(:);
[yrs, im, ih] = intersect(Year, HLW_year);
Rm = R_store(im, :);
Ht = HLW_gross(ih);
fprintf('common years: %d-%d (%d pts)\n', yrs(1), yrs(end), numel(yrs));

sse = NaN(1, numel(beta_grid));
for k = 1:numel(beta_grid)
    if ok(k), sse(k) = sum((Rm(:,k) - Ht).^2); end
end
[sm, best] = min(sse);
fprintf('>>> estimated beta = %.4f (SSE = %.3e)\n', beta_grid(best), sm);

figure; plot(beta_grid, sse, 'o-','LineWidth',1.4); xlabel('\beta'); ylabel('SSE vs HLW'); grid on;
figure; plot(yrs, Rm(:,best),'-o', yrs, Ht,'-s','LineWidth',1.3);
legend('model (\beta*)','HLW'); xlabel('Year'); ylabel('gross r*'); grid on;
title(sprintf('Fit: \\beta = %.4f', beta_grid(best)));

% ===== STANDARD ERROR of beta from SSE curvature =====
% Fit a parabola to the SSE near the minimum (+-6-point window),
% then Var(beta) = 2*sigma2 / H, with H = d2(SSE)/db2 = 2*p(1).
lo = max(1, best-6);  hi = min(numel(beta_grid), best+6);
bw = beta_grid(lo:hi)';
sw = sse(lo:hi)';
p  = polyfit(bw, sw, 2);       % p(1)*b^2 + p(2)*b + p(3)
H        = 2*p(1);
beta_hat = -p(2)/(2*p(1));     % parabola vertex (smoothed minimum)
fprintf('beta_hat (vertex) = %.5f  (grid best = %.5f)\n', beta_hat, beta_grid(best));

% Newey-West (HAC) residual variance -> autocorrelation-robust SE
res    = Rm(:,best) - Ht;
n      = numel(res);
nparam = 1;
L   = round(4*(n/100)^(2/9));  % Newey-West truncation lag
lrv = sum(res.^2)/n;
for j = 1:L
    wj  = 1 - j/(L+1);                        % Bartlett weight
    gj  = sum(res(1+j:end).*res(1:end-j))/n;  % lag-j autocovariance
    lrv = lrv + 2*wj*gj;
end
sigma2_hac = lrv * n/(n-nparam);
SE = sqrt(2*sigma2_hac / H);
fprintf('SE (Newey-West) = %.5f  ->  95%% CI = [%.5f, %.5f]\n', ...
        SE, beta_hat-1.96*SE, beta_hat+1.96*SE);

% ===== sensitivity-table betas: center +/-1,2 SE =====
band_labels = {'-2 SE','-1 SE','center','+1 SE','+2 SE'};
band_betas  = beta_hat + [-2 -1 0 1 2]*SE;
in_grid = band_betas >= min(beta_grid) & band_betas <= max(beta_grid);
fprintf('\n--- sensitivity-table betas (Newey-West, sigma = %.2f) ---\n', sigma_fix);
for q = 1:numel(band_betas)
    flag = ''; if ~in_grid(q), flag = ' [outside grid - not simulated]'; end
    fprintf('  %-7s beta = %.5f%s\n', band_labels{q}, band_betas(q), flag);
end


% ===== two panels side by side, for the appendix =====
figApp = figure('Position',[100 100 1400 500],'Color','w');

% -- left: SSE vs beta --
ax1 = subplot(1,2,1);
plot(beta_grid, sse, 'o-','LineWidth',1.5,'Color',[0 0.447 0.741], ...
     'MarkerFaceColor',[0 0.447 0.741],'MarkerSize',4); hold on;
xline(beta_hat,'--','Color',[0.3 0.3 0.3],'LineWidth',1.0);
xlabel('\beta'); ylabel('SSE vs HLW'); grid on; box on;
title(sprintf('(a) Objective: SSE vs \\beta   (\\beta^* = %.4f)', beta_hat));

% -- right: model vs HLW fit --
ax2 = subplot(1,2,2);
plot(yrs, Rm(:,best),'-o','LineWidth',1.5,'Color',[0 0.447 0.741], ...
     'MarkerFaceColor',[0 0.447 0.741],'MarkerSize',4); hold on;
plot(yrs, Ht,'-s','LineWidth',1.5,'Color',[0.85 0.325 0.098], ...
     'MarkerFaceColor',[0.85 0.325 0.098],'MarkerSize',4);
lgd = legend('model (\beta^*)','HLW','Location','northeast');
xlabel('Year'); ylabel('gross r*'); grid on; box on;
title(sprintf('(b) Fit: model vs HLW   (\\beta^* = %.4f)', beta_hat));

% -- white theme + readable fonts on BOTH panels --
xlim(ax1, [1.0200 1.0300]);
xlim(ax2, [1970 2024]);
set([ax1 ax2], 'Color','w', 'XColor','k', 'YColor','k', ...
    'GridColor',[0.8 0.8 0.8], 'FontName','Times New Roman', ...
    'GridAlpha',1, 'FontSize',11, 'LineWidth',0.8);
set(findall(figApp,'Type','text'),'Color','k');
set(lgd, 'Color','w', 'TextColor','k', 'EdgeColor',[0.7 0.7 0.7]);


% -- export --
exportgraphics(figApp, 'output_sse_beta/appendix_beta_fit.pdf', 'ContentType','vector');
exportgraphics(figApp, 'output_sse_beta/appendix_beta_fit.png', 'Resolution',300);