function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = ext_path_RoW.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 24
    T = [T; NaN(24 - size(T, 1), 1)];
end
T(20) = getPowerDeriv(y(7)/y(8),T(11),1);
T(21) = getPowerDeriv(y(8)/y(9),T(11),1);
T(22) = getPowerDeriv(y(26)/y(34),1-params(1),1);
T(23) = getPowerDeriv(y(33)*y(39),params(4)-1,1);
T(24) = getPowerDeriv(y(33)*y(40),params(4)-1,1);
end
