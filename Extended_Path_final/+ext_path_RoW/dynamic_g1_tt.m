function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = ext_path_RoW.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 22
    T = [T; NaN(22 - size(T, 1), 1)];
end
T(18) = getPowerDeriv(y(52)/y(53),T(7),1);
T(19) = getPowerDeriv(y(53)/y(54),T(7),1);
T(20) = getPowerDeriv(y(26)/y(79),1-params(1),1);
T(21) = getPowerDeriv(y(123)*y(129),params(4)-1,1);
T(22) = getPowerDeriv(y(123)*y(130),params(4)-1,1);
end
