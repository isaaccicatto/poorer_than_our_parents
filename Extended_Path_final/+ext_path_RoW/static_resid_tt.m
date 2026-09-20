function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 19
    T = [T; NaN(19 - size(T, 1), 1)];
end
T(1) = y(33)^(params(4)-1);
T(2) = params(3)^params(4);
T(3) = T(1)*T(2)*y(45);
T(4) = 1+1/y(7)*T(3);
T(5) = T(2)*(y(33)*y(39))^(params(4)-1);
T(6) = 1/y(8);
T(7) = 1+T(5)*T(6);
T(8) = T(2)*(y(33)*y(40))^(params(4)-1);
T(9) = 1/y(9);
T(10) = 1+T(8)*T(9);
T(11) = 1/(1-params(4));
T(12) = y(34)/(y(33)*y(39));
T(13) = y(10)*(y(39)-params(9))/(y(2)*(1+y(6)))+y(11)*params(9)/(1+y(5));
T(14) = y(34)/(y(33)*y(40));
T(15) = y(3)*(y(40)-params(10))/(1+y(5));
T(16) = y(11)*T(15)+y(12)*params(10)/(1+y(4));
T(17) = T(15)*y(13)+params(10)/(1+y(4))*y(14);
T(18) = y(36)^params(1);
T(19) = (y(26)/y(34))^(1-params(1));
end
