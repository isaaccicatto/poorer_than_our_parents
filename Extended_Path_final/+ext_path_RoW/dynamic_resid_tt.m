function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 17
    T = [T; NaN(17 - size(T, 1), 1)];
end
T(1) = params(3)^params(4);
T(2) = y(123)^(params(4)-1);
T(3) = T(1)*(y(123)*y(129))^(params(4)-1);
T(4) = 1+T(3)*1/y(98);
T(5) = T(1)*(y(123)*y(130))^(params(4)-1);
T(6) = 1+T(5)*1/y(99);
T(7) = 1/(1-params(4));
T(8) = y(124)/(y(123)*y(129));
T(9) = y(100)*(y(129)-params(9))/(y(47)*(1+y(96)))+params(9)/(1+y(95))*y(101);
T(10) = y(124)/(y(123)*y(130));
T(11) = y(48)*(y(130)-params(10))/(1+y(95));
T(12) = y(101)*T(11)+params(10)/(1+y(94))*y(102);
T(13) = T(11)*y(103)+params(10)/(1+y(94))*y(104);
T(14) = y(81)^params(1);
T(15) = (y(26)/y(79))^(1-params(1));
T(16) = T(2)*T(1)*y(135);
T(17) = 1+1/y(97)*T(16);
end
