function [y, T] = dynamic_7(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(66)=y(57)+y(59)+y(15)*y(43)*y(80)*params(10);
  y(65)=y(56)+y(58)+y(15)*y(80)*(params(9)*y(42)+(1-params(10))*y(43));
  y(64)=y(55)+y(15)*y(80)*(y(41)+(1-params(9))*y(42));
  y(73)=y(71)-y(26)*(1-params(2))/y(79);
  y(69)=y(54)*y(66);
  y(68)=y(53)*y(65);
  y(67)=y(52)*y(64);
  y(89)=y(70)-y(67)-y(68)-y(69)-y(73)-y(75);
end
