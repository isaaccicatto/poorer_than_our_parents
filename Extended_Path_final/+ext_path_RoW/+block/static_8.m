function [y, T] = static_8(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(21)=y(12)+y(14)+y(15)*y(43)*y(35)*params(10);
  y(20)=y(11)+y(13)+y(15)*y(35)*(params(9)*y(42)+(1-params(10))*y(43));
  y(19)=y(10)+y(15)*y(35)*(y(41)+(1-params(9))*y(42));
  y(28)=y(26)-y(26)*(1-params(2))/y(34);
  y(24)=y(9)*y(21);
  y(23)=y(8)*y(20);
  y(22)=y(7)*y(19);
  y(44)=y(25)-y(22)-y(23)-y(24)-y(28)-y(30);
end
