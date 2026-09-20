function [y, T, residual, g1] = static_2(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(4, 1);
  residual(1)=(y(1))-(x(2)*y(36)-(1-x(1))*y(38));
  residual(2)=(y(37))-((x(2)*y(36)+params(10)*y(37))/(1+y(1)));
  residual(3)=(y(36))-((params(9)*y(36)+(1-params(10))*y(37))/(1+y(1)));
  residual(4)=(y(38))-(((1-params(9))*y(36)+x(1)*y(38))/(1+y(1)));
if nargout > 3
    g1_v = NaN(12, 1);
g1_v(1)=1;
g1_v(2)=(-((-(x(2)*y(36)+params(10)*y(37)))/((1+y(1))*(1+y(1)))));
g1_v(3)=(-((-(params(9)*y(36)+(1-params(10))*y(37)))/((1+y(1))*(1+y(1)))));
g1_v(4)=(-((-((1-params(9))*y(36)+x(1)*y(38)))/((1+y(1))*(1+y(1)))));
g1_v(5)=1-params(10)/(1+y(1));
g1_v(6)=(-((1-params(10))/(1+y(1))));
g1_v(7)=(-x(2));
g1_v(8)=(-(x(2)/(1+y(1))));
g1_v(9)=1-params(9)/(1+y(1));
g1_v(10)=(-((1-params(9))/(1+y(1))));
g1_v(11)=1-x(1);
g1_v(12)=1-x(1)/(1+y(1));
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 4, 4);
end
end
