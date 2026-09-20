function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(4, 1);
  residual(1)=(y(46))-(x(2)*y(36)-(1-x(1))*y(38));
  residual(2)=(y(82))-((x(2)*y(36)+params(10)*y(37))/(1+y(46)));
  residual(3)=(y(81))-((params(9)*y(36)+(1-params(10))*y(37))/(1+y(46)));
  residual(4)=(y(83))-(((1-params(9))*y(36)+x(1)*y(38))/(1+y(46)));
if nargout > 3
    g1_v = NaN(7, 1);
g1_v(1)=1;
g1_v(2)=(-((-(x(2)*y(36)+params(10)*y(37)))/((1+y(46))*(1+y(46)))));
g1_v(3)=(-((-(params(9)*y(36)+(1-params(10))*y(37)))/((1+y(46))*(1+y(46)))));
g1_v(4)=(-((-((1-params(9))*y(36)+x(1)*y(38)))/((1+y(46))*(1+y(46)))));
g1_v(5)=1;
g1_v(6)=1;
g1_v(7)=1;
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 4, 4);
end
end
