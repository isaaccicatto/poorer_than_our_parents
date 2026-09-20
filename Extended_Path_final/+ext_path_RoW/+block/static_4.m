function [y, T, residual, g1] = static_4(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(2, 1);
  residual(1)=(y(5))-(params(9)+y(3)*(1-params(10))-1);
  residual(2)=(y(3))-((x(2)+params(10)*y(3))/(1+y(5)));
if nargout > 3
    g1_v = NaN(4, 1);
g1_v(1)=1;
g1_v(2)=(-((-(x(2)+params(10)*y(3)))/((1+y(5))*(1+y(5)))));
g1_v(3)=(-(1-params(10)));
g1_v(4)=1-params(10)/(1+y(5));
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 2, 2);
end
end
