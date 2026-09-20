function [y, T, residual, g1] = dynamic_4(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(2, 1);
  residual(1)=(y(50))-(params(9)+y(3)*(1-params(10))-1);
  residual(2)=(y(48))-((x(2)+params(10)*y(3))/(1+y(50)));
if nargout > 3
    g1_v = NaN(3, 1);
g1_v(1)=1;
g1_v(2)=(-((-(x(2)+params(10)*y(3)))/((1+y(50))*(1+y(50)))));
g1_v(3)=1;
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 2, 2);
end
end
