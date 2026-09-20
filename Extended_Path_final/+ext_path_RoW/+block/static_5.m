function [y, T, residual, g1] = static_5(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(1, 1);
  residual(1)=(y(2))-((1-params(9)+y(2)*x(1))/(1+y(5)));
if nargout > 3
    g1_v = NaN(1, 1);
g1_v(1)=1-x(1)/(1+y(5));
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 1, 1);
end
end
