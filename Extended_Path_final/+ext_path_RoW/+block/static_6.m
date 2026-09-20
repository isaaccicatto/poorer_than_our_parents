function [y, T] = static_6(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(6)=(1-params(9))/y(2)+x(1)-1;
  y(4)=params(10)+x(2)/y(3)-1;
end
