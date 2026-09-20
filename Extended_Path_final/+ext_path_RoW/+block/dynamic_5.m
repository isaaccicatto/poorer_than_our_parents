function [y, T] = dynamic_5(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(47)=(1-params(9)+y(2)*x(1))/(1+y(50));
  y(51)=(1-params(9))/y(2)+x(1)-1;
  y(49)=params(10)+x(2)/y(3)-1;
end
