function A = Matrix_D(N);

e = ones(N-1,1);
A = spdiags([e],[0],N-1,N-1);

return