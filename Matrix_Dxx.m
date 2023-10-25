function A = Matrix_Dxx(N);

e = ones(N-1,1);
A = spdiags([e -2*e e],[-1 0 1],N-1,N-1);

return