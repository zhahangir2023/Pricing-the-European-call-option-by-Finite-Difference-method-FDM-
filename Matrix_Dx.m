function A = Matrix_Dx(N);

e = ones(N-1,1);
A = spdiags([-e e],[-1 1],N-1,N-1);

return