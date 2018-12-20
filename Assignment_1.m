%% Scientific Computing Assignment 1
% Name: Aili Shao    
%% Q1 Compute A^2, A^(-1), det(A) and eigenvalues of A
A=[1 2 3 4; 1 3 2 4; 1 4 2 3; 4 1 2 3]
% Square of matrix A is
square_A=A^2
% Inverse of matix A is
inverse_A=inv(A)
% Determinant of matrix A is
d=det(A)
% Eigenvalues of matrix A
evs=eigs(A)
% Product of the eigenvalues of matrix A
product=prod(evs)
% Compare the product of eigenvalues of A and the determinant of A
err=product-det(A)

%% Q2 Compute exp(A) using expm and Taylor series

% (a) Use Matlab's expm command
exp1=expm(A)

% (b) Use taylor series to compute exp(A)
N=100; % Number of terms for the Taylor polynomial
exp2=zeros(4,4);
for i=1:N;
    exp2=exp2+A^(i-1)/factorial(i-1);
end
exp2

% Compare the two results
exp_err=norm(exp1-exp2)

%% Q3 Sparse matrices and pcg code to solve Ax=b
D=@(N)sparse(toeplitz([2 -1 zeros(1,N-2)]));
I=@(N)speye(N);
A=@(N)kron(I(N),kron(I(N),D(N)))+kron(I(N),kron(D(N),I(N)))...
    +kron(D(N),kron(I(N),I(N)));

% The dimension of A(N) is N^3 * N^3.
full(A(2))
spy(A(4)) % spy plot of A
b=@(N)ones(length(A(N)),1); % create a column vector with all entries equal to 1
x=@(N)pcg(A(N),b(N),1e-12,200);
tic;
sol40=x(40);
pcgtime=toc
sol40(1:5)





