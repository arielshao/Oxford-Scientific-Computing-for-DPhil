%% Scientific Computing Assignment 2
% Name: Aili Shao    

%% Q1 Condition number 
%Build the discretizing matrix A
D=@(n)sparse(toeplitz([2 -1 zeros(1,n-3)]));
I=@(n)speye(n-1);
A=@(n)kron(I(n),kron(I(n),D(n)))+kron(I(n),kron(D(n),I(n)))...
    +kron(D(n),kron(I(n),I(n)));

%Compute the condition number of A for n=3:12
k_A=[];
n=(3:12)';
for i=1:length(n)
    nn=n(i);
    evs=eig(full(A(nn)));
    k= max(evs)/min(evs);
    k_A=[k_A k];
end

% Print out the table 
Condition_Number_of_A=k_A';
table(n,Condition_Number_of_A)


% The results for n=6 and n=12 suggest that the condition number is
% appproximately 0.4n^2, and this is confirmed by a plot
% Display the results in a log-log plot
figure(1);clf;
loglog(n,Condition_Number_of_A,'-s');
xlabel('dimension', 'fontsize', 16)
ylabel('k(A)','fontsize', 16)
grid on;
title('loglog plot of condition number against dimension of matrix A')
hold on;
loglog(n,0.4*n.^2,'-r');
text(1e+1,2e+1,'0.4n^2','fontsize', 28)
% Alpha is around 2 while C is around 0.4
%% Q2 Solution by direct methods 
b=@(n)ones(length(A(n)),1);
time=[];n=(5:30)';
for i=1:length(n)
    nn=n(i);
    AA=A(nn);bb=b(nn); %Build the matrix and RHS first
    tic;
    x=AA\bb;
    time=[time  toc];
end
figure(2);clf;
loglog(n, time','-s');
xlabel('dimension')
ylabel('time(secs)')
grid on;
title('loglog plot of time against dimension of matrix A')
hold on;
loglog(n,2e-8*n.^5,'-r');
axis([5 50 1e-4 30]);
set(gca,'xtick',[5 10:10:50])
text(11,0.2,'2e-8n^5','fontsize', 28)

% The alpha is around 5 here. The dimension of this matrix is (n-1)^3, so
% the time would be proportional to n^9, which is hugely worse, if
% backslash did not take advantage of sparsity


%% Q3 Solution by CG 
time2=[];n=(5:40)';
for i=1:length(n)
    nn=n(i);
    AA=A(nn);
    bb=b(nn);
    tic
    x=pcg(AA,bb,1e-10,500);
    time2=[time2 toc];
end 
T=time2';
figure(3);clf;
loglog(n, T,'-s');
xlabel('dimension','fontsize', 16)
ylabel('time(secs)','fontsize',16)
grid on;
title('loglog plot of time against dimension of matrix A')
hold on;
loglog(n,1e-7*n.^4,'-r');
text(1e+1,1e-3,'1e-7n^4','fontsize', 28)
axis([5 50 1e-4 30]);
set(gca,'xtick',[5 10:10:50])

%We can explain the n^4 timing as follows.
%From problem 1, we know that the condition number is proportional to n^2,
%thus CG theory tells us that the number of steps for convergence will
%scale at worst linearly with, and this is confirmed by the numbers printed
%by pcg during execution. Each step, however, requires work proportional to
%n^3, since the biggest part of the work is a matrix-vector multiplication
%and the matrix has dimension O(n^3) with a constant number of nonzeros on
%each row.

% For n=50, following the dashed lines on the plots, we can estimate 6
% seconds for backslash and 0.6 seconds for unpreconditioned CG. For n=100,
% we get 200 seconds for backslash and 10 seconds for unpreconditioned CG.
%% Q4 Huge matrix with prime numbers on diagonal entries
format long
tic;
% List the first 100,000 prime numbers
p=primes(1299709)';
% Dimension of A
N=100000;
% Find the index matrix (the position for 1's)
% M=round(log2(100000));
M=floor(log2(100000));  % use floor instead of round
d1=2.^(0:1:M); % Note that the length of d1 is M+1
d2=-d1(end:-1:1); d=[d2 0 d1]; % create the index matrix used in spdiags
B=ones(N,2*(M+1)+1);
B(:,M+2)=p; %values corresponding to the diagonal entries of A
A=spdiags(B, d, N,N);% Create the matrix A using spdiags
b=(1:100000)';
% Use Preconditioned Conjugate Gradients Method and take the diagonal matrix 
% as preconditioner to solve the linear systems
% Note that A is symmetric and positive definite
tic;
x=pcg(A,b,1e-14,100, diag(sparse(p)));
pcg_time=toc
x(N/2)
total_time=toc 

            