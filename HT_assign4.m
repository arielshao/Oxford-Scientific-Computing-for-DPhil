%% Scientific Computing HT Assignment 4

% Name: Aili Shao    

%% Q1 Smallest eigenvalue of quartic Schrodinger problem
% -u_xx+x^4u=lambda u
format compact;
% Build the operator matrix using finite difference 
  h = 0.001;     % small gird size to ensure higher accuracy                 
  x = (-5:h:5)'; 
  m=length(x);
  e=ones(1,m)';
  Dxx=spdiags([e -2*e e],[-1 0 1],m,m);
  Dxx(1,2)=0; Dxx(m,m-1)=0; %zero boundary condition
  Dxx=Dxx/h^2;
  x4=x.^4; 
  M=-Dxx+diag(x4);
  
  % Sort out the smallest eigenvalue
  [V,D]=eigs(M,1,'SM');
  D
  
% Solve using Chebfun
dom = [-5 5];
% Assign the equation to two chebops N and B such that N(u) = lambda*B(u).
N = chebop(@(x,u) -diff(u,2)+x.^4.*u, dom);
B = chebop(@(x,u) u, dom);
% Assign boundary conditions to the chebop.
N.bc = @(x,u) [u(-5); u(5)];
options = cheboppref(); 
options.discretization = 'values';
% Number of eigenvalue and eigenmodes to compute.
k =6;
% Solve the eigenvalue problem.
[V, D] = eigs(N, B, k, options);
lambda=D(1)

% Plot of the eigenvalue
subplot(2,1,1);
plot(real(lambda), imag(lambda), '.', 'markersize', 25)
title('Eigenvalue'); xlabel('real'); ylabel('imag');
% Plot the eigenmodes.
subplot(2,1,2);
plot(real(V(:,1)), 'linewidth', 2);
title('Eigenmode'); xlabel('x'); ylabel('u');

%% Q2 Van der Pol equation with period T=10
% boundary condition is missing??? 

a=3; b=4;
if VanderPol(a)*VanderPol(b)>0;
    disp('Wrong choice')
else
    p=(a+b)/2;
    err=1;
    while err>1e-4 % to ensure higher order of accuracy
        if VanderPol(a)*VanderPol(p)<0
            b=p;
        else
            a=p;
        end
        p=(a+b)/2;
        err=abs(VanderPol(p));
    end
end
a=p
N = chebop(0,50); N.lbc = [1;0];
N.op = @(t,u) diff(u,2) - p*(1-u^2)*diff(u) + u;
u= N\0; [umin,t]=min(u,'local');
figure(2);clf;plot(u);hold on; scatter(t, umin, '*');ylim([-4,4]);
title('solution to the van der Pol equation')

type VanderPol.m

%% Q3 Allen-Cahn equation

% Note that when we use Chebfun to solve this equation it seems that the
% boundary conditions are being ignored. The solution converges to u=1 as
% t tends to infity rather than u=-1. In this case, we compute the time
% when the numeircal solution become NON-negative.
dom = [-1 1];
t = 0:0.001:3;
pdefun = @(t,x,u) .015.*diff(u,2)+u-u.^3;
% Assign boundary conditions.
bc.left = -1;
bc.right = -1;
x = chebfun(@(x) x, dom);
% initial condition.
u0 = 1-2.*x.^2;
opts = pdeset('Eps', 1e-4, 'Ylim', [-1,1.1]);
[t, u] = pde15s(pdefun, t, u0, bc, opts);
figure(3);waterfall(u, t)
xlabel('x'), ylabel('t')
tc=t(min(u)> -1e-4);
tc=tc(1)


%% Q3 Model Answer
% The critical time computed by bisection method by hand is around 3437.5
t=  0:34.375:3437.51;
pdefun= @(t,x,u) .015* diff(u,2)+u-u.^3;
bc.left =@(t,u) u+1 ; bc.right= @(t,u) u+1;
x=chebfun(@(x) x); u0= 1-2.*x.^2;
opts =pdeset('Eps', 1e-6, 'Ylim', [-1, 1.1]);
[t,u]= pde15s(pdefun, t, u0, bc, opts);
waterfall(u,t), view(110, 0), zlim([-1 1.2])

%% Q4 Inital boundary value problem
% We write the equation solver using Chebfun spin as the function
% 'chebsolve4.m', and do some test to see that for T=1, the solution
% does not blow up, but for T=1.1, the solution blows up. Now we use
% bisection method to decide t_c.
a=1; b=1.1;
chebsolve4(a);
try 
    chebsolve4(b)
catch err
end
disp(err.message)
while abs(a-b)>1e-10 %To ensure high accuracy
    midpoint=(a+b)/2
    try chebsolve4(midpoint);plot off;
        a=midpoint;
    catch
        b=midpoint;
    end
end
T=midpoint
u=chebsolve4(T)
figure(4);clf;plot(u); title(['Solution at t=' num2str(midpoint,"%1.12e")])
%To 3 digits of relative accuarcy, we can conclude that tc=1.00
type chebsolve4.m


%% Q4 Model Answer 
% Note that spin is for periodic problem, thus we should use pde15s again
% The critical time is about 6.0589 based on the following computation.
t=linspace(0, 6.058, 100);
pdefun=@(t,x,u) diff(u,2)+diff(u)+exp(u);
bc.left='dirichlet'; bc.right ='dirichlet';
x=chebfun(@(x) x); u0=chebfun(0);
opts=pdeset('Eps', 1e-7);
[t,u]=pde15s(pdefun, t, u0, bc, opts);
waterfall(u,t), xlabel('x'), ylabel('t')
max(u(:,end)), axis([-1 1 0 6.1 0 8])
%% Q5 Advection-diffusion equation 
% Use 1D finite difference to solve the problem
% Build 1D finite difference matrix
dx=0.05; dt=.25*dx^2; % to ensure at least 3 digits accuaracy 
N=8/dx+1; % number of grid points
x=(-4:dx:4)'; 
Ix = speye(N,N); Iy=speye(N,N);
e = ones(N,1);
D1xx = spdiags([e  -2*e  e], [-1 0 1], N, N);
D1xx(1,2) =0; % zero boundary condition
D1xx(N,N-1)=0;
D1xx = D1xx/dx^2;
D1xc = spdiags([-e  e], [-1 1], N, N);
D1xc(1,2)=2;D1xc(N,N-1)=-2;
D1xc = D1xc/(2*dx); % use centre differece for 1st derivative
M=D1xx-20*D1xc;
% initial condition 
u0=zeros(size(x));
y=1-abs(-1+dx:dx:1-dx);
u0((abs(x)<1))=y;
% By expriments, we show that m(t) shift to the right
Tf=.5; m=[];u=u0;
numsteps = ceil(Tf / dt);
%k = Tf / numsteps
Tf = dt*numsteps;
for n=1:numsteps
    t=dt*n;
    unew=u+dt*(M*u);
    u=unew;
    m=[m; x(u==max(u))];
    plot(x,u); %pause  % use pause to see the shift
end
 max(m)

%% Q6 Heat equation on a square
% We use finite difference to solve the problem
% Build 1D finite difference matrix
dx=0.0025; dt=.25*dx^2; %fine grid points to ensure high order accuracy
N=1/dx+1; % number of grid points
x1d=0:dx:1; y1d=0:dx:1;
[xx,yy]=meshgrid(x1d,y1d);
Ix = speye(N,N); Iy=speye(N,N);
e = ones(N,1);
D1xx = spdiags([e  -2*e  e], [-1 0 1], N, N);
D1xx(1,2) = 0; % zero boundary condition
D1xx(N,N-1) = 0;
D1xx = D1xx/dx^2; D1yy=D1xx; %grid size is the same on both x and y directions
Dxx = kron(Iy, D1xx);
Dyy = kron(D1yy, Ix);

u0=zeros(size(xx));
u0(xx>= 0.25 & xx<=0.75 & yy>=0.25 & yy<=0.75)=1; 
u=u0(:); t=0;
while max(u)>0.5
    unew=u+dt*(Dxx*u+Dyy*u);
    u=unew;
    t=t+dt;
end
t
max(u)
figure(61); clf; surf(xx, yy, u0);
title('initial condition');
xlabel('x'); ylabel('y');
zlabel('u(t,x,y)');
uplot=reshape(u,size(xx));
figure(62); clf; surf(xx, yy, uplot);
title(['Solution at t=' num2str(t)]);
xlabel('x'); ylabel('y');
zlabel('u(t,x,y)');

% Note that the grid is not fine enough. The critical time should be around
% t=0.0281

%% Model Answer (use CrankNicolson)
J=1000; h=1/J; s=(h:h:1)'; k=.0001;
[xx,yy]=meshgrid(s,s);
x=xx(:); y=yy(:);
u=double(abs(x-.5)<.25 & abs(y-.5)<.25);
I=speye(J); II=speye(J^2);
D=h^(-2)*toeplitz([-2 1 zeros(1, J-2)]);
L=kron(I,D)+kron(D, I);
A= II+k*L/2; B=II-k*L/2;
t=0;
while max(u)>.5
    t= t+k; u= B\(A*u);
end
t
