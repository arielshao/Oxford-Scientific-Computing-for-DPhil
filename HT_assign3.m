%% Scientific Computing HT Assignment 3
% Name: Aili Shao    

%% Q1 Advection-Diffusion Equation in 1D
% explicit solution of the advection-diffusion equqation  u_t = u_xx+10u_x, u(-1)=u(1)=0.
% u(x,0)=exp(-10x^4/(1-x^2))

% Grid and initial data:
  h = .005;                    % small space step to achieve at least 4 digits accuaracy
  k = .4*h^2;                  % time step
  x = (-1+h:h:1-h)';                 
  u = exp(-10*x.^4./(1-x.^2)); % initial data

% Set-up for plot:
  hold off, shg
  plt = plot(x,u,'linewidth',4);
  axis([-1 1 -.01 1.01]), grid
  hold on;
   plt2 = plot(2, 0.5, 'r*', 'markersize', 20);
% Sparse matrix to execute finite difference operation:
  L = length(x);
  a = (1-2*k/h^2);
  b = k/h^2+10*k/(2*h); % add the 10u_x term
  c = k/h^2+-10*k/(2*h);
  main = a*sparse(ones(L,1));
  off1  = b*sparse(ones(L-1,1));
  off2  = c*sparse(ones(L-1,1));
  A = diag(main) + diag(off1,1) + diag(off2,-1);
 
Tf = 0.125; %final time
numsteps = ceil(Tf/k);
Tf = k*numsteps;
disp('type <return> to see solution')
pause
% Time-stepping:
for n=1:numsteps
  t = n*k;
  unew =A*u;
  u = unew;
  set(plt,'ydata',u)
  drawnow
  %pause
end
format long;
uf=u((-0.75+1)/h)
plot(-0.75, uf,'r*','markersize',20)

%% Q2 Leap frog for the heat equation 
%  See the paper at the back 
%  (b)
    h = 0.05;
    x = (-1+h:h:1-h)';                
    k =0.001;
    uold = exp(-10*x.^4./(1-x.^2));
    L = length(x);
    sigma=k/h^2;
    Aoff=sigma/2*sparse(ones(L-1,1));
    Adiag=(1-sigma)*sparse(ones(L,1));
    Acn= diag(Adiag) + diag(Aoff,1) + diag(Aoff,-1);
    Boff=-sigma/2*sparse(ones(L-1,1));
    Bdiag=(1+sigma)*sparse(ones(L,1));
    Bcn=diag(Bdiag) + diag(Boff,1)+ diag(Boff, -1);
    u=Bcn\(Acn*uold);
   
    a = (-2*2*k/h^2); b =2* k/h^2;
    main = a*sparse(ones(L,1));
    off  = b*sparse(ones(L-1,1));
    A = diag(main) + diag(off,1) + diag(off,-1);
  
    nsteps = 40;
    G20=[];
    for step = 2:nsteps
      unew = A*u + uold;                     % leap frog
      uold = u; u = unew;
      if (mod(step, 20) == 0)
          t=step*k;
          figure(step); axis([-1 1 -2 2]); plot(x,unew,'linewidth',2); grid on;
          title(['n = ' num2str(step) ', t = ' num2str(t)])
          G20=[G20; max(unew)]
          pause;
      end
      drawnow
    end
    G=(G20(end)/G20(end-1))^(1/20)
    g=4*k/h^2+sqrt(16*k^2/h^4+1)
    diff=g-G
    

%% Q3 The Gray-Scott equations
%(a)
figure(1);gs(32)
figure(2);gs(64)
figure(3);gs(128)
% figure(4);gs(256)

%(b) This code solves the Gray-Scott equations
 
%         u_t = 2e-5*laplacian(u) + F*(1-u) - u*v^2,
%         v_t = 1e-5*laplacian(v) - (F+K)*v + u*v^2,
%  
%      on [0 1]^2 from t=0 to t=5000, with initial condition
%  
%         u0(x,y) = 1 - exp(-100*((x-1/2.05)^2 + (y-1/2.05)^2)),
%         v0(x,y) = exp(-100*((x-1/2)^2 + 2*(y-1/2)^2)),
%  
%      with F=0.030 and K=0.057.
% It aims to compute u(0,0) with N grid points in each direction
% and time step dt=4,2,1. 
% Using the 'zebra' plot, we get negative values in black and positive values in white.

% For each N, the images show the zebra plots of u-0.5 for N grid points in each direction for
% dt =4, 2, 1 (from left to right). 
% For each N, the output values show the value u(0,0) computed from N grid
% points in each direction with time step dt=4, 2, 1 (from left to right).

% (c)   N\ dt     4            2            1
%       ------------------------------------------
%      32      0.785681     0.785690    0.785690
%      64      0.517885     0.518282    0.518279
%      128     0.516123     0.516523    0.516556
%      256     0.516120     0.516521    0.516554
%  From the table above, we obseve that the rate of convergence is
%  O((1/N)^2,dt^2) 
% u00 should be 0.5165( or say 0.517). N=128 and dt=2.
% (d) 
figure(4);gsspots(32)
figure(5);gsspots(64)
figure(6);gsspots(128)
%         u = spin2('GSspots');
%  
%      solves the Gray-Scott equations
%  
%         u_t = 2e-5*laplacian(u) + F*(1-u) - u*v^2,
%         v_t = 1e-5*laplacian(v) - (F+K)*v + u*v^2,
%  
%      on [0 1]^2 from t=0 to t=5000, with initial condition
%  
%         u0(x,y) = 1 - exp(-100*((x-1/2.05)^2 + (y-1/2.05)^2)),
%         v0(x,y) = exp(-100*((x-1/2)^2 + 2*(y-1/2)^2)), 
%  
%      with F=0.026 and K=0.059.
% 
%(e)
figure(7);gss(32)
figure(8);gss(64)
figure(9);gss(128)
%figure(10);gss(256)

% The images now show a mixture of spots and strips
% More grid points and smaller time step are needed to achieve high order
% accuaracy 

function u00=gss(N);
tic;
dom = [0 1 0 1]; x = chebfun('x',dom(1:2)); tspan = [0 5000];
S = spinop2(dom,tspan); F=0.027; K=0.0585;

S.lin = @(u,v) [2e-5*lap(u); 1e-5*lap(v)];
S.nonlin = @(u,v) [F*(1-u)-u.*v.^2;-(F+K)*v+u.*v.^2];
%S.init = chebfun2v(@(x,y)1 - exp(-100*((x-1/2.05)^2 + (y-1/2.05)^2)), ...
                  % @(x,y) exp(-100*((x-1/2)^2 + 2*(y-1/2)^2)), dom);
T=spinop2('gs');
S.init=T.init;
tic,

for i=1:3
    dt=2^(3-i);
    u = spin2(S,N,dt,'plot','off');
  subplot(1,3,i)
    plot(u{1}-.5,'zebra'), axis square off
    u00=u{1}(0,0);
    s=sprintf('u(0,0)= %8.6f\n', u00);
    title(s,'fontsize',8), drawnow
end
toc
end



