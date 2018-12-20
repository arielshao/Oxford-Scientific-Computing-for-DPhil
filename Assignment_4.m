%% Scientific Computing Assignment 4
% Name: Aili Shao    
%% Q1 A degenerate eigenvalue problem 
format long
%Build the matrices
D=diag(1:10);
e=ones(10,1);
T=spdiags([e  2*e  e], [-1 0 1], 10, 10);
S=zeros(10,10);
for i=1:10
    for j=1:10
        S(i,j)=sin(i*j);
    end
end
G=S'*S;

%Define the function A(g,t)
A=@(x) D+x(1)*G+x(2)*T;

%Define the minima function sep 
 sep=@(x)min(diff(sort(eig(A(x)))))
 g=[0:.1:5]; t= g;
 [gg,tt] = meshgrid(g,t);
 ff=zeros(size(gg));
 
% Plot the function sep
 for i=1:length(g)
    for j=1:length(t)
     ff(i,j)=sep([g(i),t(j)]);
    end
 end
 hold off, surf(gg,tt,ff)
 colorbar, grid on, hold on
 LW = 'linewidth'; MS = 'markersize';

% Search for the minima for different initial points x0

lambda=Inf; %assign lambda to be infinity first 
for i=1:length(gg)
    for j=1:length(tt)
        x0=[g(i),t(j)]; 
        [x,fval]=fminsearch(sep,x0);
        if x(1)>0 & x(2)>0   % only keep those positve x's (positve g and t)
           eigs=sort(eig(A(x)));  
           index=find(diff(eigs)<1e-4);
           lambdanew=eigs(index); % find the corresponding lambda)
           if lambdanew<lambda  
              lambda=lambdanew;  %update lambda if it is smaller
           end
        end
    end
end

lambda



%% Q2 A maximization problem 

% Plot of the function 
f=@(x,y)x.*exp(-(x.^2+y.^2)).*sin(5*(atan2(y,x)+sqrt(x.^2+y.^2)))
x = [-3:.05:3]; y = x;
[xx,yy] = meshgrid(x,y);
hold off;surf(xx,yy,f(xx,yy))          
colorbar, grid on, hold on

% Use fmincon to minimize the reformulated problem with constraint
% Note that the function here has 4 variables
fun=@(x)x(1)*exp(-x(1)^2-x(3)^2)*sin(5*(atan(x(3)/x(1))+...
      sqrt(x(1)^2+x(3)^2)))-...
        x(2)*exp(-x(2)^2-x(4)^2)*sin(5*(atan(x(4)/x(2))+...
      sqrt(x(2)^2+x(4)^2)))  
  
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];
% an initial guess satisfies the constraint
x0=[0.48,0,0.08,sqrt(1-0.48^2-0.08^2)];  % a good initial guess
nonlcon=@circlecon; % the nonlinear constraint
x=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon) %minimization 

% Find the max difference
max=abs(f(x(2),x(4))-f(x(1),x(3)))
    
function [c,ceq]=circlecon(x)
c=(x(1)-x(2))^2+(x(3)-x(4))^2-1;
ceq=[];
end



%% Model Answer

% eigenvalue problem
N=10; 
S=zeros(N);
for j= 1:N
    for i=1:N
        S(i,j)=sin(i*j);
    end
end
G=S'*S;
T=toeplitz([2 1 zeros(1,N-2)]);D=diag(1:N);
sep=@(g,t) min(diff(sort(eig(D+g*G+t*T))));

%Plot of seg(g,t)
np=120; g=linspace(0,3,np); t=linspace(0,4,np);
[gg,tt]=meshgrid(g,t);levels = 0:.025: .2; ss=[];
for i= 1:np
    for j=1:np
        ss(i,j)=sep(g(i),t(j));
    end
end
clf;contour(g,t,ss,levels)
colorbar
grid on, xlabel g, ylabel t

% The plots suggest that seg(g,t) has many zeros lying at discrete points
% In order to get the smallest eigenvalue, we try the following code:

% while 1
%     [g,t]=ginput(1); A=D+g*G+t*T;
%     e=sort(eig(A)); de=diff(e);
%     i=find(de==min(de)); e=e(i:i+1); 
%     disp([i e', diff(e)])
% end
% This enables us to click at a point on the plot and find the
% corresponding nearly-degenerate eigenvalue.

% Experiments suggest that the smallest lambda will lie near 
% g=0.95,t= 0.77


% Get the local minima
sep=@(x) min(diff(sort(eig(D+x(1)*G+x(2)*T))));
opts=optimset('TolFun',1e-15);
[x,fval]=fminsearch(sep,[.95 .78],opts)
A=D+x(1)*G+x(2)*T;
eigA=eig(A)



%% A maximisation problem 
% View this problem as an optimization problem in 3 variables: the two
% coordinates of one endpoint of the needle, and its angle.

x0=linspace(-1.5, 1.5,90)'; y0=x0;
[x,y]=meshgrid(x0,y0);
r =sqrt(x.^2+y.^2); t=atan2(y,x);
f=x.*exp(-r.^2).*sin(5*(t+r));
hold off, contour(x0,y0,f),colorbar
axis square, hold on
h=plot([0 1],[0 0],'.-r');

% First throw  10^5 random needles and pick the best one
dfbest=0;
for i=1:100000
    x1=3*rand-1.5; y1=3*rand-1.5;
    r1=sqrt(x1^2+y1^2);t1=atan2(y1,x1);
    f1=x1.*exp(-r1.^2).*sin(5*(t1+r1));
    tt=2*pi*rand;
    x2=x1+cos(tt); y2=y1+sin(tt);
    r2=sqrt(x2^2+y2^2);t2=atan2(y2,x2);
    f2=x2.*exp(-r2.^2).*sin(5*(t2+r2));
    df1=abs(f2-f1);
    if df1>dfbest
        set(h,'xdata',[x1, x2], 'ydata', [y1 y2]);
        dfbest=df1;
        x1best=x1; y1best=y1; ttbest=tt;
        disp(df1)
        pause(.1)
        title(num2str(df1))
    end
end
x1best, y1best, ttbest

% This gives us a good approximate solution x1=0.4884 y=-.0974,
% theta=3.2504. We then use this as initial guess to use fminunc




a4sol(x1best, y1best, ttbest)


