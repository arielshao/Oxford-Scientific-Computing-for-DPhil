%% Scientific Computing HT Assignment 2
% Name: Aili Shao    

%% Exercise 8.1 Exploring resonance to increase amplitude
% $y^{''}+y=1-\cos(\nu t)$, $y(0)=1, y'(0)=0.$ 
%(a)
nu=1;
L=chebop(0,100);L.lbc=[1;0];
L.op=@(t,y)diff(y,2)+y;
t=chebfun('t',[0 100]);
f=1-cos(nu*t);subplot(2,1,1);plot(f); xlabel t, ylabel f; title('Forcing function f=1-cos(t)');
y=L\f; subplot(2,1,2);plot(y); xlabel t, ylabel y; title('Response y(t), resonant frequency w=1');hold on;
t10=t(y==10);
tc=t10.ends(2)
%check
y(tc)

%(b) 
nu=1; y=chebfun('100',[0 100]);
while max(y)>10;
    nu=0.99*nu; % decrease nu by 1% each time
    L=chebop(0,100);L.lbc=[1;0];
    L.op=@(t,y)diff(y,2)+y;
    t=chebfun('t',[0 100]);
    f=1-cos(nu*t);y=L\f; %solve the equation with the new nnu value 
end
nu=nu/0.99

%% Exercise 11.1 Cleve Moler's favorite ODE.
% We can write this equation as two equations $y'=\pm\sqrt{1-y^2}=f(y)$ 
% Note that $\frac{\partial f}{\partial y}=\mp \frac{y}{\sqrt{1-y^2}}$ is
% bounded for $-1<y<1$ and thus $f(y)$ is Lipschitz continuous with respect
% to $y$. By Picard's theorem, we know each equation has a unique solution on
% for $-1<y<1.$  By inspection (or solve it explicitly), we know that $\pm
% \sin(t)$ is a unique solution for each equation on $t=(0, \pi/2)$.
% Thus, there are exactly two solutions on $t=[0,1]$, which are
% $\pm\sin(t).$ For $t=[0,2]$, $y(t)=\pm\sin(t)$ on $[0,\pi/2)$ and $y(t)=\pm\sin(c+t)$
% for any $c$ on $[\pi/2, 2]$ are solutions of this problem. Since $c$ is
% arbitrary, there are infinitely many solutions on this interval.

%% Exercise 13.3 Alternative choices of the Lorenz coefficient 28
N = chebop(0,100); N.lbc = [-15; -15; 20];
for r=20:2:24
N.op = @(t,u,v,w) [diff(u)-10*(v-u); ...
diff(v)-u*(r-w)+v; diff(w)-u*v+(8/3)*w];
[u,v,w] = N\0; 
subplot(3,2,(r-20)+1);plot(u,w);title(['Trajectory on u-w plane for r=',num2str(r)])
subplot(3,2,(r-20)+2);plot(u);title(['Plot of u(t) for r=',num2str(r)])
end
% The case r=24 seems to be chaotic since the u-w plots shows the 'butterfly' structure of the strange attractor.
% The case r=20 gives the clearest example of transient chaos since the plot of u(t) shows a smooth solution after t=60 (need to zoom in to see the plot).
%% Exercise 13.4 Lorenz equations with a breeze
N = chebop(0,100); N.lbc = [-15; -15; 20];
for a=20:5:30
N.op = @(t,u,v,w) [diff(u)-10*(v-u)-a; ...
diff(v)-u*(28-w)+v; diff(w)-u*v+(8/3)*w];
[u,v,w] = N\0; 
subplot(3,2,(a-20)/5*2+1);plot(u,w);title(['Trajectory on u-w plane for a=',num2str(a)])
subplot(3,2,(a-20)/5*2+2);plot(u);title(['Plot of u(t) for a=',num2str(a)])
end
% The case a=20 seems to be chaotic since the u-w plots shows the 'butterfly' structure of the strange attractor.
% The case a=30 gives the clearest example of transient chaos since the plot of u(t) shows a smooth solution after t=50 (need to zoom in to see the plot).
%% Exercise 15.5 A cyclic system of three ODEs
% $u'=u(1-u^2-bv^2-cw^2)$
% $v'=v(10v^2-bw^2-cu^2)$
% $w'=w(1-w^2-bu^2-cv^2)$
b=0.55; c=1.5;
% (a) 
N = chebop(0,800); %up tp t=800
N.op=@(t,u,v,w) [diff(u)-u*(1-u^2-b*v^2-c*w^2);
                 diff(v)-v*(1-v^2-b*w^2-c*u^2);
                 diff(w)-w*(1-w^2-b*u^2-c*v^2)];
N.lbc=[0.5;0.49;0.49]; %initial conditions
[u,v,w]=N\0;
subplot(2,2,1);plot(u);xlabel t, ylabel u;title('plot of u(t)');
subplot(2,2,2);plot(v);xlabel t, ylabel v;title('plot of v(t)');
subplot(2,2,3);plot(w);xlabel t, ylabel w;title('plot of w(t)');
subplot(2,2,4);plot3(u,v,w);xlabel u, ylabel v, zlabel w; title('plot of trajectory in u-v-w space');
% The plots of u(t), v(t) and w(t) oscillate bewtween 0 and 1 with
% increasing amplitude but decreasing frequency.
% The trajectory in u-v-w plane moves approximately in a cycle for large t.

% (b) See at the back of the paper

%% Exercise 16.1 Fisher equation
% $y^{''}+y-y^2=0$ for $x\in [-1,1]$ with $y(-1)=1, y(1)=0$.

N=chebop(-1,1); N.op=@(x,y) diff(y,2)+y-y^2;
N.lbc=1; N.rbc=0;
x=chebfun('x',[-1 1]);
% Constant initial guesses consistent with $y(0)\approx 0.6$ and $y(0)\approx -2.5$
N.init=0.6; y1=N\0; subplot(2,1,1);plot(y1);xlabel x, ylabel y, y1(0.5) 
N.init=-2.5; y2=N\0;subplot(2,1,2);plot(y2);xlabel x, ylabel y, y2(0.5)