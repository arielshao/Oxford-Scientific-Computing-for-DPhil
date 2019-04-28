%% Scientific Computing HT Assignment 1
% Name: Aili Shao    
%% We formulate the ODE first: 
% Let the position of unit mass point p be zp=xp+iyp 
% and the positions of the other three fixed unit masses be
% $z_1=1, z_2=\cos(2\pi/3)+i\sin(2\pi/3)$, and
% $z_3=\cos(4\pi/3)+i\sin(4\pi/3)$
% then the ODE is $\frac{dz_p^2}{dt^2}=-\sum_{i=1}^3\frac{z_p-z_i}{|z_p-z_i|^3}$
% with initial conditions $z_p(0)=2-i2$ and $\frac{dz_p}{dt}(0)=0.$
% Note that we assume $G=1$ in the formulation

%% Solve the ODE using Chebfun 
cheboppref.setDefaults('ivpAbsTol',1e-14) %14-digits accurate
z1=1; %three fixed unit masses
z2=cos(2*pi/3)+i*sin(2*pi/3);
z3=cos(4*pi/3)+i*sin(4*pi/3);
N = chebop(0,40); % compute up to t=40
N.op=@(zp) diff(zp,2)+(zp-z1)/abs(zp-z1)^3 ...
+(zp-z2)/abs(zp-z2)^3+(zp-z3)/abs(zp-z3)^3; %define the differential operator
N.lbc=@(zp) [zp-(2-i*2); diff(zp)]; %set initial conditions
zp=N\0; %solve the ODE
figure(1);clf;plot(zp, 'interval' ,[0 40]);
xlabel x; ylabel y; title('Plot of orbit at t=40')
% The position of p at t=40
p=[real(zp(end)) imag(zp(end))];
hold on; scatter(p(1),p(2),'filled');
text(p(1),p(2)+0.2,num2str([p(1),p(2)]),'fontsize',7)
theta=atan2(imag(zp),real(zp));
theta=unwrap(theta);
figure(2);clf;;plot(theta)
xlabel t; ylabel theta; title('Plot of \theta (t)')

