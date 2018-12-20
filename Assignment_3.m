%% Scientific Computing Assignment 3
% Name: Aili Shao    
%% Q1 Fitting ellipse via least squares
%(a)  The matrix is  
%   |x_1^2 | x_1y_1 | y_1^2|
%   |x_2^2 | x_2y_2 | y_2^2|
%   |x_3^2 | x_3y_3 | y_3^2|
%   | ...  |  ...   |  ... |
%   |x_n^2 | x_ny_n | y_n^2|


%(b)
x=[3,1,0,-1,-2,0,-2,2]';
y=[3,-2,3,2,-2,-4,0,0]';
figure(1);clf;
[b,c,d]=ellipse(x,y) % The function is at the back

% (c) 
figure(2);clf
axis([-3 3 -3 3]), axis manual, hold on, grid on
x=[]; y=[]; button=1;
disp('input points with mouse, button >=2 for final point')
while button==1 
[xx,yy,button]=ginput(1);
x=[x; xx]; 
y=[y; yy];
plot(xx,yy,'x');
end
[b,c,d]=ellipse(x,y)


%% Q2 Solution of the Laplace equation on a square - polynomial best fitting
cheb=@(n) cos(pi*((1:n)-.5)/n);
alpha=1:7;
for k=2.^alpha
    n=2*k+1;
% Compute the chebshev points on each boundary
yleft=cheb(2*n); zleft=yleft'*1i-ones(2*n,1);
yright=cheb(2*n) ;zright=yright'*1i+ones(2*n,1);
xbottom=cheb(2*n); zbottom=xbottom'-ones(2*n,1)*1i;
xup=cheb(2*n);zup=xup'+ones(2*n,1)*1i;

% Build the vector z
z=[zleft; zbottom; zright; zup];
A=zeros(8*n, n);
% Construct the least square matrix
A(:,1)=ones(8*n,1);
for j=1:k
    A(:,2*j)=real(z.^j);
    A(:,2*j+1)=imag(z.^j);
end

%Build the RHS 

fleft=zeros(2*n,1); fright=zeros(2*n,1);
fup=1-xup'.^2; fbottom=zeros(2*n,1);

f=[fleft; fbottom; fright; fup];

% Solve the least square problem using backslash (which uses QR
% factorization)
c= A\f;

x=0.99; y=0.99;
% Compute u(0.99,0.99)
u=c(1);
newz=x+y*1i;
for j=1:k
    unew=u+real(newz^j)*c(2*j)+imag(newz^j)*c(2*j+1);
    u=unew;
end
fprintf('k= %3d  u=%10.8f\n' , k, u)
end

% The best estimate is  u=0.01923180 when k=64
%% Q2 Solution of the Laplace equation on a square -rational functions best fitting
cheb=@(n) cos(pi*((1:n)-.5)/n);

alpha=1:6;
for k=2.^alpha
    n=8*k+1;
% Compute the chebshev points on each boundary
yleft=cheb(2*n); zleft=yleft'*1i-ones(2*n,1);
yright=cheb(2*n) ;zright=yright'*1i+ones(2*n,1);
xbottom=cheb(2*n); zbottom=xbottom'-ones(2*n,1)*1i;
xup=cheb(2*n);zup=xup'+ones(2*n,1)*1i;

% Build the vector z
z=[zleft; zbottom; zright; zup];

% Create the rays from four corner points
z_northwest=zeros(k,1);
z_northeast=zeros(k,1);
z_southwest=zeros(k,1);
z_southeast=zeros(k,1);
d=zeros(k,1);
for j=1:k
    d(j)=2*exp(-sqrt(j-1));
    z_northwest(j)=-1-d(j)/sqrt(2) +(1+d(j)/sqrt(2))*1i;
    z_northeast(j)=1+d(j)/sqrt(2) +(1+d(j)/sqrt(2))*1i;
    z_southwest(j)=-1-d(j)/sqrt(2) +(-1-d(j)/sqrt(2))*1i;
    z_southeast(j)=1+d(j)/sqrt(2) +(-1-d(j)/sqrt(2))*1i;
end
z_corner=[z_northwest ;z_northeast;z_southeast;z_southwest];
% Check the position of these points -> consistent with the fig in the
% question
figure(3);clf;plot(z_corner,'r.')

A=zeros(8*n, n);
% Construct the least square matrix
A(:,1)=ones(8*n,1);
for j=1:k
    A(:,8*(j-1)+2)=real(d(j)./(z-z_northwest(j)));
    A(:,8*(j-1)+3)=imag(d(j)./(z-z_northwest(j)));
    A(:,8*(j-1)+4)=real(d(j)./(z-z_northeast(j)));
    A(:,8*(j-1)+5)=imag(d(j)./(z-z_northeast(j)));
    A(:,8*(j-1)+6)=real(d(j)./(z-z_southwest(j)));
    A(:,8*(j-1)+7)=imag(d(j)./(z-z_southwest(j)));
    A(:,8*j)=real(d(j)./(z-z_southeast(j)));
    A(:,8*j+1)=imag(d(j)./(z-z_southeast(j)));
end
%Build the RHS 

fleft=zeros(2*n,1); fright=zeros(2*n,1);
fup=1-xup'.^2; fbottom=zeros(2*n,1);

f=[fleft; fbottom; fright; fup];

% Solve the least square problem using backslash (which uses QR
% factorization)
c= A\f;

x=0.99; y=0.99;
% Compute u(0.99,0.99)
u=c(1);
newz=x+y*1i;
for j=1:k
    unew=u+c(8*(j-1)+2)*real(d(j)./(newz-z_northwest(j)))+...
       +c(8*(j-1)+3)*imag(d(j)./(newz-z_northwest(j)))+...
    c(8*(j-1)+4)*real(d(j)./(newz-z_northeast(j)))+...
    c(8*(j-1)+5)*imag(d(j)./(newz-z_northeast(j)))+...
    c(8*(j-1)+6)*real(d(j)./(newz-z_southwest(j)))+...
    c(8*(j-1)+7)*imag(d(j)./(newz-z_southwest(j)))+...
    c(8*j)*real(d(j)./(newz-z_southeast(j)))+...
    c(8*j+1)*imag(d(j)./(newz-z_southeast(j)));
    u=unew;
end
fprintf('k= %3d  u=%10.8f\n' , k, u)
end
% The best estimate is u=0.01923508 when k=16
%% The ellipse function
function [b,c,d ]= ellipse(x,y)
% The function produce an ellipse centred at (0,0) with the equation 
% bx^2+cxy+dy^2=1 which best fits the input data (x,y)

% Build the n*3 least square matrix
A=[x.^2 x.*y y.^2];
% Build the RHS vector b
b=ones(length(x),1);
% Solve the least square problem using backslash (which uses QR
% factorization)
xnew= A\b;
b=xnew(1);
c=xnew(2);
d=xnew(3);

% work out the semi-major and semi-minor axes using b,c,d
% formula adapted from http://mathworld.wolfram.com/Ellipse.html
r1=sqrt((-c^2/2+2*b*d)/((c^2/4-b*d)*(sqrt((b-d)^2+c^2)-(b+d))));
r2=sqrt((-c^2/2+2*b*d)/((c^2/4-b*d)*(-sqrt((b-d)^2+c^2)-(b+d))));

% Build plotting points
theta=linspace(0,2*pi,100);
xx=r1*cos(theta);
yy=r2*sin(theta);
% Calculate the counterclockwise angle of rotation
% from the x-axis to the major axis of the ellipse 
if c==0 & b<d
    phi=0;
elseif c==0 & b>d
    phi=pi/2;
elseif c~=0 & b<d
    phi=acot((b-d)/c)/2;
else   c~=0 & b>d
    phi=pi/2+acot((b-d)/c)/2;
end
    

% Build the rotation matrix
R  = [cos(phi) -sin(phi); ...
      sin(phi)  cos(phi)];
% Compute the rotated coordinates
rCoords = R*[xx; yy];   
xr = rCoords(1,:)';      
yr = rCoords(2,:)'; 

% Plot the ellipse
hold off;
plot(xr,yr,'b');
hold on;
plot(x,y,'x')
title('The best fit ellipse')
end
