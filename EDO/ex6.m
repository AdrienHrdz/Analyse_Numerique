clear variables;
clc;
close all;

%% param
G = 4*pi^2;
M = 1;
a = 1;
T = 1;
tmin = 0;
tmax = 3*T;
h = 1/365;

%% cond init
x0 = 0.5;
y0 = 0;
xp0 = 0;
yp0 = 11.5;

f =@(x,y,z,w)( z );
g =@(x,y,z,w)( -G*M*x/((x.^2 + y.^2).^(3/2)) );
l =@(x,y,z,w)( w );
m =@(x,y,z,w)( -G*M*y/((x.^2 + y.^2).^(3/2)) );

[x,xp,t] = RK4_2D(x0,xp0,tmin,tmax,h,f,g);
[y,yp,t] = RK4_2D(y0,yp0,tmin,tmax)





