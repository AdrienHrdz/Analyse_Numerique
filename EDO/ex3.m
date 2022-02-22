clear all;
clear variables;
clc;
close all;

%% fonctions
f=@(t,y,z,w)(z);
g=@(t,y,z,w)(w);
L=@(t,y,z,w)(t*exp(-t^2)-1/y);
%% initialisation
tmin= 0;
tmax= 8;
h=0.01; % pas de calcul
t=tmin:h:tmax;
y=zeros(size(t)); 
z=zeros(size(t));
w=zeros(size(t));
y(1) = 1;
z(1) = 0;
w(1) = 1;
%% boucle iterative
for k=1: (tmax-tmin)/h
    y(k+1) = y(k) + h*f(t(k),y(k),z(k),w(k));
    z(k+1) = z(k) + h*g(t(k),y(k),z(k),w(k));
    w(k+1) = w(k) + h*L(t(k),y(k),z(k),w(k));
end
figure(1);hold on;
plot(t,y)
plot(t,z)
plot(t,w)
