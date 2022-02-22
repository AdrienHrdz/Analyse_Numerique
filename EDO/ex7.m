clear variables;
close all;
clc;

m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
g = 9.81;

theta1_0 = pi/3;
theta1p_0 = 0;
theta2_0 = 2*pi/3;
theta2p_0 = 0;

tmin = 0;
tmax = 4;
h = 0.01;

f1 =@(t,theta1,theta2,z1,z2)( z1 );
f2 =@(t,theta1,theta2,z1,z2)( z2 );
g1 =@(t,theta1,theta2,z1,z2)( -(g*(2*m1 + m2)*sin(theta1) + m2*(g*sin(theta1-2*theta2) + 2*(l2*z2^2 + l1*z1^2*cos(theta1-theta2)))*sin(theta1-theta2) )/(2*l1*(m1 + m2*sin(theta1 - theta2).^2)) );
g2 =@(t,theta1,theta2,z1,z2)( sin(theta1-theta2)*((m1+m2)*(l1*z1^2+g*cos(theta1)) + l2*m2*z2^2*cos(theta1-theta2))/(l2*(m1+m2*sin(theta1-theta2)^2)) );


[theta1,theta1p,theta2,theta2p,t] = fct_RK4_4D(theta1_0,theta1p_0,theta2_0,theta2p_0,tmin,tmax,h,f1,f2,g1,g2);

figure(1)
for k=1:65:length(theta1)
    x= l1*sin(theta1(k));
    y= 1-l1*cos(theta1(k));
    plot([0,x],[0,y],'Marker','o','MarkerFacecolor','r','MarkerSize',10);
    axis('equal');
    axis(1.1*[xmin,xmax,ymin,ymax]);
    grid 'on';
    t1=title(strTitle);
    set(t1,'interpreter','latex');
    
    pause(0.0001);
end


