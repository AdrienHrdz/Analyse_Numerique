function [x,y,t]=RK2_2D(x0,y0,tmin,tmax,pas,beta,F,G)
    t = tmin:pas:tmax;
    x = zeros(size(t));
    y = zeros(size(t));
    x(1) = x0;
    y(1) = y0;
    for k=1: (tmax-tmin)/pas
        k1x = F(t(k),x(k),y(k));
        k1y = G(t(k),x(k),y(k));
        k2x = F(t(k) + pas/(2*beta),x(k) + k1x*pas/(2*beta),y(k) + k1y*pas/(2*beta));
        k2y = G(t(k) + pas/(2*beta),x(k) + k1x*pas/(2*beta),y(k) + k1y*pas/(2*beta));
        x(k+1) = x(k) + pas*((1-beta)*k1x + beta*k2x);
        y(k+1) = y(k) + pas*((1-beta)*k1y + beta*k2y);
    end
end