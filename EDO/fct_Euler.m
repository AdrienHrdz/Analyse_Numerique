function [yEuler1,t1]=fct_Euler(y0,tmin,tmax,h,f)
    t1 = tmin:h:tmax;
    
    y = zeros(1,length(t1));
    y(1) = y0;
    for k = 1:length(t1)-1
        y(k+1) = y(k) + h*f(t1(k),y(k));
    end
    yEuler1 = y;
end 