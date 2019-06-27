function [x y] = drawCurve(s, e, r)
    x0 = 0;
    y0 = 0;
    a1 = 2*pi*degtorad(s);  %start
    a2 = a1 + degtorad(e);
    t = linspace(a1,a2);
    x = x0 + r*cos(t);
    y = y0 + r*sin(t);