function y = wrap(x,x0,x1)

 v = (x-x0)/(x1-x0);
 v = v - floor(v);
 y = v*(x1-x0) + x0;
