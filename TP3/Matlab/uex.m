function uex = uex(x,y)
%UEX
%    UEX = UEX(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    20-Jun-2023 08:14:06

t2 = x.*2.0;
t3 = -y;
uex = t2+t3+1.0./(x+y.^3+1.0)+t2.*y+t3.*y+x.^2.*3.0+cos(x.*5.0).*sin(y.*7.0).*5.0-1.0;
