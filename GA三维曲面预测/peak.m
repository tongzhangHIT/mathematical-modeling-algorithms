function  [xz,y,z] = peak(arg1,arg2)


if nargin == 0
    dx = 1/8;
    [x,y] = meshgrid(-3:dx:3);
elseif nargin == 1
    if length(arg1) == 1
        [x,y] = meshgrid(linspace(-3,3,arg1));
    else
        [x,y] = meshgrid(arg1,arg1);     
    end
else
    x = arg1; y = arg2;
end

%z =  3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ...
 %  - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ...
  % - 1/3*exp(-(x+1).^2 - y.^2);
z=-20*exp(-0.2*sqrt((x.^2+y.^2)/2))-exp((cos(2*pi*x)+cos(2*pi*y))/2)+20+2.71282;
if nargout > 1
    xz = x;
elseif nargout == 1
    xz = z;
else
    % Self demonstration
   
    surf(x,y,z)
    axis('tight')
    xlabel('x'), ylabel('y')
end
