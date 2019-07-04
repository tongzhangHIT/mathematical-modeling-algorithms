function z=test_func(in)
nn=size(in);
x=in(:,1);
y=in(:,2);
nx=nn(1);
for i=1:nx
    temp=sin(sqrt(x(i)^2+y(i)^2))/sqrt(x(i)^2+y(i)^2)+exp((cos(2*pi*x(i))+cos(2*pi*y(i)))/2)-2.71289;
z(i,:)=temp;
end
