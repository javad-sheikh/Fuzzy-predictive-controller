function dx=Wei(t,x)
dx = zeros(3,1);  
dx(1)=-x(2,1);
dx(2)=x(3,1)+x(1,1);
dx(3)=2*x(2,1).*x(2,1)+x(1,1).*x(3,1)-0.35; 