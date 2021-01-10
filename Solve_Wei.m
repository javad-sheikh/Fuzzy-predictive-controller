function x=Solve_Wei(a,b,initials)
[T,x] = ode45(@Wei,[a:0.01:b],initials);