clc, clear all
close all
warning off

initials=[-1.6 0.82 1.9];
% N=2000;
N=400;
% runtime=4000;
runtime=700;

k=N+1;
xd=.5;
sumE=0;
emax=5;
br=1;
ss=0.1;
% ss=1;
% % endtime=0.01;
endtime=0.05;
m=Solve_Wei(0,endtime,initials);
hi(1)=length(m);
x(1)=m(1,1);
initials=m(end,:);

%% Step 1

for i=1:N
    y=Solve_Wei(i*endtime,i*endtime+endtime,initials);
    x(i+1)=y(end,1);
    hi(i+1)=hi(i)+length(y(:,1));
    m(hi(i)+1:hi(i+1),:)=y;
    initials=y(end,:);
end
xx=x;
sigma=var(x)*ss;

%% Step2

while br==1
    sum1=0;
    sum2=0;
    xts = [x(1:k-1);x(2:k)];
    ab = 1;
    for j=2:k-ab
        sum1=sum1+x(j).*exp(-norm(xts(:,k-ab)-xts(:,j-1))^2/(2*sigma^2));
        sum2=sum2+exp(-norm(xts(:,k-ab)-xts(:,j-1))^2/(2*sigma^2));
    end

    xx(k+1)=sum1/sum2;
    y=Solve_Wei(k*endtime,k*endtime+endtime,initials);
    x(k+1)=y(end,1);
    hi(k+1)=hi(k)+length(y(:,1));
    m(hi(k)+1:hi(k+1),:)=y;
    initials=y(end,:);
    
    %% Step3
    
    error(k)=xx(k+1)-x(k+1);
    if abs(error(k))<emax
        %% step4
        ggg = 1;
        U(k)=xd-xx(k+1)+((1/(k-1))*(sumE))*ggg;
        sumE=sumE+error(k);
        initials=initials+U(k);
        k=k+1;
    else
        %% step5
        N=N+1;
        initials=[-1.6 0.82 1.9];
        clear hi
        clear x
        clear xx
        sumE=0;
        clear error
        m=Solve_Wei(0,endtime,initials);
        hi(1)=length(m);
        x(1)=m(1,1);
        initials=m(end,:);
        for i=1:N
            y=Solve_Wei(i*endtime,i*endtime+endtime,initials);
            x(i+1)=y(end,1);
            hi(i+1)=hi(i)+length(y(:,1));
            m(hi(i)+1:hi(i+1),:)=y;
            initials=y(end,:);
            lll = hi(i);
        end
        xx=x;
        sigma=var(x)*ss;
        %% Step6
        k=N+1;
    end
    if k>runtime
        br=0;
    end
end

low=1;
for r=1:k-6
    for j=low:hi(r)
        X(j,1)=m(j,1)+U(r);
    end
    low=hi(r);
end
figure
plot([X(1:lll)],'b')
hold on
plot(m(1:lll,1),'r')
title('Controlled value')
legend('Controlled value of Wei','Real Value of Wei','Location','northwest')
