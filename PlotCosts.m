
function PlotCosts(PopFit,Fbest,Title)
figure(1)
plot(PopFit(:,1),PopFit(:,2),'ro')
hold on
plot(Fbest(:,1),Fbest(:,2),'b*')
xlabel('目标函数1：运行能耗')
ylabel('目标函数2：运行时间')
grid on
hold off
title(Title)
legend('NSDBO的可行解' ,'存档库内非占优解')%,'location','best')


figure(10)
plot(Fbest(:,1),Fbest(:,2),'m*');
legend('NSDBO');
xlabel('运行能耗')
ylabel('运行时间')
title('pareto前沿解集')
end