
function PlotCosts(PopFit,Fbest,Title)
figure(1)
plot(PopFit(:,1),PopFit(:,2),'ro')
hold on
plot(Fbest(:,1),Fbest(:,2),'b*')
xlabel('Ŀ�꺯��1�������ܺ�')
ylabel('Ŀ�꺯��2������ʱ��')
grid on
hold off
title(Title)
legend('NSDBO�Ŀ��н�' ,'�浵���ڷ�ռ�Ž�')%,'location','best')


figure(10)
plot(Fbest(:,1),Fbest(:,2),'m*');
legend('NSDBO');
xlabel('�����ܺ�')
ylabel('����ʱ��')
title('paretoǰ�ؽ⼯')
end