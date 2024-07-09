close all;
clear ; 
clc;
addpath('./NSDBO/')%添加非支配排序的蜣螂优化算法路径  
addpath('./NSGA3/')%添加非支配遗传算法3路径 
addpath('./MOGWO/')%添加灰狼优化算法路径
addpath('./NSWOA/')%添加非支配排序的鲸鱼优化算法路径
addpath('./MOPSO/')%添加多目标粒子群优化算法路径

global SS;%主轴转速
global FR;%进给速度
global CF; %切削力

%%
TestProblem=1;
MultiObj = GetFunInfo(TestProblem);
MultiObjFnc=MultiObj.name;%问题名
% Parameters
params.Np =200;        %  种群大小(可以修改)
params.Nr =200 ; % （外部存档的大小）
params.maxgen =200;    % 最大迭代次数(可以修改)
[Xbest{1},Fbest{1}] = NSDBO(params,MultiObj);
[Xbest{2},Fbest{2}] = NSGA3(params,MultiObj);
[Xbest{3},Fbest{3}]= MOGWO(params,MultiObj);
[Xbest{4},Fbest{4}] = NSWOA(params,MultiObj);
[Xbest{5},Fbest{5}] = MOPSO(params,MultiObj);

%% 比较不同目标函数寻优对调度结果的影响:
% idxn=1:第1种.将两个目标函数值归一化相加，取相加后最小的目标值的粒子，即寻找折衷解
% idxn=2:第2种寻找总成本最低时的解
% idxn=3:第3种寻找运行能耗最低时的解
% idxn=4:第4种运行时间最低时的解
idxn=3;

%% 结果处理
for i=1:size(Xbest,2)
PG{i}=DealData(Xbest{i},Fbest{i},idxn);
end
strColor={'r*','go','b<','k>','mp','c.','y*'};
strColor1={'r*-','go--','b<-','k>-','mp-','c-.','y-*'};
AlgorithmName={'NSDBO','NSGA3','MOGWO','NSWOA','MOPSO'};%算法名称


%% 画结果图
figure(1)
for  i=1:size(Fbest,2)
plot(Fbest{1,i}(:,1),Fbest{1,i}(:,2),strColor{i});
hold on
end
legend(AlgorithmName);
xlabel('time/s')
ylabel('energy/J')
saveas(gcf,'./Picture/ParetoFont.jpg') %将图片保存到Picture文件夹下面

%% 画图
figure
for i=1:size(PG,2)
    plot(PG{1,i}.pg_SS,strColor1{i})
    hold on;
end
plot(SS,'-c*')
xlim([1 24])
legend([AlgorithmName 'actural power' ]);
xlabel('time/ms')
ylabel('energy/J')
title([PG{1,1}.Title  'SS'])
saveas(gcf,'./Picture/SS.jpg') %将图片保存到Picture文件夹下面

%% 画图
figure
for i=1:size(PG,2)
    plot(PG{1,i}.pg_FR,strColor1{i})
    hold on;
end
plot(FR,'-c*')
xlim([1 24])
legend([AlgorithmName 'actural power' ]);
xlabel('time/s')
ylabel('energy/J')
title([PG{1,1}.Title  'FR'])
saveas(gcf,'./Picture/FR.jpg') %将图片保存到Picture文件夹下面

%% 画图
figure
for i=1:size(PG,2)
    plot(PG{1,i}.pg_CD,strColor1{i})
    hold on;
end
xlim([1 24])
legend(AlgorithmName);
xlabel('time/s')
ylabel('energy/J')
title([PG{1,1}.Title  'CD'])
saveas(gcf,'./Picture/CD.jpg') %将图片保存到Picture文件夹下面

%% 画图
figure
for i=1:size(PG,2)
    plot(PG{1,i}.pg_CF,strColor1{i})
    hold on;
end
xlim([1 24])
legend(AlgorithmName);
xlabel('time/ms')
ylabel('energy/J')
title([PG{1,1}.Title  'CF'])
saveas(gcf,'./Picture/CF.jpg') %将图片保存到Picture文件夹下面

%% 画图
figure
for i=1:size(PG,2)
    plot(PG{1,i}.pg_CW,strColor1{i})
    hold on;
end
xlim([1 24])
legend(AlgorithmName);
xlabel('time/ms')
ylabel('energy/J')
title([PG{1,1}.Title  'CW'])
saveas(gcf,'./Picture/CW.jpg') %将图片保存到Picture文件夹下面

%% 画图
figure
for i=1:size(PG,2)
    plot(PG{1,i}.pg_FA,strColor1{i})
    hold on;
end
xlim([1 24])
legend(AlgorithmName);
xlabel('time/s')
ylabel('energy/J')
title([PG{1,1}.Title  'FA'])
saveas(gcf,'./Picture/FA.jpg') %将图片保存到Picture文件夹下面

%% 画图
for i=1:size(PG,2)
figure
plot(PG{1,i}.pg_FA,'-k<')
hold on
plot(PG{1,i}.pg_CD,'-ro')
hold on
plot(PG{1,i}.pg_CF,'-mp');
hold on
plot(PG{1,i}.pg_CW,'-c>')
legend('FA','CD','CF','CW');
xlabel('time/s')
ylabel('energy/J')
title([PG{1,i}.Title  AlgorithmName{i}])
saveas(gcf,['./Picture/' AlgorithmName{i} '1.jpg']) %将图片保存到Picture文件夹下面

figure
plot(PG{1,i}.pg_SS,'-rd')
hold on
plot(PG{1,i}.pg_FR,'-g*');
hold on
plot(CF,'-bo');
legend('SS','FR','CF');
xlabel('time/s')
ylabel('energy/J')
title([PG{1,i}.Title  AlgorithmName{i}])
saveas(gcf,['./Picture/' AlgorithmName{i} '2.jpg']) %将图片保存到Picture文件夹下面
end


%% 增加优化时间统计
 tic;
 elapsedTime = toc;
 fprintf('优化时间：%f 秒\n', elapsedTime);

