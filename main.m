close all;
clear ; 
clc;
addpath('./NSDBO/')%��ӷ�֧������������Ż��㷨·��  
addpath('./NSGA3/')%��ӷ�֧���Ŵ��㷨3·�� 
addpath('./MOGWO/')%��ӻ����Ż��㷨·��
addpath('./NSWOA/')%��ӷ�֧������ľ����Ż��㷨·��
addpath('./MOPSO/')%��Ӷ�Ŀ������Ⱥ�Ż��㷨·��

global SS;%����ת��
global FR;%�����ٶ�
global CF; %������

%%
TestProblem=1;
MultiObj = GetFunInfo(TestProblem);
MultiObjFnc=MultiObj.name;%������
% Parameters
params.Np =200;        %  ��Ⱥ��С(�����޸�)
params.Nr =200 ; % ���ⲿ�浵�Ĵ�С��
params.maxgen =200;    % ����������(�����޸�)
[Xbest{1},Fbest{1}] = NSDBO(params,MultiObj);
[Xbest{2},Fbest{2}] = NSGA3(params,MultiObj);
[Xbest{3},Fbest{3}]= MOGWO(params,MultiObj);
[Xbest{4},Fbest{4}] = NSWOA(params,MultiObj);
[Xbest{5},Fbest{5}] = MOPSO(params,MultiObj);

%% �Ƚϲ�ͬĿ�꺯��Ѱ�ŶԵ��Ƚ����Ӱ��:
% idxn=1:��1��.������Ŀ�꺯��ֵ��һ����ӣ�ȡ��Ӻ���С��Ŀ��ֵ�����ӣ���Ѱ�����Խ�
% idxn=2:��2��Ѱ���ܳɱ����ʱ�Ľ�
% idxn=3:��3��Ѱ�������ܺ����ʱ�Ľ�
% idxn=4:��4������ʱ�����ʱ�Ľ�
idxn=3;

%% �������
for i=1:size(Xbest,2)
PG{i}=DealData(Xbest{i},Fbest{i},idxn);
end
strColor={'r*','go','b<','k>','mp','c.','y*'};
strColor1={'r*-','go--','b<-','k>-','mp-','c-.','y-*'};
AlgorithmName={'NSDBO','NSGA3','MOGWO','NSWOA','MOPSO'};%�㷨����


%% �����ͼ
figure(1)
for  i=1:size(Fbest,2)
plot(Fbest{1,i}(:,1),Fbest{1,i}(:,2),strColor{i});
hold on
end
legend(AlgorithmName);
xlabel('time/s')
ylabel('energy/J')
saveas(gcf,'./Picture/ParetoFont.jpg') %��ͼƬ���浽Picture�ļ�������

%% ��ͼ
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
saveas(gcf,'./Picture/SS.jpg') %��ͼƬ���浽Picture�ļ�������

%% ��ͼ
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
saveas(gcf,'./Picture/FR.jpg') %��ͼƬ���浽Picture�ļ�������

%% ��ͼ
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
saveas(gcf,'./Picture/CD.jpg') %��ͼƬ���浽Picture�ļ�������

%% ��ͼ
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
saveas(gcf,'./Picture/CF.jpg') %��ͼƬ���浽Picture�ļ�������

%% ��ͼ
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
saveas(gcf,'./Picture/CW.jpg') %��ͼƬ���浽Picture�ļ�������

%% ��ͼ
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
saveas(gcf,'./Picture/FA.jpg') %��ͼƬ���浽Picture�ļ�������

%% ��ͼ
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
saveas(gcf,['./Picture/' AlgorithmName{i} '1.jpg']) %��ͼƬ���浽Picture�ļ�������

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
saveas(gcf,['./Picture/' AlgorithmName{i} '2.jpg']) %��ͼƬ���浽Picture�ļ�������
end


%% �����Ż�ʱ��ͳ��
 tic;
 elapsedTime = toc;
 fprintf('�Ż�ʱ�䣺%f ��\n', elapsedTime);

