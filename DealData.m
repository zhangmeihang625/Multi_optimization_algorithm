
function PG=DealData(Xbest,Fbest,idxn)
%% 比较不同目标函数寻优对调度结果的影响
if idxn==1
%% 第1种.将两个目标函数值归一化相加，取相加后最小的目标值的粒子，即寻找折衷解并画图
object=sum(Fbest./max(Fbest),2);
[~,idx]=min(object);
pg=Xbest(idx,:);
Title = sprintf('折衷解情况下:');
elseif idxn==2
%% 第2种寻找总成本最低时的解并画图
[~,idx]=min(sum(Fbest,2));
pg=Xbest(idx,:);
Title = sprintf('总成本最低情况下:');
elseif idxn==3
%% 第3种寻找运行成本最低时的解并画图
[~,idx]=min(Fbest(:,1));
pg=Xbest(idx,:);
Title = sprintf('运行能耗最低情况下:');
else
%% 第4种运行时间最低时的解并画图
[~,idx]=min(Fbest(:,2));
pg=Xbest(idx,:);
Title = sprintf('运行时间本最低情况下:');
end
%% 不同情况下的解赋值
 for i=1:24
   pg_SS(i)=pg(i);
 end  
 
 for m=25:48
    pg_FR(m-24)=pg(m);
end
for m=49:72
    pg_CF(m-48)=pg(m);
end
for m=73:96
    pg_CD(m-72)=pg(m);
end
for m=97:120
    pg_CW(m-96)=pg(m);
end
for m=121:144
    pg_FA(m-120)=pg(m);
end
PG.pg=pg;
PG.pg_SS=pg_SS;
PG.pg_FR=pg_FR;
PG.pg_CF=pg_CF;
PG.pg_CD=pg_CD;
PG.pg_CW=pg_CW;
PG.pg_FA=pg_FA;
PG.Title=Title;
end

