
function PG=DealData(Xbest,Fbest,idxn)
%% �Ƚϲ�ͬĿ�꺯��Ѱ�ŶԵ��Ƚ����Ӱ��
if idxn==1
%% ��1��.������Ŀ�꺯��ֵ��һ����ӣ�ȡ��Ӻ���С��Ŀ��ֵ�����ӣ���Ѱ�����ԽⲢ��ͼ
object=sum(Fbest./max(Fbest),2);
[~,idx]=min(object);
pg=Xbest(idx,:);
Title = sprintf('���Խ������:');
elseif idxn==2
%% ��2��Ѱ���ܳɱ����ʱ�ĽⲢ��ͼ
[~,idx]=min(sum(Fbest,2));
pg=Xbest(idx,:);
Title = sprintf('�ܳɱ���������:');
elseif idxn==3
%% ��3��Ѱ�����гɱ����ʱ�ĽⲢ��ͼ
[~,idx]=min(Fbest(:,1));
pg=Xbest(idx,:);
Title = sprintf('�����ܺ���������:');
else
%% ��4������ʱ�����ʱ�ĽⲢ��ͼ
[~,idx]=min(Fbest(:,2));
pg=Xbest(idx,:);
Title = sprintf('����ʱ�䱾��������:');
end
%% ��ͬ����µĽ⸳ֵ
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

