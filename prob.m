
%�õ���Ŀ������Ľ�
function [y,c] = prob(x)  %c=1��xΪ�ǿ��н�

%Ŀ�꺯��1�������ܺĸ�ֵ��y��1��
[c,y(1)] = fitness(x);  %�ó����ӵ���Ӧ��  

%% Ŀ�꺯��2:����ʱ�丳ֵ��y��2��
C_CD_en=0;C_FA_en=0;C_CW_en=0;
for i=1:144
    if i>72&&i<97
      C_CD_en=C_CD_en+(0.023*680+6*0.306+8*10.09)*x(i); 
    elseif  i>96&&i<121
     C_CW_en=C_CW_en+((0.023*889+6*1.8+8*1.6)*6*x(i));
    end  
end
for i=121:144
    if x(i)>0
        C_FA_en=C_FA_en+(0.023*724+6*0.0036+8*0.2)*x(i);
    else
        C_FA_en=C_FA_en+0;
    end
end
C_en= C_CD_en+ C_CW_en+ C_FA_en;
y(2) =C_en; %���������ɱ���ֵ��y��2��

end

% %Ŀ�꺯��3:���гɱ���ֵ��y��3��
% C_Ec=0;C_Tp=0;C_Ad=0;
% for i=1:144
%     if i>72&&i<97
%       C_Ec=C_Ec+(0.023*680+6*0.306+8*10.09)*x(i); 
%     elseif  i>96&&i<121
%      C_Tp=C_Tp+((0.023*889+6*1.8+8*1.6)*6*x(i));
%     end  
% end
% for i=121:144
%     if x(i)>0
%         C_Ad=C_Ad+(0.023*724+6*0.0036+8*0.2)*x(i);
%     else
%         C_Ad=C_Ad+0;
%     end
% end
% C_total= C_Ec+ C_Tp+ C_Ad;
% y(3) =C_total; %���������ɱ���ֵ��y��3��
% 
% end