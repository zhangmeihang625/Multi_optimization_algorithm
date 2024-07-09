
function MultiObj = GetFunInfo(TestProblem) %46����Ŀ����Ժ���
switch TestProblem
    case 1
        global CF; %������
        global FR;%��������
        global PV;%
        global Ec;%�ӹ��ܺ�
        global Tp;%�ӹ�ʱ��
        %��ȡ����
        data=xlsread('CFRPdata.xlsx');
        SS=data(:,1);
        FR=data(:,2);
        CF=data(:,3);
        Ec=data(:,4);
        Tp=data(:,5);
        
        %CF���ֵ
        CFMax_dischar=40;
        %CF��Сֵ
        CFMax_char=-40;
        %����������
        CDMax=2.5;
        %��С�������
        CDMin=0;
        %����������
        CWMax=8;
        %��С�������
        CWMin=0;
        %��ά��������ֵ
        FAMax=90;
        %��ά�������Сֵ
        FAMin=0;
        
        %���豸����Լ��
        for n=1:144 %���ӳ���Ϊ144��SS��FR��CF��CD,CW��FA��
            if n<25
                lower_bound(n)=0;
                upper_bound(n) =SS(n);
            end
            if n>24&&n<49
                lower_bound(n)=0;
                upper_bound(n) =FR(n-24);
            end
            if n>48&&n<73
                lower_bound(n)=CFMax_char;
                upper_bound(n) =CFMax_dischar;
            end
            if n>72&&n<97
                lower_bound(n)=CDMin;
                upper_bound(n) =CDMax;
            end
            if n>96&&n<121
                lower_bound(n)=CWMin;
                upper_bound(n) =CWMax;
            end
            if n>120
                lower_bound(n)=FAMin;
                upper_bound(n) =FAMax;
            end
        end
        CostFunction = @Fun;
        nVar = 144;
        VarMin = lower_bound;
        VarMax = upper_bound;
        name='parameter optimization';
        numOfObj = 2;
end
MultiObj.nVar=nVar;
MultiObj.var_min = VarMin;
MultiObj.var_max =VarMax;
MultiObj.fun=CostFunction;
MultiObj.numOfObj=numOfObj;
MultiObj.name=name;
end
function o=Fun(x)
o=prob(x) ;
end