
function MultiObj = GetFunInfo(TestProblem) %46个多目标测试函数
switch TestProblem
    case 1
        global CF; %切削力
        global FR;%进给速率
        global PV;%
        global Ec;%加工能耗
        global Tp;%加工时间
        %获取数据
        data=xlsread('CFRPdata.xlsx');
        SS=data(:,1);
        FR=data(:,2);
        CF=data(:,3);
        Ec=data(:,4);
        Tp=data(:,5);
        
        %CF最大值
        CFMax_dischar=40;
        %CF最小值
        CFMax_char=-40;
        %最大切削深度
        CDMax=2.5;
        %最小切削深度
        CDMin=0;
        %最大切削宽度
        CWMax=8;
        %最小切削宽度
        CWMin=0;
        %纤维方向角最大值
        FAMax=90;
        %纤维方向角最小值
        FAMin=0;
        
        %各设备出力约束
        for n=1:144 %粒子长度为144（SS，FR，CF，CD,CW，FA）
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