function [Xbest,Fbest] = NSDBO(params,MultiObj)
name=MultiObj.name;%问题名
numOfObj=MultiObj.numOfObj;%目标函数个数
evaluate_objective=MultiObj.fun;
D=MultiObj.nVar;
LB=MultiObj.var_min;
UB=MultiObj.var_max;
Max_iteration = params.maxgen;  % Set the maximum number of generation (GEN)
SearchAgents_no = params.Np;      % Set the population size (Search Agent) 
ishow = 1;
Nr=params.Nr;
chromosome = initialize_variables(SearchAgents_no, numOfObj, D, LB, UB,evaluate_objective);
intermediate_chromosome = non_domination_sort_mod(chromosome, numOfObj, D);
Pop = replace_chromosome(intermediate_chromosome, numOfObj,D,Nr);  
M=numOfObj;
K = D+M;
POS = Pop(:,1:K+1);
POS_ad = POS(:,1:K);
newPOS=zeros(SearchAgents_no,K);
DOMINATED= checkDomination(POS(:,D+1:D+M));
Pop  = POS(~DOMINATED,:);
ERP=Pop(:,1:K+1); % 
%% Optimization Circle
Iteration = 1;
pNum1=floor(SearchAgents_no*0.2);
pNum2=floor(SearchAgents_no*0.4);
pNum3=floor(SearchAgents_no*0.63);
while Iteration<=Max_iteration % for each generation   
    leng=size(ERP,1);
    r2=rand;
    for i = 1 : pNum1
        if(r2<0.9)
            r1=rand;
            a=rand;
            if (a>0.1)
                a=1;
            else
                a=-1;
            end
            worse=ERP(randperm(leng,1),1:D);
            newPOS(i,1:D) =  POS(i,1:D)+0.3*abs(POS(i,1:D)-worse)+a*0.1*POS_ad(i,1:D); % Equation (1)
        else
            aaa= randperm(180,1);
            if ( aaa==0 ||aaa==90 ||aaa==180 )
                newPOS(i,1:D) = POS(i,1:D);
            end
            theta= aaa*pi/180;
            newPOS(i,1:D) = POS(i,1:D)+tan(theta).*abs(POS(i,1:D)-POS_ad(i,1:D));    % Equation (2)
        end
    end
    R=1-Iteration/Max_iteration;    
    bestXX=ERP(randperm(leng,1),1:D);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xnew1 = bestXX.*(1-R);
    Xnew2 =bestXX.*(1+R);                    %%% Equation (3)
    Xnew1=bound(Xnew1,UB,LB); %越界判断
    Xnew2=bound(Xnew2,UB,LB); %越界判断
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bestX=ERP(randperm(leng,1),1:D);
    Xnew11 = bestX.*(1-R);
    Xnew22 =bestX.*(1+R);                     %%% Equation (5)
    Xnew11=bound(Xnew11,UB,LB); %越界判断
    Xnew22=bound(Xnew22,UB,LB); %越界判断
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = ( pNum1 + 1 ) :pNum2                  % Equation (4)
        newPOS(i,1:D)=bestXX+((rand(1,D)).*(POS(i,1:D)-Xnew1)+(rand(1,D)).*(POS(i,1:D)-Xnew2));
    end
    for i = pNum2+1: pNum3                  % Equation (6)
        newPOS(i,1:D)=POS(i,1:D)+((randn(1)).*(POS(i,1:D)-Xnew11)+((rand(1,D)).*(POS(i,1:D)-Xnew22)));
    end
    for j = pNum3+1 : SearchAgents_no                 % Equation (7)
        newPOS(i,1:D)=bestX+randn(1,D).*((abs(( POS(i,1:D)-bestXX)))+(abs(( POS(i,1:D)-bestX))))./2;
    end
    %% 计算函数值
    for i=1:SearchAgents_no
        newPOS(i,1:D)=bound(newPOS(i,1:D),UB,LB); %越界判断
        newPOS(i,D + 1: K) = evaluate_objective( newPOS(i,1:D));% 计算函数值
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if ( newPOS(i,D+k)<POS(i,D+k))
                dom_less = dom_less + 1;
            elseif (newPOS(i,D+k)== POS(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M
            POS_ad(i,1:K) = POS(i,1:K);
            POS(i,1:K) = newPOS(i,1:K);
        else
            POS_ad(i,1:K)= newPOS(i,1:K);
        end % end if
    end
  %%        
     pos_com = [POS(:,1:K) ; POS_ad];
     intermediate_pos = non_domination_sort_mod(pos_com, M, D);
      POS=replace_chromosome(intermediate_pos, M,D,Nr);
     DOMINATED= checkDomination(POS(:,D+1:D+M));
    Pop  = POS(~DOMINATED,:);
    ERP=Pop(:,1:K+1); %     
    %% 画图
    pl_data= ERP(:,D+1:D+M); % extract data to plot
    Fbest=sortrows(pl_data,2);
    PopFit=POS(:,D+1:D+M);
    Title = sprintf('迭代第 %d 次 , 存档库内非支配解个数 = %d \n可行解的个数=%d',Iteration,size(Fbest,1),size(PopFit,1));
%     PlotCosts(PopFit,Fbest,Title);
    disp(Title);
    %% 
    Iteration = Iteration+1;
end 
Xbest=ERP(:,1:D);
Fbest= ERP(:,D+1:D+M); % extract data to plot
end
% Check the boundary limit
function a=bound(a,ub,lb)
flagub=a>ub;
flaglb=a<lb;
idx=flaglb+flagub;
a=(~idx).*a+idx.*(rand(size(ub)).*(ub-lb)+lb);
end
function d = dominates(x,y)
    d = all(x<=y,2) & any(x<y,2);
end
function dom_vector = checkDomination(fitness)
    Np = size(fitness,1);
    dom_vector = zeros(Np,1);
    all_perm = nchoosek(1:Np,2);    % Possible permutations
    all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];
    
    d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
    dominated_particles = unique(all_perm(d==1,2));
    dom_vector(dominated_particles) = 1;
end
