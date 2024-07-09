function [Xbest,Fbest]=NSGA3(params,MultiObj)
name=MultiObj.name;%问题名
numOfObj=MultiObj.numOfObj;%目标函数个数
CostFunction=MultiObj.fun;
nVar=MultiObj.nVar;
VarSize = [1 nVar]; % Size of Decision Variables Matrix
VarMin=MultiObj.var_min;
VarMax=MultiObj.var_max;
% Number of Objective Functions
nObj = numel(CostFunction(unifrnd(VarMin, VarMax, VarSize)));
%% NSGA-II Parameters

% Generating Reference Points
nDivision = 30;
Zr = GenerateReferencePoints(nObj, nDivision);

MaxIt = params.maxgen;  % Maximum Number of Iterations

nPop = params.Np;  % Population Size

pCrossover = 0.5;       % Crossover Percentage
nCrossover = 2*round(pCrossover*nPop/2); % Number of Parnets (Offsprings)

pMutation = 0.5;       % Mutation Percentage
nMutation = round(pMutation*nPop);  % Number of Mutants

mu = 0.02;     % Mutation Rate
sigma=0.2;% Mutation Step Size

%% Colect Parameters

params.nPop = nPop;
params.Zr = Zr;
params.nZr = size(Zr,2);
params.zmin = [];
params.zmax = [];
params.smin = [];

%% Initialization
empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.NormalizedCost = [];
empty_individual.AssociatedRef = [];
empty_individual.DistanceToAssociatedRef = [];

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    pop(i).Cost = CostFunction(pop(i).Position);
end

% Sort Population and Perform Selection
[pop, ~, params] = SortAndSelectPopulation(pop, params);


%% NSGA-II Main Loop

for it = 1:MaxIt
 
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for k = 1:nCrossover/2

        i1 = randi([1 nPop]);
        p1 = pop(i1);

        i2 = randi([1 nPop]);
        p2 = pop(i2);

        [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position);
        
      Flag4ub= popc(k, 1).Position>VarMax;
      Flag4lb= popc(k, 1).Position<VarMin;
     popc(k, 1).Position= popc(k, 1).Position.*(~(Flag4ub+Flag4lb))+(Flag4ub+Flag4lb).*(rand.*(VarMax-VarMin)+VarMin);
        popc(k, 1).Cost = CostFunction(popc(k, 1).Position);
        
      Flag4ub= popc(k, 2).Position>VarMax;
      Flag4lb= popc(k, 2).Position<VarMin;
     popc(k, 2).Position= popc(k, 2).Position.*(~(Flag4ub+Flag4lb))+(Flag4ub+Flag4lb).*(rand.*(VarMax-VarMin)+VarMin);
        popc(k, 2).Cost = CostFunction(popc(k, 2).Position);

    end
    popc = popc(:);

    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation

        i = randi([1 nPop]);
        p = pop(i);

        popm(k).Position = Mutate(p.Position, mu, sigma);
        
      Flag4ub= popm(k).Position>VarMax;
      Flag4lb= popm(k).Position<VarMin;
      popm(k).Position= popm(k).Position.*(~(Flag4ub+Flag4lb))+(Flag4ub+Flag4lb).*(rand.*(VarMax-VarMin)+VarMin);
      popm(k).Cost = CostFunction(popm(k).Position);

    end

    % Merge
    pop = [pop
           popc
           popm]; %#ok
    
    % Sort Population and Perform Selection
    [pop, F, params] = SortAndSelectPopulation(pop, params);
    
    % Store F1
    F1 = pop(F{1});

    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Repository size:' num2str(numel(F1))]);

    % Plot
   h_fig=figure(1);
   if it>1
       delete(h_par)
       if(isfield(MultiObj,'truePF'))
       delete(h_pf)
       end
   end
   for i=1:size(F1,1)
        POS_fit(i,:)=F1(i).Cost;
   end

         if(numOfObj==2)    
            h_par = plot(POS_fit(:,1),POS_fit(:,2),'or'); hold on;
            if(isfield(MultiObj,'truePF'))
                h_pf = plot(MultiObj.truePF(:,1),MultiObj.truePF(:,2),'.','color',0.8.*ones(1,3)); hold on;
            end
            title(name);
            grid on; xlabel('f1'); ylabel('f2');
            axis square;
        end
        if(numOfObj==3)
            h_par = plot3(POS_fit(:,1),POS_fit(:,2),POS_fit(:,3),'or'); hold on;
            if(isfield(MultiObj,'truePF'))
                h_pf = plot3(MultiObj.truePF(:,1),MultiObj.truePF(:,2),MultiObj.truePF(:,3),'.','color',0.8.*ones(1,3)); hold on;
            end
                       title(name);
            grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
            axis square;
        end

end
        delete(h_fig);
for i=1:size(F1,1)
    Xbest(i,:)=F1(i).Position;
    Fbest(i,:)=F1(i).Cost';
end
end



