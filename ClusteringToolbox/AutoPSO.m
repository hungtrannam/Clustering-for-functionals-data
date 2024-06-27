%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML101
% Project Title: Evolutionary Automatic Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
clc;
clear;
close all;
%% Problem Definition
tic

addpath Data Package/Clust Package/ToolBox Package/Val

[Data, param.x, param.truelabels] = SimPDFAbnormal( ...
    { ...
    linspace(0, 1, 10*10), ...
    linspace(3, 4, 10)}, ...
    sqrt([.5, .5]));


X = Data';

% X= pdf';
k = 2;
CostFunction=@(s) ClusteringCost(s, X);     % Cost Function
VarSize=[k size(X,2)+1];  % Decision Variables Matrix Size
nVar=prod(VarSize);     % Number of Decision Variables
VarMin= repmat([min(X) 0],k,1);      % Lower Bound of Variables
VarMax= repmat([max(X) 1],k,1);      % Upper Bound of Variables


MaxIt=300;      % Maximum Number of Iterations
nPop=50;        % Population Size (Swarm Size)
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient
% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;
%% Initialization
empty_particle.Position      = [];
empty_particle.Cost          = [];
empty_particle.Out           = [];
empty_particle.Velocity      = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost     = [];
empty_particle.Best.Out      = [];
particle                     = repmat(empty_particle,nPop,1);
BestSol.Cost                 = inf;

%%%%%%%%%%%%%%%%%%%%%
% Clustering

for i=1:nPop
    
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    [particle(i).Cost, particle(i).Out]=CostFunction(particle(i).Position);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    particle(i).Best.Out=particle(i).Out;
    
    % Update Global Best
    if particle(i).Best.Cost<BestSol.Cost
        
        BestSol=particle(i).Best;
        
    end
    
end
BestCost=zeros(MaxIt,1);
%% PSO Main Loop
for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(BestSol.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        [particle(i).Cost, particle(i).Out] = CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            particle(i).Best.Out=particle(i).Out;
            
            % Update Global Best
            if particle(i).Best.Cost<BestSol.Cost
                
                BestSol=particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it)=BestSol.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    % Plot Solution
    figure(1);
    PlotSolution(X, BestSol);
    pause(0.01);
    
    w=w*wdamp;
    
end



%%%%%%%%%%%%%%%%%%%%%
% Results

figure;
plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


h = PlotPDFeachIteration(Data,BestSol.Out.ind', param.x);