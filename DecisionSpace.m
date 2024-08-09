%% DecisionSpace
% Run through the apprentice algorithm for the MDP "decision space"
% formulation of the problem. This uses MATLAB's RL toolbox for Q learning.
% A helper function at the bottom manages the Q learning step. 
% This does not call external functions except for loading the data below. 

%data_all_trials = Load_Analyzed_Data();
load("data_all_trials.mat")

%% setup MDP structure without reward
StateNames = ["Choice"; "Force Low"; "Force High"]; 
ActNames = ["To Low"; "To High"; "Other"]; 
MDP = createMDP(StateNames, ActNames);

% Choice -> Low -> Choice or Force 
MDP.T(1,1:3,1) = [.5, .25, .25]; 

% Choice -> High -> Choice or Force
MDP.T(1,1:3,2) = [.5, .25, .25]; 

% Force Low -> Low -> Choice or Force 
MDP.T(2,1:3,1) = [.5, .25, .25]; 
% Force Low -> High -> Force Low 
MDP.T(2,2,2) = 1; 

% Force High -> High -> Choice or Force 
MDP.T(3,1:3,2) = [.5, .25, .25]; 
% Force High -> Low -> Force High 
MDP.T(3,3,1) = 1;

% Any -> Other -> Same 
MDP.T(1,1,3) = 1;
MDP.T(2,2,3) = 1;
MDP.T(3,3,3) = 1;

% MDP.TerminalStates = [];

%% setup 
gamma = .5;
epochdur = 5; % s
ind1 = randi(length(data_all_trials)); 

%% get muE of actual data 
dtas = data_all_trials{ind1}; 
muE = zeros(length(StateNames), length(dtas)); 
Ssess = []; Asess = [];
for ind2 = 1:length(dtas)
    dta = dtas(ind2);
    dtadur = cellfun(@(t)max(t)-min(t), dta.time);
    Srec = []; Arec = []; 
    choicetrl = dta.task_cond == 0; 
    for trl = 1:(length(dta.task_cond))
        nEpoch = round(dtadur(trl)/epochdur);
        if choicetrl(trl)
            cueChosen = dta.choice(trl); 
            Atrl = cueChosen+1;
            Strl = [1; 0; 0]; % one-hot encoding 
        else
            cueForced = dta.tgt_cond(trl);
            Atrl = cueForced+1;
            Strl = [0; ~cueForced; cueForced];
        end
        % add epochs of other action + no state change
        Arec = [Arec, 3*ones(1,nEpoch-1), Atrl];
        Srec = [Srec, repmat(Strl,1,nEpoch)];
    end
    Gamma = gamma.^(0:(size(Srec,2)-1));
    muE(:,ind2) = Srec * Gamma';
    Ssess = [Ssess, Srec]; Asess = [Asess, Arec];
end
muE = mean(muE, 2); 

%% init: randomize pi0 and get mu0
mu0 = zeros(length(StateNames),length(dtas));
for ind2 = 1:1000
    Sapp = strcmp(MDP.CurrentState, MDP.States); Strl = Sapp; % start 
    StrlIdx = find(Strl);
    Aapp = [];
    for trl = 1:100
        Atrl = randi(length(ActNames)); Aapp = [Aapp, Atrl]; 
        StrlOpts = MDP.T(StrlIdx,:,Atrl);
        StrlOpts = cumsum(StrlOpts);
        StrlIdx = rand < StrlOpts; 
        StrlIdx = find(StrlIdx); StrlIdx = StrlIdx(1); 
        Strl = zeros(size(StateNames)); Strl(StrlIdx) = 1;
        Sapp = [Sapp, Strl];
    end
    Gamma = gamma.^(0:(size(Sapp,2)-1));
    mu0(:,ind2) = Sapp * Gamma';
end
mu0 = mean(mu0, 2);

%% iterate until policy convergence 
Del = inf; theta0 = .001; theta = theta0;
MT = [muE, mu0]; YT = [1, 0];
Qfuns = {}; Qtbls = {}; Dels = [Del]; Ssim = {};

    while Del > theta

        % step 2: get w, Del 
        SvmMdl = fitclinear(MT', YT'); 
        wT = SvmMdl.Beta'; wT = wT./norm(wT); % unit row vector 
        Del = abs(wT*(muE - MT(:,2:end))); Del = min(Del) 

        % step 4: get pi, mu 
        [Qtbl, Sind] = doRL(MDP, wT, gamma, .8, 0.9, 0.1);
        if Del >= min(Dels)
            theta = theta*2;
        else
            theta = theta0;
        end
        %if Del < min(Dels)
            Dels = [Dels, Del];
            Qtbls = [Qtbls; Qtbl];
        %else
        Srl = zeros(length(StateNames), length(Sind));
        for si = 1:length(Sind)
            Srl(Sind(si),si) = 1;
        end
        Ssim = [Ssim; {Srl}];
        Gamma = gamma.^(0:(size(Srl,2)-1));
        MT = [MT, Srl*Gamma']; YT = [YT, 0];

    end

%% analysis 

% pick best of tested policies 
m = quadprog(eye(length(muE)), 0*muE, eye(length(muE)), muE);
mu2 = muE - m;
M = MT(:,2:end);
lambda = quadprog(M'*M, M'*muE, -eye(size(M,2)), zeros(size(M,2),1), ones(1,size(M,2)), 1);
mu = MT(:,2:end)*lambda;
Qtbl = Qtbl{1};
Qtbls = [rand(size(Qtbl)); Qtbls];
Qtbl = zeros(size(Qtbl));
for i = 1:length(Qtbls)
    Qtbl = Qtbl + lambda(i)*Qtbls{i};
end
Qtbls = [Qtbls; Qtbl];

%%
function [QTable, StateSim] = doRL(MDP, w, gamma, alpha, epsilon, depsilon)
% add reward to MDP 
% since phi is a one-hot encoding of state, w(s) = R(s) 
for s = 1:length(w)
    for s_orig = 1:length(MDP.States)
        for a = 1:length(MDP.Actions)
            if MDP.T(s_orig,s,a)
                MDP.R(s_orig,s,a) = w(s); 
            end
        end
    end
end

% Q learn 
env = rlMDPEnv(MDP);
env.ResetFcn = @() 1;

obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);
qTable = rlTable(obsInfo, actInfo);
qRepresentation = rlQValueRepresentation(qTable, obsInfo, actInfo);
qRepresentation.Options.LearnRate = gamma;

agentOpts = rlQAgentOptions;
agentOpts.DiscountFactor = gamma;
agentOpts.EpsilonGreedyExploration.Epsilon = epsilon;
agentOpts.EpsilonGreedyExploration.EpsilonDecay = depsilon;
% agentOpts.CriticOptimizerOptions = qOptions;
qAgent = rlQAgent(qRepresentation,agentOpts);

trainOpts = rlTrainingOptions;
trainOpts.MaxStepsPerEpisode = 50;
trainOpts.MaxEpisodes = 200;
trainOpts.StopTrainingCriteria = "AverageReward";
trainOpts.StopTrainingValue = 13;
trainOpts.ScoreAveragingWindowLength = 30;

trainingStats = train(qAgent,env,trainOpts);

Data = sim(qAgent,env);
StateSim = squeeze(Data.Observation.MDPObservations.Data);
QTable = getLearnableParameters(getCritic(qAgent));
end