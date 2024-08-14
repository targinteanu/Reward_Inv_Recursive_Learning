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

%% setup reward based on observations 

ind1 = randi(length(data_all_trials)); 
dtas = data_all_trials{ind1}; 
ind2 = randi(length(dtas));
dta = dtas(ind2);

Rcue = dta.tgt_cond; 
Rtru = dta.rew_cond; 
pRc0 = mean(Rtru(Rcue == 0));
pRc1 = mean(Rtru(Rcue == 1));

% Force Low -> Low -> Low rew
MDP.R(2,1:3,1) = pRc0;
% Force High -> High -> High rew 
MDP.R(3,1:3,2) = pRc1;

% Choice -> Low -> Low rew
MDP.R(1,1:3,1) = pRc0;
% Choice -> High -> High rew 
MDP.R(1,1:3,2) = pRc1;

%% Q learning 
gamma = .9; epsilon = .9; depsilon = .1;

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
QTable = QTable{1};

%% visualize Q table 

    figure; 
    imagesc(QTable); colorbar;
    ylabel('State')
    yticks(1:4)
    yticklabels(StateNames)
    xlabel('Action')
    xticks(1:3)
    xticklabels(ActNames)
    title('Q Table')