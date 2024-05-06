%% setup MDP structure without reward
MDP = createMDP(["Choice"; "At Low"; "At High"], ["To Low"; "To High"]);
MDP.T(1,2,1) = 1; 
MDP.T(1,3,2) = 1; 
MDP.TerminalStates = ["At Low"; "At High"];

%% setup 
gamma = .5;
ind1 = randi(length(data_all_trials)); 

%% get muE of actual data 
dtas = data_all_trials{ind1}; 
muE = zeros(3, length(dtas)); 
for ind2 = 1:length(dtas)
    dta = dtas(ind2);
    Srec = []; Arec = []; 
    choicetrl = dta.task_cond == 0; % choice trials only 
    for trl = 1:(length(dta.task_cond))
        if choicetrl(trl)
            cueChosen = dta.choice(trl); 
            Arec = [Arec, cueChosen+1];
            Srec = [Srec, [0; ~cueChosen; cueChosen]]; % one-hot encoding 
        end
    end
    Gamma = gamma.^(0:(size(Srec,2)-1));
    muE(:,ind2) = Srec * Gamma';
end
muE = mean(muE, 2); 

%% init: randomize pi0 and get mu0
mu0 = zeros(3,length(dtas));
for ind2 = 1:length(dtas)
    Sapp = [1;0;0]; Strl = Sapp; % start 
    Aapp = [];
    for trl = 1:min(size(Srec,2),100)
        Atrl = randi(2); Aapp = [Aapp, Atrl]; 
        % 1 = low (cueChosen = 0); 2 = high (cueChosen = 1)
        cueChosen = Atrl-1; 
        if find(Strl) == 1
            % at start 
            Sapp = [Sapp, [0; ~cueChosen; cueChosen]];
        else
            % goes back to start 
            Sapp = [Sapp, 1; 0; 0];
        end
    end
    Gamma = gamma.^(0:(size(Sapp,2)-1));
    mu0(:,ind2) = Sapp * Gamma';
end
mu0 = mean(mu0, 2);

%% iterate until policy convergence 
Del = inf; theta = .001;
MT = [muE, mu0]; YT = [1, 0];
Qfuns = {}; Qtbls = {}; Dels = [Del];

    while Del > theta

        % step 2: get w, Del 
        SvmMdl = fitclinear(MT', YT'); 
        wT = SvmMdl.Beta'; wT = wT./norm(wT); % unit row vector 
        Del = wT*(muE - MT(:,2:end)); Del = min(Del) 

        % step 4: get pi, mu 
        %[Qfun, Qtbl, Srl] = ReinforcementLearnGrid(@(s) wT*unwrapPhi(phiGrid(s)), gamma, .8, 100);
        [Qtbl, Sind] = doRL(MDP, wT, gamma, .8, 0.1, 0.1);
        if Del <= min(Dels)
            Dels = [Dels, Del];
            Qtbls = [Qtbls; Qtbl];
        end
        Srl = zeros(3, length(Sind));
        for si = 1:length(Sind)
            Srl(Sind(si),si) = 1;
        end
        Gamma = gamma.^(0:(size(Srl,2)-1));
        MT = [MT, Srl*Gamma']; YT = [YT, 0];

    end

%%
function [QTable, StateSim] = doRL(MDP, w, gamma, alpha, epsilon, depsilon)
%% add reward to MDP 
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

%% Q learn 
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
trainOpts.MaxEpisodes = 500;
trainOpts.StopTrainingCriteria = "AverageReward";
trainOpts.StopTrainingValue = 13;
trainOpts.ScoreAveragingWindowLength = 30;

trainingStats = train(qAgent,env,trainOpts);

Data = sim(qAgent,env);
StateSim = squeeze(Data.Observation.MDPObservations.Data);
QTable = getLearnableParameters(getCritic(qAgent));
end