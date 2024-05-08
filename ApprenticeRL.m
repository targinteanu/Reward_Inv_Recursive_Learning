%% ApprenticeRL 
% Run through the apprentice algorithm for the 2D grid-world state space
% formulation of the problem. Does not use MATLAB's RL toolbox. 
% This calls various external functions with "State" or "Grid" in the
% function name, and uses the external function "ReinforcementLearnGrid" 
% to manage the Q learning step. 

disp('Accessing data...')
data_all_trials = Load_Analyzed_Data();

gamma = .5; % discount factor 

ind1 = randi(length(data_all_trials)); 
    
%% get muE of actual data 
dtas = data_all_trials{ind1};
muE = zeros(8, length(dtas)); 
for ind2 = 1:(length(dtas))
    dta = dtas(ind2);
    disp(['Analyzing ',dta.name,' ...']);
    Srec = []; Arec = []; % recorded state/action 
    trlct = 0;
    for trl = 1:(length(dta.task_cond))
        [Strl, Atrl] = getStateSpace(dta,trl,false);
        Strl = gridifyState(Strl); Atrl = gridifyState(Atrl);
        Srec = [Srec; Strl]; 
        Arec = [Arec; Atrl];

        if ~isempty(Strl)
            Phi = phiGrid(Strl); Phi = Phi{:,:}';
            Gamma = gamma.^(0:(height(Strl)-1));
            muE(:,ind2) = muE(:,ind2) + Phi*Gamma'; 
            trlct = trlct + 1;
        end
    end

    if trlct
        muE(:,ind2) = muE(:,ind2)/trlct;
    end

    Srec.Time = Srec.Time - Srec.Time(1);
    Arec.Time = Arec.Time - Arec.Time(1); 
end
muE = mean(muE,2);
clear dta trlct Phi Gamma 

%% load items to construct MDP and Q table 
load('StateQT.mat');

%% init apprentice: randomize pi0 and get mu0 
disp('Initializing Apprentice...')
mu0 = zeros(8, length(dtas));
for ind2 = 1:(length(dtas))
    Sapp = Srec(randi(height(Srec)),:); Strl = Sapp; % apprentice state
    Aapp = []; Atrl = Arec(1,:); % apprentice action 
    for trl = 1:min(height(Srec), 100)
        Atrl.eye_px_filt_trl = rand; 
        Atrl.eye_py_filt_trl = rand;
        Atrl = phiGridInv(Atrl);
        Strl = updateGridState(Strl, Atrl); 
        Sapp = [Sapp; Strl]; Aapp = [Aapp; Atrl];
    end
    Phi = phiGrid(Sapp); Phi = Phi{:,:}';
    Gamma = gamma.^(0:(height(Sapp)-1));
    mu0(:,ind2) = Phi*Gamma';
end
mu0 = mean(mu0,2);
clear Sapp Aapp Phi Gamma 

%% iterate until policy convergence 
Del = inf; theta = .001; 
MT = [muE, mu0]; YT = [1, 0]; 
Qtbls = []; Dels = [Del];

while Del > theta

    % step 2: get w, Del 
    SvmMdl = fitclinear(MT', YT'); 
    wT = SvmMdl.Beta'; wT = wT./norm(wT); % unit row vector 
    Del = abs(wT*(muE - MT(:,2:end))); Del = min(Del) 
    clear SvmMdl

    % init Q table 
    Qtbl = rand(length(Qstate), length(Qaction)); 
    Qtbl = single(Qtbl);

    % step 4: get pi, mu 
    [~, Qtbl2, Srl] = ReinforcementLearnGrid3(wT, gamma, .8, 10, L, Qstate, Qaction, Qtbl);
    if Del >= min(Dels)
        theta = theta*5;
    end
    %if Del <= min(Dels)
        Dels = [Dels, Del];
        Qtbls = cat(3,Qtbls,Qtbl2);
    %end
    Phi = phiGrid(Srl); Phi = Phi{:,:}';
    Gamma = gamma.^(0:(height(Srl)-1));
    MT = [MT, Phi*Gamma']; YT = [YT, 0];

end

%% pick best of tested policies 
m = quadprog(eye(length(muE)), 0*muE, eye(length(muE)), muE);
mu2 = muE - m;
M = MT(:,2:end);
lambda = quadprog(M'*M, M'*muE, -eye(size(M,2)), zeros(size(M,2),1), ones(1,size(M,2)), 1);
mu = MT(:,2:end)*lambda;
Qtbl2 = zeros(size(Qtbl2));
for ind3 = 1:size(Qtbls,3)
    Qtbl2 = Qtbl2 + lambda(ind3)*Qtbls(:,:,ind3);
end

%% pick best action at each state under this policy 
Abest = zeros(size(Qstate)); 
for ind3 = 1:length(Abest)
    [~,Abest(ind3)] = max(Qtbl2(ind3,:));
end

%% test accuracy of this policy against training data 
disp('Computing training accuracy...')
train_result = false(height(Srec),1);
for t = 1:height(Srec)
    St = Srec(t,:); At = Arec(t,:);
    StInd = state2ind(St,L); 
    if ~isEqualState(St, Qstate{StInd})
        error('State does not match index of table.')
    end
    Apred = Qaction{Abest(StInd)};
    train_result(t) = isEqualState(Apred, At);
end
training_accuracy = mean(train_result)

%% test accuracy of this policy against other data
disp('Computing testing accuracy...')
indTest = randi(length(data_all_trials)); 
dtasTest = data_all_trials{ind1};
test_result = logical([]);
for ind2 = 1:(length(dtasTest))
    dta = dtasTest(ind2);
    disp(['Testing ',dta.name,' ...']);
    for trl = 1:(length(dta.task_cond))
        [Strl, Atrl] = getStateSpace(dta,trl,false);
        Strl = gridifyState(Strl); Atrl = gridifyState(Atrl);
        if ~isempty(Strl)
            for t = 1:height(Strl)
                St = Strl(t,:); At = Atrl(t,:);
                StInd = state2ind(St,L);
                if ~isEqualState(St, Qstate{StInd})
                    error('State does not match index of table.')
                end
                Apred = Qaction{Abest(StInd)};
                test_result = [test_result; isEqualState(Apred, At)];
            end
        end
    end
end
testing_accuracy = mean(test_result)
clear dta 