%% ApprenticeRL 
% Run through the apprentice algorithm for the 2D grid-world state space
% formulation of the problem. Does not use MATLAB's RL toolbox. 
% This calls various external functions with "State" or "Grid" in the
% function name, and uses the external function "ReinforcementLearnGrid" 
% to manage the Q learning step. 

data_all_trials = Load_Analyzed_Data();

gamma = .5; % discount factor 
%{
phi = @(s) 1./(1 + exp(-s/20)); % sigmoid mapping to [0,1]
inv_phi = @(ph) 20*log(ph./(1-ph)); % inv sigmoid mapping to [-inf, inf]
%}
unwrapPhi = @(phiTbl) phiTbl{:,:}';

% loop through all 
ind1 = randi(length(data_all_trials)); 
%for ind1 = 1:(length(data_all_trials))
    
    %% get muE of actual data 
    dtas = data_all_trials{ind1};
%    ind2 = randi(length(dtas));
%    muE = 0;
    muE = zeros(8, length(dtas)); 
    for ind2 = 1:(length(dtas))
        dta = dtas(ind2);
        Srec = []; Arec = []; 
        trlct = 0;
        for trl = 1:(length(dta.task_cond))
            [Strl, Atrl] = getStateSpace(dta,trl,false);
            Strl = gridifyState(Strl); Atrl = gridifyState(Atrl);
            Srec = [Srec; Strl]; 
            Arec = [Arec; Atrl];

            if ~isempty(Strl)
                % Phi = phi(Strl{1:height(Strl), 1:width(Strl)}');
                Phi = phiGrid(Strl); Phi = Phi{:,:}';
                Gamma = gamma.^(0:(height(Strl)-1));
%               muE = muE + Phi*Gamma';
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

    %% init: randomize pi0 and get mu0 
    mu0 = zeros(8, length(dtas));
    for ind2 = 1:(length(dtas))
        Sapp = Srec(randi(height(Srec)),:); Strl = Sapp; 
        Aapp = []; Atrl = Arec(1,:);
        for trl = 1:min(height(Srec), 100)
            Atrl.eye_px_filt_trl = rand; 
            Atrl.eye_py_filt_trl = rand;
            Atrl = phiGridInv(Atrl);
            Strl = updateGridState(Strl, Atrl); 
            Sapp = [Sapp; Strl]; Aapp = [Aapp; Atrl];
        end
        % Phi = phi(Sapp{1:height(Sapp), 1:width(Sapp)}');
        Phi = phiGrid(Sapp); Phi = Phi{:,:}';
        Gamma = gamma.^(0:(height(Sapp)-1));
        mu0(:,ind2) = Phi*Gamma';
    end
    mu0 = mean(mu0,2);
    clear Sapp Aapp Phi Gamma 

    %% iterate until policy convergence 
    ind3 = 1; Del = inf; theta = .01; 
    MT = [muE, mu0]; YT = [1, 0]; 
    Qfuns = {}; Qtbls = {}; Dels = [Del];

    while Del > theta

        % step 2: get w, Del 
        SvmMdl = fitclinear(MT', YT'); 
        wT = SvmMdl.Beta'; wT = wT./norm(wT); % unit row vector 
        Del = wT*(muE - MT(:,2:end)); Del = min(Del) 
        clear SvmMdl

        % step 4: get pi, mu 
        [Qfun, Qtbl, Srl] = ReinforcementLearnGrid(wT, gamma, .8, 100);
        if Del <= min(Dels)
            Dels = [Dels, Del];
            Qfuns = {Qfuns; Qfun}; Qtbls = [Qtbls; {Qtbl}];
        end
        Phi = phiGrid(Srl); Phi = Phi{:,:}';
        Gamma = gamma.^(0:(height(Srl)-1));
        MT = [MT, Phi*Gamma']; YT = [YT, 0];

    end
%end