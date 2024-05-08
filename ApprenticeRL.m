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
    Srec = []; Arec = []; 
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

%% init apprentice: randomize pi0 and get mu0 
disp('Initializing Apprentice...')
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
    Phi = phiGrid(Sapp); Phi = Phi{:,:}';
    Gamma = gamma.^(0:(height(Sapp)-1));
    mu0(:,ind2) = Phi*Gamma';
end
mu0 = mean(mu0,2);
clear Sapp Aapp Phi Gamma 

%% init empty Q table 

[~,posOpts] = gridifyState();
posOpts = [-inf,posOpts,inf];
L = length(posOpts); % position takes L possible values from 0 to L-1 (base L digit) 

A0 = timetable(seconds(0), 0, 0); % null action 
A0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl'};
A0 = gridifyState(A0); 

S0 = timetable(seconds(0), 0,0, 0,0, inf,inf, inf,inf); % default state
S0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl', ...
    'tgt_px_fx', 'tgt_py_fx', 'tgt_px_lo', 'tgt_py_lo', 'tgt_px_hi', 'tgt_py_hi'};
S0 = gridifyState(S0);

disp('Getting action space...')

Qaction = repmat({A0}, 1, L^2); % all possible actions 
for eyeX = 0:(L-1)
    for eyeY = 0:(L-1)
        A1 = A0; 
        A1.eye_px_filt_trl = eyeX; 
        A1.eye_py_filt_trl = eyeY; 
        ind = act2ind(A1,L);
        if ind > length(Qaction)
            error('Action exceeds table size.')
        end
        Qaction{ind} = A1;
    end
end

disp('Getting state space...')

Qstate = repmat({S0}, 1, L^6); % all possible states 
%for eyeX = 0:(L-1)
%    for eyeY = 0:(L-1)
        for fixX = 0:(L-1)
            for fixY = 0:(L-1)
                for loX = 0:(L-1)
                    for loY = 0:(L-1)
                        for hiX = 0:(L-1)
                            for hiY = 0:(L-1)
                                S1 = S0; 
                                %S1.eye_px_filt_trl = eyeX;
                                %S1.eye_py_filt_trl = eyeY;
                                S1.tgt_px_fx = fixX; 
                                S1.tgt_py_fx = fixY; 
                                S1.tgt_px_lo = loX; 
                                S1.tgt_py_lo = loY;
                                S1.tgt_px_hi = hiX; 
                                S1.tgt_py_hi = hiY;
                                ind = state2ind(S1,L); 
                                if ind > length(Qstate)
                                    error('State exceeds table size.');
                                end
                                Qstate{ind} = S1;
                            end
                        end
                    end
                end
            end
        end
%    end
%end

clear eyeX eyeY fixX fixY loX loY hiX hiY A1 S1

disp('State-Action space initialized.')

%% iterate until policy convergence 
ind3 = 1; Del = inf; theta = .001; 
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
    [~, Qtbl2, Srl] = ReinforcementLearnGrid(wT, gamma, .8, 200, L, Qstate, Qaction, Qtbl);
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