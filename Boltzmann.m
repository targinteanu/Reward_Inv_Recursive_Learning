%% load data to array 

f_name_end = 'ANALYZED.mat';
f_name_str = 'reward_prediction';

basedirs = {'Y:\Ephys\data_132F\2024-03\', 'Y:\Ephys\data_132F\2024-04\'}; 

data_all_trials = {};

for bd = basedirs
    data_this_dir = {};
    basedir = bd{:};
    sessions = dir(basedir);
    sessions = sessions(3:end); % remove '.', '..' 
    for session = sessions'
        if ~isempty(session) & session.isdir
            data_this_session = [];
            recs = dir([session.folder,filesep,session.name]);
            recs = recs(3:end); % remove '.', '..'
            for rec = recs'
                if ~isempty(rec) & rec.isdir
                    filepath = [rec.folder,filesep,rec.name,filesep,'analyzed_data',filesep,'behavior_data',filesep,'eye'];
                    if isfolder(filepath)
                        files = dir(filepath);
                        for f = files'
                            if ~isempty(f) & ~f.isdir
                                if strcmpi(f_name_end, f.name((end-length(f_name_end)+1):end) ) & ...
                                   strcmpi(f_name_str, f.name(1:length(f_name_str)) )
                                    load([f.folder,filesep,f.name], 'trials_data');
                                    tr_d = struct('task_cond', trials_data.task_cond, ...
                                                  'tgt_cond', trials_data.tgt_cond, ...
                                                  'rew_cond', trials_data.rew_cond, ...
                                                  'jump_cond', trials_data.jump_cond, ...
                                                  'choice', trials_data.choice, ...
                                                  'cue_x', trials_data.cue_x, ...
                                                  'cue_y', trials_data.cue_y, ...
                                                  'tgt_px', trials_data.tgt_px, ...
                                                  'tgt_py', trials_data.tgt_py, ...
                                                  'eye_px_filt', trials_data.eye_px_filt, ...
                                                  'eye_py_filt', trials_data.eye_py_filt, ...
                                                  'eye_vx_filt', trials_data.eye_vx_filt, ...
                                                  'eye_vy_filt', trials_data.eye_vy_filt ...
                                                  );
                                    data_this_session = [data_this_session, tr_d];
                                    clear trials_data sac_data meta_data tr_d
                                end
                            end
                        end
                    end
                end
            end
            data_this_dir = [data_this_dir, data_this_session];
            clear data_this_session
        end
    end
    data_all_trials = [data_all_trials, data_this_dir];
    clear data_this_dir
end

%% stateless Boltzmann decision-making: determine rational belief function

ind1 = randi(length(data_all_trials));
dta = data_all_trials{ind1};
ind2 = randi(length(dta)); 
dta = dta(ind2);

Rcue = dta.tgt_cond; 
Rtru = dta.rew_cond; 
a_all = dta.choice; 

choicetrl = dta.task_cond == 0;
forcedtrl = dta.task_cond == 1;

% belief from ratio only 
Brat = nan(2,length(Rtru)+1); % belief that (Rtru | [Rcue=0; Rcue=1]) = 1 
Brat0 = [.5; .5]; % prior 
Brat(:,1) = Brat0; 

% bayesian belief using beta/binomial distribution 
% B's dim 3 params: [alpha, beta, denom] 
% B(r|obs) = (1/denom) * r^alpha * (1-r)^beta
% r_ML = alpha/(alpha+beta)
b = @(theta,alpha,beta) theta.^alpha .* (1-theta).^beta; % unscaled beta distribution
Bbay = nan(2,length(Rtru)+1,3);  
alpha0 = 1; beta0 = 1; 
denom0 = integral(@(r) b(r,alpha0,beta0), 0, 1);
Bbay0 = cat(3, [alpha0;alpha0], [beta0;beta0], [denom0;denom0]); % prior 
Bbay(:,1,:) = Bbay0;

for t = 1:length(Rtru)
    % Brat_prev = Brat(:,t); Bbay_prev = Bbay(:,t); % prior 
    if dta.task_cond(t)
        % forced trial 
        cueChosen = Rcue(t);
    else
        % choice trial
        cueChosen = a_all(t);
    end
    pastDataInd = find(Rcue(1:t) == cueChosen); 
    n = length(pastDataInd); % # of trials 
    h = sum(Rtru(pastDataInd)); % # of hi reward (r=1)
    cueChosen = cueChosen + 1; % row index of B: 1 or 2

    % cue not chosen 
    cueNotChosen = 3-cueChosen; 
    Brat(cueNotChosen,t+1) = Brat(cueNotChosen,t);
    Bbay(cueNotChosen,t+1,:) = Bbay(cueNotChosen,t,:);

    % ratio 
    Brat(cueChosen,t+1) = h/n;

    % bayesian 
    alpha0 = Bbay(cueChosen,1,1);
    beta0  = Bbay(cueChosen,1,2);
    denom0 = Bbay(cueChosen,1,3); 
    alpha = h + alpha0; beta = n - h + beta0; 
    %denom = denom0 * factorial(h)*factorial(n-h)/factorial(n); 
    denom = integral(@(r) b(r,alpha,beta), 0, 1);
    Bbay(cueChosen,t+1,1) = alpha; 
    Bbay(cueChosen,t+1,2) = beta;
    Bbay(cueChosen,t+1,3) = denom; 
end

%% visualize rational belief function 
%{
numTrl = size(Brat,2); numToView = 8;
trlToView = sort(randperm(numTrl,numToView)); 
%trlToView = 1:numToView;
lgd = reshape([trlToView; trlToView], 1,2*numToView); 
lgd = arrayfun(@(n) ['Trial ',num2str(n)], lgd, 'UniformOutput',false);
figure('Units', 'normalized', 'Position', [.1,.1,.5,.8]); 
rplot = linspace(0,1,1000);
for t = trlToView
    alpha = Bbay(:,t,1); beta = Bbay(:,t,2); denom = Bbay(:,t,3); 
    rat = Brat(:,t);
    clr = colorwheel(t/numTrl); % color 
    lwd = 3-2.5*t/numTrl; % line width
    for q = 1:2
        Bbayplot = b(rplot,alpha(q),beta(q))/denom(q);
        Bratplot = rat(q);
        subplot(2,1,q);
        plot(rplot,Bbayplot, 'Color',clr, 'LineWidth',lwd); 
        hold on; 
        plot([Bratplot,Bratplot],[0,max(Bbayplot)], 'Color',clr, 'LineWidth',lwd);
    end
end
subplot(2,1,1); title('Cue 0 rational belief'); 
xlabel('reward'); ylabel('probability density'); 
legend(lgd, 'Location','eastoutside'); grid on; 
subplot(2,1,2); title('Cue 1 rational belief'); 
xlabel('reward'); ylabel('probability density'); 
legend(lgd, 'Location','eastoutside'); grid on; 

figure('Units', 'normalized', 'Position', [.1,.1,.5,.8]); 
for q = 1:2
    BbayImg = nan(length(rplot),numTrl);
    for t = 1:numTrl
        alpha = Bbay(:,t,1); beta = Bbay(:,t,2); denom = Bbay(:,t,3);
        Bbayplot = b(rplot,alpha(q),beta(q))/denom(q);
        BbayImg(:,t) = Bbayplot;
    end
    subplot(2,1,q);
    imagesc(1:numTrl, rplot, BbayImg); colorbar; 
end
subplot(2,1,1); title('Cue 0 rational belief'); 
ylabel('reward'); xlabel('trial');
subplot(2,1,2); title('Cue 1 rational belief'); 
ylabel('reward'); xlabel('trial');
%}

%% find max likelihood Boltzmann rationality 

a = a_all(choicetrl); aind = a+1; % index of choice 
Rrat = Brat(:,[choicetrl,false],:);
Rbay = Bbay(:,[choicetrl,false],:);
Rbay1 = Rbay(:,:,1)./(Rbay(:,:,1) + Rbay(:,:,2)); % ML reward 
Rbay2 = nan(size(Rbay1)); % expected reward
for t = 1:size(Rbay2,2)
    for q = 1:size(Rbay2,1)
        alpha = Rbay(q,t,1); 
        beta  = Rbay(q,t,2);
        denom = Rbay(q,t,3);
        Rbay2(q,t) = integral(@(r) r.*b(r,alpha,beta), 0,1)/denom;
    end
end

LHSrat = sum(diag(Rrat(aind,:))); 
LHSbay1 = sum(diag(Rbay1(aind,:)));
LHSbay2 = sum(diag(Rbay2(aind,:)));

dLL_RHS = @(boltz, R1, R2) sum( (R1.*exp(boltz.*R1) + R2.*exp(boltz.*R2)) ...
                               ./ (exp(boltz.*R1) + exp(boltz.*R2)) , 2);
LL_RHS = @(boltz, R1, R2) sum( log(exp(boltz.*R1) + exp(boltz.*R2)) , 2);
prb = @(a, boltz, R) exp(boltz.*R(a,:)) ./ (exp(boltz.*R(1,:)) + exp(boltz.*R(2,:))); 

bmax = 50;
bvals = linspace(-bmax,bmax,1000)';
figure('Units', 'normalized', 'Position', [.4,.1,.5,.8]); 

subplot(2,1,1);
LL = exp( bvals*LHSrat - LL_RHS(bvals, Rrat(1,:), Rrat(2,:)) );
plot(bvals, LL); grid on; hold on; 
LL = exp( bvals*LHSbay1 - LL_RHS(bvals, Rbay1(1,:), Rbay1(2,:)) );
plot(bvals, LL); 
LL = exp( bvals*LHSbay2 - LL_RHS(bvals, Rbay2(1,:), Rbay2(2,:)) );
plot(bvals, LL);
title('likelihood'); xlabel('Boltzmann rationality');
%ylabel('probability');
legend('ratio', 'bayes max like', 'bayes expected', 'Location','westoutside')

subplot(2,1,2);
dLL = dLL_RHS(bvals, Rrat(1,:), Rrat(2,:)) - LHSrat;
plot(bvals, dLL); grid on; hold on; 
dLL = dLL_RHS(bvals, Rbay1(1,:), Rbay1(2,:)) - LHSbay1;
plot(bvals, dLL); 
dLL = dLL_RHS(bvals, Rbay2(1,:), Rbay2(2,:)) - LHSbay2;
plot(bvals, dLL);
title('derivative log likelihood'); xlabel('Boltzmann rationality');
legend('ratio', 'bayes max like', 'bayes expected', 'Location','westoutside')

boltz_rat = fzero(@(boltz) dLL_RHS(boltz, Rrat(1,:), Rrat(2,:)) - LHSrat, [-bmax,bmax])
L_rat = exp( boltz_rat*LHSrat - LL_RHS(boltz_rat, Rrat(1,:), Rrat(2,:)) )
boltz_bay1 = fzero(@(boltz) dLL_RHS(boltz, Rbay1(1,:), Rbay1(2,:)) - LHSbay1, [-bmax,bmax])
L_bay1 = exp( boltz_bay1*LHSbay1 - LL_RHS(boltz_bay1, Rbay1(1,:), Rbay1(2,:)) )
boltz_bay2 = fzero(@(boltz) dLL_RHS(boltz, Rbay2(1,:), Rbay2(2,:)) - LHSbay2, [-bmax,bmax])
L_bay2 = exp( boltz_bay2*LHSbay2 - LL_RHS(boltz_bay2, Rbay2(1,:), Rbay2(2,:)) )

figure('Units', 'normalized', 'Position', [.4,.1,.5,.3]); 
PRB = nan(3,size(Rbay2,2)); 
for t = 1:size(PRB,2)
    PRB(1,t) = prb(aind(t), boltz_rat, Rrat(:,t)); 
    PRB(2,t) = prb(aind(t), boltz_bay1, Rbay1(:,t)); 
    PRB(3,t) = prb(aind(t), boltz_bay2, Rbay2(:,t)); 
end
plot(PRB'); grid on; 
title('Each Trial Likelihood'); 
ylabel('probability'); xlabel('trial'); 
legend('ratio', 'bayes max like', 'bayes expected', 'Location','westoutside')