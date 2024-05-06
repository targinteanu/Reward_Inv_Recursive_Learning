function [mean_acc_rat, mean_acc_bay1, mean_acc_bay2, ...
          alpha_end, beta_end, denom_end, h_end, n_end] = ...
    giveBoltzmann(dta, boltz, depict, alpha0, beta0)
% stateless Boltzmann decision-making 
% Determine Boltzmann rationality during a single rec data (dta). 
% If depict is true, details will be shown. 
% alpha0 and beta0 define the initial prior pdf, with h ~ alpha0 
%   and n ~ alpha0 + beta0 (i.e. MLE = alpha0/(alpha0+beta0)

boltz_rat = boltz(1); 
boltz_bay1 = boltz(2); 
boltz_bay2 = boltz(3);

if nargin < 5
    beta0 = [];
    if nargin < 4
        alpha0 = [];
        if nargin < 3
            depict = false;
        end
    end
end

if isempty(alpha0)
    if ~isempty(beta0)
        alpha0 = beta0;
    else
        alpha0 = 1; beta0 = 1;
    end
else
    if isempty(beta0)
        beta0 = alpha0;
    end
end

if numel(beta0) == 1
    beta0 = [beta0; beta0];
end
if numel(alpha0) == 1
    alpha0 = [alpha0; alpha0];
end

%% determine rational belief function

Rcue = dta.tgt_cond; 
Rtru = dta.rew_cond; 
a_all = dta.choice; 

choicetrl = dta.task_cond == 0;
forcedtrl = dta.task_cond == 1;

% belief from ratio only 
Brat = nan(2,length(Rtru)+1); % belief that (Rtru | [Rcue=0; Rcue=1]) = 1 
Brat0 = alpha0./(alpha0 + beta0); % prior 
Brat(:,1) = Brat0; 

% bayesian belief using beta/binomial distribution 
% B's dim 3 params: [alpha, beta, denom] 
% B(r|obs) = (1/denom) * r^alpha * (1-r)^beta
% r_ML = alpha/(alpha+beta)
b = @(theta,alpha,beta) theta.^alpha .* (1-theta).^beta; % unscaled beta distribution
Bbay = nan(2,length(Rtru)+1,3);  
denom0 = [integral(@(r) b(r,alpha0(1),beta0(1)), 0, 1); ...
          integral(@(r) b(r,alpha0(2),beta0(2)), 0, 1)];
Bbay0 = cat(3, alpha0, beta0, denom0); % prior 
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

alpha_end = Bbay(:,end,1); 
beta_end  = Bbay(:,end,2);
denom_end = Bbay(:,end,3);

cueL = Rcue==0; cueH = Rcue==1; 
n_end = [sum(cueL), sum(cueH)];
h_end = [sum(Rtru(cueL) == 1), sum(Rtru(cueH) == 1)];

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

prb = @(a, boltz, R) exp(boltz.*R(a,:)) ./ (exp(boltz.*R(1,:)) + exp(boltz.*R(2,:))); 

PRB = nan(3,size(Rbay2,2));
for t = 1:size(PRB,2)
    PRB(1,t) = prb(aind(t), boltz_rat, Rrat(:,t));
    PRB(2,t) = prb(aind(t), boltz_bay1, Rbay1(:,t));
    PRB(3,t) = prb(aind(t), boltz_bay2, Rbay2(:,t));
end
meanAcc = mean(PRB,2);
mean_acc_rat = meanAcc(1); 
mean_acc_bay1 = meanAcc(2);
mean_acc_bay2 = meanAcc(3);

if depict
    figure('Units', 'normalized', 'Position', [.4,.1,.5,.3]); 
    plot(PRB'); grid on; 
    title('Each Trial Likelihood'); 
    ylabel('probability'); xlabel('trial'); 
    legend('ratio', 'bayes max like', 'bayes expected', 'Location','westoutside')
end

%% helpers
    function x = fzero_wrapper(fun, x0)
        try 
            x = fzero(fun, x0);
        catch ME
            warning(ME.message)
            x = nan; 
        end
    end

end