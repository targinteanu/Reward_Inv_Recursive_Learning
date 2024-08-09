function [boltz, L, meanAcc, Q] = getQ(dta, depict, Qlo0, Qhi0, lrnrt, frgrt)
% stateless Boltzmann decision-making 
% Determine Boltzmann rationality during a single rec data (dta). 
% If depict is true, details will be shown. 
% alpha0 and beta0 define the initial prior pdf, with h ~ alpha0 
%   and n ~ alpha0 + beta0 (i.e. MLE = alpha0/(alpha0+beta0)
% try a few different methods for the expected reward: ratio h/n (rat),
% bayesian maximum likelihood (bay1), bayesian expected value (bay2). 

if nargin < 6
    frgrt = [];
    if nargin < 5
        lrnrt = [];
        if nargin < 4
            Qhi0 = [];
            if nargin < 3
                Qlo0 = [];
                if nargin < 2
                    depict = false;
                end
            end
        end
    end
end

if isempty(Qlo0)
    if ~isempty(Qhi0)
        Qlo0 = Qhi0;
    else
        Qlo0 = 0; Qhi0 = 0;
    end
else
    if isempty(Qhi0)
        Qhi0 = Qlo0;
    end
end

if isempty(lrnrt)
    lrnrt = 1;
end
if isempty(frgrt)
    frgrt = 1;
end

%% determine rational belief function

Rcue = dta.tgt_cond; 
Rtru = dta.rew_cond; 
a_all = dta.choice; 

choicetrl = dta.task_cond == 0;
forcedtrl = dta.task_cond == 1;

Q = zeros(2, (length(Rcue)+1)); 
Q(:,1) = [Qlo0; Qhi0];

for t = 1:length(Rtru)
    if dta.task_cond(t)
        % forced trial 
        cueChosen = Rcue(t);
    else
        % choice trial
        cueChosen = a_all(t);
    end

    % cue chosen 
    cueChosen = cueChosen + 1; % row index of Q: 1 or 2
    Q(cueChosen,t+1) = frgrt*Q(cueChosen,t) + lrnrt*(Rtru(t) - Q(cueChosen,t));

    % cue not chosen 
    cueNotChosen = 3-cueChosen; 
    Q(cueNotChosen,t+1) = frgrt*Q(cueNotChosen,t);
end

%% find max likelihood Boltzmann rationality 

a = a_all(choicetrl); aind = a+1; % index of choice 
Qchoice = Q(:,[choicetrl,false]);

if length(aind) > 1
    LHS = sum(diag(Qchoice(aind,:)));
else
    LHS = sum((Qchoice(aind,:)));
end

dLL_RHS = @(boltz, R1, R2) sum( (R1.*exp(boltz.*R1) + R2.*exp(boltz.*R2)) ...
                               ./ (exp(boltz.*R1) + exp(boltz.*R2)) , 2);
LL_RHS = @(boltz, R1, R2) sum( log(exp(boltz.*R1) + exp(boltz.*R2)) , 2);
prb = @(a, boltz, R) exp(boltz.*R(a,:)) ./ (exp(boltz.*R(1,:)) + exp(boltz.*R(2,:))); 

bmax = 50;

if depict

bvals = linspace(-bmax,bmax,1000)';
figure('Units', 'normalized', 'Position', [.4,.1,.5,.8]); 

subplot(2,1,1);
LL = exp( bvals*LHS - LL_RHS(bvals, Qchoice(1,:), Qchoice(2,:)) );
plot(bvals, LL); grid on; 
title('likelihood'); xlabel('Boltzmann rationality');
%ylabel('probability');

subplot(2,1,2);
dLL = dLL_RHS(bvals, Qchoice(1,:), Qchoice(2,:)) - LHS;
plot(bvals, dLL); grid on;  
title('derivative log likelihood'); xlabel('Boltzmann rationality');

end

boltz = fzero_wrapper(@(boltz) dLL_RHS(boltz, Qchoice(1,:), Qchoice(2,:)) - LHS, [-bmax,bmax]);
L = exp( boltz*LHS - LL_RHS(boltz, Qchoice(1,:), Qchoice(2,:)) );

PRB = nan(1,size(Qchoice,2));
for t = 1:size(PRB,2)
    PRB(1,t) = prb(aind(t), boltz, Qchoice(:,t));
end
meanAcc = mean(PRB,2);

if depict
    figure('Units', 'normalized', 'Position', [.4,.1,.5,.3]); 
    plot(PRB'); grid on; 
    title('Each Trial Likelihood'); 
    ylabel('probability'); xlabel('trial'); 
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