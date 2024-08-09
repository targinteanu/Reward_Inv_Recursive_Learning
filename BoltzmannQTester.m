%% load all data 
%data_all_trials = Load_Analyzed_Data({'Y:\Ephys\data_132F\2024-02\', 'Y:\Ephys\data_132F\2024-03\', 'Y:\Ephys\data_132F\2024-04\'});
load("data_all_trials.mat")

%% test random 
ind1 = randi(length(data_all_trials));
dtas = data_all_trials{ind1};
ind2 = randi(length(dtas)); 
dta = dtas(ind2);

zeta = .7; alpha = .9;
getQ(dta, true, 0, 0, alpha, zeta);

%% test all 

NT = 0;
for ind1 = 1:length(data_all_trials)
    dta = data_all_trials{ind1};
    NT = NT + length(dta);
end
Ni = 1;
difftesttrain = 8;

boltz1 = zeros(NT,1); 
L1 = zeros(NT,1); 
P1 = zeros(NT,1); 
P4 = P1; 
newsess = false(NT,1); newsess(1) = true;
Qlo0 = []; Qhi0 = [];

for ind1 = 1:length(data_all_trials)
    disp(['Testing Series ',num2str(ind1),' of ',num2str(length(data_all_trials))])
    dtas = data_all_trials{ind1};
    progdisptick = .1; progdisp = progdisptick; 
    for ind2 = 1:length(dtas)
        if ind2/length(dtas) >= progdisp
            disp(['Testing Rec ',num2str(ind2),' of ',num2str(length(dtas))])
            progdisp = progdisp + progdisptick;
        end
        dta = dtas(ind2);

        [boltz, L, meanAcc, Q] = getQ(dta, false, Qlo0, Qhi0, alpha, zeta);

        if Ni <= difftesttrain
            b1 = boltz; 
        else
            b1 = boltz1(Ni-difftesttrain); 
        end

        boltz1(Ni) = boltz;
        L1(Ni) = L;
        P1(Ni) = meanAcc; 

        Qlo0 = Q(1,end); Qhi0 = Q(2,end);

        Ni = Ni+1;
    end
    newsess(Ni-1) = true;
end

P4(1:difftesttrain) = nan; P5(1:difftesttrain) = nan; P6(1:difftesttrain) = nan;
%% show results 

%{
figure('Units', 'normalized', 'Position', [.1,.1,.5,.8]); 
subplot(2,1,1); 
plot([L1, L2, L3]); grid on;
legend('ratio', 'bayes max like', 'bayes expected', 'Location','westoutside')
title('Likelihood'); xlabel('rec'); 
subplot(2,1,2); 
plot([boltz1, boltz2, boltz3]); grid on;
legend('ratio', 'bayes max like', 'bayes expected', 'Location','westoutside')
title('Boltzmann Rationality'); xlabel('rec'); 
%}
sessInd = find(newsess);
figure('Units', 'normalized', 'Position', [.1,.1,.5,.8]); 
subplot(2,1,1); 
plot(P1); grid on; 
title('Mean Accuracy'); 
ylabel('probability of trial'); xlabel('rec'); 
xticks(sessInd);
subplot(2,1,2); 
plot(boltz1); grid on; 
title('Boltzmann Rationality'); 
ylabel('\beta'); xlabel('rec');
xticks(sessInd);

%% display (v2) 
sesssel = ~isnan(boltz1); 
sessInd = find(newsess(sesssel));
figure('Units', 'normalized', 'Position', [.1,.1,.5,.5]); 

subplot(2,1,1); 
plot(P1(sesssel), 'LineWidth',1); grid on; 
title('Mean Accuracy of Proposed Probability Distribution'); 
ylabel('probability of trial'); xlabel('recording'); 
xticks(sessInd);

subplot(2,1,2); 
plot(boltz1(sesssel), 'LineWidth',1); grid on; 
title('Boltzmann Rationality Parameter'); 
ylabel('\beta'); xlabel('recording');
xticks(sessInd);

%{
subplot(2,1,1); hold on;
plot(P6(sesssel), 'LineWidth',1); grid on; 
legend({'Current \beta', 'Previous \beta'}, 'Location','best')
%}

%% stats 
x = boltz1; x = x(~isnan(x));
figure; histogram(x)
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM                       % Confidence Intervals