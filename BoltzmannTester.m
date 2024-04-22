%% load all data 
data_all_trials = Load_Analyzed_Data({'Y:\Ephys\data_132F\2024-02\', 'Y:\Ephys\data_132F\2024-03\', 'Y:\Ephys\data_132F\2024-04\'});

%% test random 
ind1 = randi(length(data_all_trials));
dtas = data_all_trials{ind1};
ind2 = randi(length(dtas)); 
dta = dtas(ind2);

Boltzmann(dta, true);

%% test all 

NT = 0;
for ind1 = 1:length(data_all_trials)
    dta = data_all_trials{ind1};
    NT = NT + length(dta);
end
Ni = 1;

boltz1 = zeros(NT,1); boltz2 = zeros(NT,1); boltz3 = zeros(NT,1);
L1 = zeros(NT,1); L2 = zeros(NT,1); L3 = zeros(NT,1);
P1 = zeros(NT,1); P2 = zeros(NT,1); P3 = zeros(NT,1);
h = zeros(NT,2); n = zeros(NT,2);
newsess = false(NT,1); newsess(1) = true;
alpha0 = []; beta0 = [];

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

        [boltz_rat, boltz_bay1, boltz_bay2, ...
          L_rat, L_bay1, L_bay2, ...
          mean_acc_rat, mean_acc_bay1, mean_acc_bay2, ...
          alpha_end, beta_end, denom_end, h_end, n_end] = ...
        Boltzmann(dta, false, alpha0, beta0);

        boltz1(Ni) = boltz_rat;
        boltz2(Ni) = boltz_bay1;
        boltz3(Ni) = boltz_bay2;

        L1(Ni) = L_rat;
        L2(Ni) = L_bay1;
        L3(Ni) = L_bay2;

        P1(Ni) = mean_acc_rat; 
        P2(Ni) = mean_acc_bay1; 
        P3(Ni) = mean_acc_bay2; 

        h(Ni,:) = h_end; n(Ni,:) = n_end; 

        Ni = Ni+1;
    end
    newsess(Ni-1) = true;
end

pRew = h./n;

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
subplot(3,1,1); 
plot(P3); grid on; 
title('Mean Accuracy'); 
ylabel('probability of trial'); xlabel('rec'); 
xticks(sessInd);
subplot(3,1,2); 
plot(pRew); grid on; hold on; 
%plot(find(newsess), pRew(newsess,:), '*'); 
ylim([0,1]);
title('Reward Likelihood'); 
ylabel('probability of high'); xlabel('rec'); 
xticks(sessInd);
legend('Cue Low', 'Cue High', 'Location','best')
subplot(3,1,3); 
plot(boltz3); grid on; 
title('Boltzmann Rationality'); 
ylabel('\beta'); xlabel('rec');
xticks(sessInd);