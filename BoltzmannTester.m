%% test random 
ind1 = randi(length(data_all_trials));
dta = data_all_trials{ind1};
ind2 = randi(length(dta)); 
dta = dta(ind2);

Boltzmann(dta, true, [1;2], [2;1]);

%% test all 

NT = 0;
for ind1 = 1:length(data_all_trials)
    dta = data_all_trials{ind1};
    NT = NT + length(dta);
end
Ni = 1;

boltz1 = zeros(NT,1); boltz2 = zeros(NT,1); boltz3 = zeros(NT,1);
L1 = zeros(NT,1); L2 = zeros(NT,1); L3 = zeros(NT,1);
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

        [boltz_rat, boltz_bay1, boltz_bay2, L_rat, L_bay1, L_bay2, ...
          alpha0, beta0] = ...
        Boltzmann(dta, false, alpha0, beta0);

        boltz1(Ni) = boltz_rat;
        boltz2(Ni) = boltz_bay1;
        boltz3(Ni) = boltz_bay2;

        L1(Ni) = L_rat;
        L2(Ni) = L_bay1;
        L3(Ni) = L_bay2;

        Ni = Ni+1;
    end
end

figure('Units', 'normalized', 'Position', [.1,.1,.5,.8]); 
subplot(2,1,1); 
plot([L1, L2, L3]); grid on;
legend('ratio', 'bayes max like', 'bayes expected', 'Location','westoutside')
title('Likelihood'); xlabel('rec'); 
subplot(2,1,2); 
plot([boltz1, boltz2, boltz3]); grid on;
legend('ratio', 'bayes max like', 'bayes expected', 'Location','westoutside')
title('Boltzmann Rationality'); xlabel('rec'); 