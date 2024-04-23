%% test random 
ind1 = randi(length(data_all_trials));
dtas = data_all_trials{ind1};
ind2 = randi(length(dtas)); 
dta = dtas(ind2);
trl = randi(length(dta.task_cond));
%trl = trl+1;

getStateSpace(dta, trl, true);