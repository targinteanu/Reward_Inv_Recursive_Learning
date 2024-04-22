%% test random 
ind1 = randi(length(data_all_trials));
dtas = data_all_trials{ind1};
ind2 = randi(length(dtas)); 
dta = dtas(ind2);

figure('Units', 'normalized', 'Position', [.1,.1,.5,.5]); plot(dta.tgt_px, dta.tgt_py, '-x')
hold on; plot(dta.eye_px_filt, dta.eye_py_filt);