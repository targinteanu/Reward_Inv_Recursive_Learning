%% test random 
ind1 = randi(length(data_all_trials));
dtas = data_all_trials{ind1};
ind2 = randi(length(dtas)); 
dta = dtas(ind2);
%%
%trl = randi(length(dta.task_cond));
trl = trl+1;

taskcondtype = {'choice', 'forced'}; 
figure('Units', 'normalized', 'Position', [.1,.1,.5,.5]); 
%%{
plot(dta.tgt_px{trl}, dta.tgt_py{trl}, 'b-o')
hold on; 
plot(dta.tgt_px{trl}(1), dta.tgt_py{trl}(1), 'bx', 'MarkerSize',10);
plot(dta.tgt_px{trl}(end), dta.tgt_py{trl}(end), 'b+', 'MarkerSize',10);
%}
plot(dta.eye_px_filt{trl}, dta.eye_py_filt{trl}, 'r');
hold on;
plot(dta.eye_px_filt{trl}(1), dta.eye_py_filt{trl}(1), 'rd', 'MarkerSize',8);
plot(dta.eye_px_filt{trl}(end), dta.eye_py_filt{trl}(end), 'rs', 'MarkerSize',8);
plot(dta.cue_x_high_rew(trl), dta.cue_y_high_rew(trl), '^k', 'MarkerSize',10,'LineWidth',3);
plot(dta.cue_x_low_rew(trl),  dta.cue_y_low_rew(trl),  'vk', 'MarkerSize',10,'LineWidth',3);
plot(dta.cue_x(trl),  dta.cue_y(trl),  'ok', 'MarkerSize',10,'LineWidth',2);
plot(dta.start_x(trl), dta.start_y(trl), 'xk', 'MarkerSize',8,'LineWidth',2);
plot(dta.end_x(trl), dta.end_y(trl), '+k', 'MarkerSize',8,'LineWidth',2);
grid on; 
legend(...
    'tgt p', 'tgt start', 'tgt end', ...
    'eye p filt', 'eye start', 'eye end', ...
    'cue hi', 'cue lo', 'cue', 'start', 'end', ...
    'Location', 'eastoutside');
ttl = [taskcondtype{dta.task_cond(trl)+1},' trial; ']; 
if dta.task_cond(trl)
    ttl = [ttl,' target cond = ',num2str(dta.tgt_cond(trl)),'; '];
else
    ttl = [ttl,' choice = ',num2str(dta.choice(trl)),'; '];
end
if dta.jump_cond(trl)
    ttl = [ttl,' with jump'];
else
    ttl = [ttl,' no jump'];
end
title(ttl);

figure('Units', 'normalized', 'Position', [.4,.4,.5,.5]); 
plot(dta.tgt_px{trl}, '-o'); hold on; plot(dta.tgt_py{trl}, '-o'); 
plot(dta.eye_px_filt{trl}); plot(dta.eye_py_filt{trl});
grid on; 
xlabel('time (sample)'); ylabel('pos');
title(ttl);