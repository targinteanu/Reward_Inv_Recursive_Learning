% loop through all 
chose_hi = [];
chose_lo = [];
force_hi = [];
force_lo = [];
chose_hi_jumped = [];
chose_lo_jumped = [];
force_hi_jumped = [];
force_lo_jumped = [];
dist_jumped = [];
dist_nojump = [];
dist = @(x1,y1,x2,y2) norm([x1-x2,y1-y2]);
for ind1 = 1:(length(data_all_trials))
    dtas = data_all_trials{ind1};
    for ind2 = 1:(length(dtas))
        dta = dtas(ind2);
        dist_from_prim = zeros(length(dta.task_cond),1); 
        for trl = 1:(length(dta.task_cond))
            dist_from_hi = dist(dta.eye_px_filt{trl}(end), dta.eye_py_filt{trl}(end), dta.cue_x_high_rew(trl), dta.cue_y_high_rew(trl));
            dist_from_lo = dist(dta.eye_px_filt{trl}(end), dta.eye_py_filt{trl}(end), dta.cue_x_low_rew(trl), dta.cue_y_low_rew(trl));
            dist_from_prim(trl) = min(dist_from_hi, dist_from_lo);
        end
            forced_hi = (dta.task_cond == 1) & (dta.tgt_cond == 1);
            forced_lo = (dta.task_cond == 1) & (dta.tgt_cond == 0);
            chosen_hi = (dta.task_cond == 0) & (dta.choice == 1);
            chosen_lo = (dta.task_cond == 0) & (dta.choice == 0);
            jumped = dta.jump_cond == 1;

            chose_hi = [chose_hi; mean(chosen_hi)]; 
            chose_lo = [chose_lo; mean(chosen_lo)];
            force_hi = [force_hi; mean(forced_hi)];
            force_lo = [force_lo; mean(forced_lo)]; 

            dist_jumped = [dist_jumped; mean(dist_from_prim(jumped))];
            dist_nojump = [dist_nojump; mean(dist_from_prim(~jumped))];

            chose_hi_jumped = [chose_hi_jumped; mean( jumped(chosen_hi) )];
            chose_lo_jumped = [chose_lo_jumped; mean( jumped(chosen_lo) )];
            force_hi_jumped = [force_hi_jumped; mean( jumped(forced_hi) )];
            force_lo_jumped = [force_lo_jumped; mean( jumped(forced_lo) )];
    end
end
%%
figure('Units', 'normalized', 'Position', [.1,.1,.5,.8]);
subplot(2,1,1); 
plot([dist_jumped, dist_nojump]); grid on;
legend('jump', 'no jump');
title('distance');

subplot(2,1,2);
plot(chose_hi_jumped, 'LineWidth',2.5); hold on; grid on;
plot(chose_lo_jumped, 'LineWidth',2);
plot(force_hi_jumped, 'LineWidth',1.5);
plot(force_lo_jumped, 'LineWidth',1);
legend('chose hi', 'chose lo', 'force hi', 'force lo');
title('jump likelihood');