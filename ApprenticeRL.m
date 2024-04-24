% loop through all 
chose_hi_cue_hi = [];
chose_lo_cue_lo = [];
force_hi_cue_hi = [];
force_lo_cue_lo = [];
chose_hi = [];
chose_lo = [];
force_hi = [];
force_lo = [];
cue_hi = [];
cue_lo = [];
dist = @(x1,y1,x2,y2) norm([x1-x2,y1-y2]);
for ind1 = 1:(length(data_all_trials))
    dtas = data_all_trials{ind1};
    for ind2 = 1:(length(dtas))
        dta = dtas(ind2);
        goes_to_hi = false(length(dta.task_cond),1); 
        for trl = 1:(length(dta.task_cond))
            dist_from_hi = dist(dta.eye_px_filt{trl}(end), dta.eye_py_filt{trl}(end), dta.cue_x_high_rew(trl), dta.cue_y_high_rew(trl));
            dist_from_lo = dist(dta.eye_px_filt{trl}(end), dta.eye_py_filt{trl}(end), dta.cue_x_low_rew(trl), dta.cue_y_low_rew(trl));
            goes_to_hi(trl) = dist_from_hi < dist_from_lo;
        end
            forced_hi = (dta.task_cond == 1) & (dta.tgt_cond == 1);
            forced_lo = (dta.task_cond == 1) & (dta.tgt_cond == 0);
            chosen_hi = (dta.task_cond == 0) & (dta.choice == 1);
            chosen_lo = (dta.task_cond == 0) & (dta.choice == 0);

            chose_hi = [chose_hi; mean(chosen_hi)]; 
            chose_lo = [chose_lo; mean(chosen_lo)];
            force_hi = [force_hi; mean(forced_hi)];
            force_lo = [force_lo; mean(forced_lo)]; 
            cue_hi = [cue_hi; mean(goes_to_hi)];
            goes_to_lo = ~goes_to_hi;

            chose_hi_cue_hi = [chose_hi_cue_hi; mean( goes_to_hi(chosen_hi) )];
            chose_lo_cue_lo = [chose_lo_cue_lo; mean( goes_to_lo(chosen_lo) )];
            force_hi_cue_hi = [force_hi_cue_hi; mean( goes_to_hi(forced_hi) )];
            force_lo_cue_lo = [force_lo_cue_lo; mean( goes_to_lo(forced_lo) )];
    end
end
%%
figure('Units', 'normalized', 'Position', [.1,.1,.5,.8]);
subplot(2,1,1); 
plot([chose_hi, chose_lo, force_hi, force_lo]); grid on;
legend('chose hi', 'chose lo', 'force hi', 'force lo');
subplot(2,1,2);

plot(chose_hi_cue_hi, 'LineWidth',2.5); hold on; grid on;
plot(chose_lo_cue_lo, 'LineWidth',2);
plot(force_hi_cue_hi, 'LineWidth',1.5);
plot(force_lo_cue_lo, 'LineWidth',1);
legend('chose hi', 'chose lo', 'force hi', 'force lo');

%plot([cue_hi, cue_lo]); grid on; legend('cue hi', 'cue lo')