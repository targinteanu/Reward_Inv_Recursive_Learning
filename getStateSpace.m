function SS = getStateSpace(dta, trl, depict)

if nargin < 2
    depict = false;
end

t = dta.time{trl};
eye_px_filt_trl = dta.eye_px_filt{trl};
eye_py_filt_trl = dta.eye_py_filt{trl};
tgt_px_trl = dta.tgt_px{trl};
tgt_py_trl = dta.tgt_py{trl};

movethresh = .25;
newState = (abs(diff(tgt_px_trl)) >= movethresh) | (abs(diff(tgt_py_trl)) >= movethresh);
newObs = (abs(diff(eye_px_filt_trl)) >= movethresh) | (abs(diff(eye_py_filt_trl)) >= movethresh);
newState = [false; newState]; newObs = [false; newObs]; % ignore first position
sel = newState | newObs;

SS = timetable(seconds(t(sel)), ...
    eye_px_filt_trl(sel), eye_py_filt_trl(sel), ...
    tgt_px_trl(sel), tgt_py_trl(sel));

if depict
    taskcondtype = {'choice', 'forced'};
    figure('Units', 'normalized', 'Position', [.1,.1,.5,.5]);
    %%{
    plot(tgt_px_trl, tgt_py_trl, 'b-o')
    hold on;
    plot(tgt_px_trl(1), tgt_py_trl(1), 'bx', 'MarkerSize',10);
    plot(tgt_px_trl(end), tgt_py_trl(end), 'b+', 'MarkerSize',10);
    %}
    plot(eye_px_filt_trl, eye_py_filt_trl, 'r');
    hold on;
    plot(eye_px_filt_trl(1), eye_py_filt_trl(1), 'rd', 'MarkerSize',8);
    plot(eye_px_filt_trl(end), eye_py_filt_trl(end), 'rs', 'MarkerSize',8);
    plot(dta.cue_x_high_rew(trl), dta.cue_y_high_rew(trl), '^k', 'MarkerSize',10,'LineWidth',3);
    plot(dta.cue_x_low_rew(trl),  dta.cue_y_low_rew(trl),  'vk', 'MarkerSize',10,'LineWidth',3);
    plot(dta.cue_x(trl),  dta.cue_y(trl),  'ok', 'MarkerSize',10,'LineWidth',2);
    plot(dta.start_x(trl), dta.start_y(trl), 'xk', 'MarkerSize',8,'LineWidth',2);
    plot(dta.end_x(trl), dta.end_y(trl), '+k', 'MarkerSize',8,'LineWidth',2);
    plot(eye_px_filt_trl(newState), eye_py_filt_trl(newState), 'or', 'MarkerSize',10);
    grid on;
    legend(...
        'tgt p', 'tgt start', 'tgt end', ...
        'eye p filt', 'eye start', 'eye end', ...
        'cue hi', 'cue lo', 'cue', 'start', 'end', ...
        'cue moved', ...
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

    figure('Units', 'normalized', 'Position', [.5,.5,.4,.4]);
    plot(t, tgt_px_trl, '-om'); hold on; plot(t, tgt_py_trl, '-og');
    plot(t, eye_px_filt_trl, 'm'); plot(t, eye_py_filt_trl, 'g');
    grid on;
    legend('tgt x', 'tgt y', 'eye x', 'eye y', 'Location','eastoutside');
    xlabel('time (s)'); ylabel('pos');
    title(ttl);
end

end