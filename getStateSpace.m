function [S,A] = getStateSpace(dta, trl, depict)

if nargin < 2
    depict = false;
end

t = dta.time{trl};
eye_px_filt_trl = dta.eye_px_filt{trl};
eye_py_filt_trl = dta.eye_py_filt{trl};
tgt_px_trl = dta.tgt_px{trl};
tgt_py_trl = dta.tgt_py{trl};
jump_dx = dta.end_x - dta.cue_x; 
jump_dy = dta.end_y - dta.cue_y;

% in choice 1 trials with jump, tgt_p fails to follow the jump
% in forced trials, tgt_p inserts jump where there is none 
% forced 1 trials with jump OK
% choice no jump OK

movethresh = .25;
newState = (abs(diff(tgt_px_trl)) >= movethresh) | (abs(diff(tgt_py_trl)) >= movethresh);
newObs = (abs(diff(eye_px_filt_trl)) >= movethresh) | (abs(diff(eye_py_filt_trl)) >= movethresh);
newState = [false; newState]; newObs = [false; newObs]; % ignore first position
newObs(newState) = true; % new state is always observed, even if no action

tgt_px_lo = nan(size(tgt_px_trl)); tgt_py_lo = nan(size(tgt_py_trl)); 
tgt_px_hi = nan(size(tgt_px_trl)); tgt_py_hi = nan(size(tgt_py_trl));
if dta.task_cond(trl)
    % forced 
    if dta.tgt_cond(trl)
        % forced hi 
        tgt_px_hi = tgt_px_trl; 
        tgt_py_hi = tgt_py_trl;
    else
        % forced lo 
        tgt_px_lo = tgt_px_trl; 
        tgt_py_lo = tgt_py_trl;
    end
else
    % choice 
    if dta.choice(trl)
        % chose hi 
        %tgt_px_hi = 0;
    else
        % chose lo
    end
end

S = timetable(seconds(t), ...
    eye_px_filt_trl, eye_py_filt_trl, ...
    tgt_px_lo, tgt_py_lo, ...
    tgt_px_hi, tgt_py_hi);
S = S(newObs,:);
% state is current eye and target positions 

A = [S.eye_px_filt_trl, S.eye_py_filt_trl]; % action uses eye position
A = A(2:end,:) - A(1:(end-1),:); % action is change in position 
S = S(1:(end-1),:); % last observation will be counted in the next trial

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