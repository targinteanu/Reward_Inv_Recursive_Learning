%{
function [S,A] = getStateSpace(dta, trl, depict)

if nargin < 2
    depict = false;
end
%}

t = dta.time{trl};
eye_px_filt_trl = dta.eye_px_filt{trl};
eye_py_filt_trl = dta.eye_py_filt{trl};
tgt_px_trl = dta.tgt_px{trl};
tgt_py_trl = dta.tgt_py{trl};
jump_dx = dta.end_x - dta.cue_x; 
jump_dy = dta.end_y - dta.cue_y;

%% construct new version of tgt_p
% in choice 1 trials with jump, tgt_p fails to follow the jump
% in forced trials, tgt_p inserts jump where there is none 
% forced 1 trials with jump OK
% choice no jump OK

tgt_px_fx = nan(size(tgt_px_trl)); tgt_py_fx = nan(size(tgt_py_trl)); 
tgt_px_lo = nan(size(tgt_px_trl)); tgt_py_lo = nan(size(tgt_py_trl)); 
tgt_px_hi = nan(size(tgt_px_trl)); tgt_py_hi = nan(size(tgt_py_trl));

% only fixation target visible until fixation 
tgtFix = [dta.start_x(trl), dta.start_y(trl)];
reward_radius = dta.reward_area(trl);
reward_radius = sqrt(reward_radius/pi);
[tgt_px_fx, tgt_py_fx, eye_p1, tInd1] = ...
    unicycle(tgtFix, inf, 0, eye_px_filt_trl, eye_py_filt_trl, tgt_px_fx, tgt_py_fx, reward_radius, length(t));

tgtHi = [dta.cue_x_high_rew(trl), dta.cue_y_high_rew(trl)];
tgtLo = [dta.cue_x_low_rew(trl), dta.cue_y_low_rew(trl)]; 

if dta.task_cond(trl)
    % forced 

    if dta.tgt_cond(trl)
        % forced hi 
        tgt2 = tgtHi; 
        if dta.jump_cond(trl)
            tgt3 = tgt2 + [jump_dx(trl), jump_dy(trl)]; 
        else
            tgt3 = tgt2;
        end
        [tgt_px_hi, tgt_py_hi, tgt_px_fx, tgt_py_fx] = tricycle(tgtFix,tgt2,tgt3, tInd1,eye_p1, ...
            eye_px_filt_trl, eye_py_filt_trl, tgt_px_hi, tgt_py_hi, tgt_px_fx, tgt_py_fx, reward_radius, length(t));

    else
        % forced lo 
        tgt2 = tgtLo; 
        if dta.jump_cond(trl)
            tgt3 = tgt2 + [jump_dx(trl), jump_dy(trl)]; 
        else
            tgt3 = tgt2;
        end
        [tgt_px_lo, tgt_py_lo, tgt_px_fx, tgt_py_fx] = tricycle(tgtFix,tgt2,tgt3, tInd1,eye_p1, ...
            eye_px_filt_trl, eye_py_filt_trl, tgt_px_lo, tgt_py_lo, tgt_px_fx, tgt_py_fx, reward_radius, length(t));
    end

else
    % choice 

    if dta.choice(trl)
        % chose hi 
        tgt2 = tgtHi; 
        if dta.jump_cond(trl)
            tgt3 = tgt2 + [jump_dx(trl), jump_dy(trl)]; 
        else
            tgt3 = tgt2;
        end
        [tgt_px_hi, tgt_py_hi, tgt_px_fx, tgt_py_fx] = tricycle(tgtFix,tgt2,tgt3, tInd1,eye_p1, ...
            eye_px_filt_trl, eye_py_filt_trl, tgt_px_hi, tgt_py_hi, tgt_px_fx, tgt_py_fx, reward_radius, length(t));
        [tgt_px_lo, tgt_py_lo, tgt_px_fx, tgt_py_fx] = tricycle(tgtFix,tgtLo,tgtLo, tInd1,eye_p1, ...
            eye_px_filt_trl, eye_py_filt_trl, tgt_px_lo, tgt_py_lo, tgt_px_fx, tgt_py_fx, reward_radius, length(t));
    
    else
        % chose lo
        tgt2 = tgtLo; 
        if dta.jump_cond(trl)
            tgt3 = tgt2 + [jump_dx(trl), jump_dy(trl)]; 
        else
            tgt3 = tgt2;
        end
        [tgt_px_lo, tgt_py_lo, tgt_px_fx, tgt_py_fx] = tricycle(tgtFix,tgt2,tgt3, tInd1,eye_p1, ...
            eye_px_filt_trl, eye_py_filt_trl, tgt_px_lo, tgt_py_lo, tgt_px_fx, tgt_py_fx, reward_radius, length(t));
        [tgt_px_hi, tgt_py_hi, tgt_px_fx, tgt_py_fx] = tricycle(tgtFix,tgtHi,tgtHi, tInd1,eye_p1, ...
            eye_px_filt_trl, eye_py_filt_trl, tgt_px_hi, tgt_py_hi, tgt_px_fx, tgt_py_fx, reward_radius, length(t));
    end
end

TGT = [tgt_px_fx, tgt_py_fx, ...
       tgt_px_hi, tgt_py_hi, ...
       tgt_px_lo, tgt_py_lo]; 
TGT(isnan(TGT)) = inf;

%%

movethresh = .25;
%newState = (abs(diff(tgt_px_trl)) >= movethresh) | (abs(diff(tgt_py_trl)) >= movethresh);
newState = abs(diff(TGT)) >= movethresh; newState = ~~sum(newState, 2);
newObs = (abs(diff(eye_px_filt_trl)) >= movethresh) | (abs(diff(eye_py_filt_trl)) >= movethresh);
newState = [false; newState]; newObs = [false; newObs]; % ignore first position
newObs(newState) = true; % new state is always observed, even if no action

S = timetable(seconds(t), ...
    eye_px_filt_trl, eye_py_filt_trl, ...
    tgt_px_fx, tgt_py_fx, ...
    tgt_px_lo, tgt_py_lo, ...
    tgt_px_hi, tgt_py_hi);
S = S(newObs,:);
% state is current eye and target positions 

A = [S.eye_px_filt_trl, S.eye_py_filt_trl]; % action uses eye position
A = A(2:end,:) - A(1:(end-1),:); % action is change in position 
S = S(1:(end-1),:); % last observation will be counted in the next trial

%%

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

%% functions 

function [tgtX, tgtY, eye_p1, tInd1] = unicycle(tgt1, dist2tgt, tInd0, eyeX, eyeY, tgtX, tgtY, reward_radius, T)
    tInd1 = tInd0; tInd0 = tInd0+1;
    eye_p1 = [nan, nan];
    while (tInd1 < T) & (dist2tgt > reward_radius)
        tInd1 = tInd1+1;
        eye_p1 = [eyeX(tInd1), eyeY(tInd1)];
        dist2tgt = norm(eye_p1 - tgt1);
    end
    tgtX(tInd0:tInd1) = tgt1(1); tgtY(tInd0:tInd1) = tgt1(2);
end


function [tgtX, tgtY, fixX, fixY] = tricycle(tgt1, tgt2, tgt3, tInd1, eye_p1, ...
    eyeX, eyeY, tgtX, tgtY, fixX, fixY, reward_radius, T)

    tInd2 = tInd1; tInd3 = tInd2;

    while (tInd2 < T) & (tInd1 < T) & (tInd3 < T)

        tInd2 = tInd1+1; 
        eye_p2 = [eyeX(tInd2), eyeY(tInd2)]; 
        while (tInd2 < T) & (norm(eye_p2 - tgt2) >= norm(eye_p1 - tgt2))
            % sac not yet started 
            tInd2 = tInd2+1;
            eye_p2 = [eyeX(tInd2), eyeY(tInd2)];
        end
        % sac started 
        tgtX((tInd1+1):tInd2) = tgt2(1); tgtY((tInd1+1):tInd2) = tgt2(2); 

        [tgtX, tgtY, eye_p2, tInd3] = ...
            unicycle(tgt3, norm(eye_p2 - tgt3), tInd2, eyeX, eyeY, tgtX, tgtY, reward_radius, T); 
        % reached target 

        [fixX, fixY, eye_p1, tInd1] = ... 
            unicycle(tgt1, norm(eye_p2 - tgt1), tInd3, eyeX, eyeY, fixX, fixY, reward_radius, T);
        % back to fixation 
    end
end

%end