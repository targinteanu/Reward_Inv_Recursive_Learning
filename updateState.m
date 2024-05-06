function s2 = updateState(s1, a1)
% EXPECTED new state given current state and current action 

% mean, SD: 
% reward_area ~ 4.8228, 0.4951
% tgt_px_fx = .5
% tgt_py_fx = -1
% tgt_px_lo ~ 2.292, 3.8302
% tgt_py_lo ~ -1.7106, 4.3061
% tgt_px_hi ~ 1.6949, 3.8401
% tgt_py_hi ~ -1.0202, 3.9836

% actual: 
% tgt_hi and tgt_lo are at ~const radius but random angle, and there is
% some probability they won't appear at all; estimate does not account for
% this. 

% estimate: 
% fixation target is set; assume both other targets appear at zero (mean
% val is insignificant compared to variance); assume reward area of 4.8

reward_area = 4.8; 
reward_radius = sqrt(reward_area/pi);

%pos_new = [s1.eye_px_filt_trl, s1.eye_py_filt_trl] + [a1.eye_px_filt_trl, a1.eye_py_filt_trl];
pos_new = [a1.eye_px_filt_trl, a1.eye_py_filt_trl];
tgt_fx  = [s1.tgt_px_fx, s1.tgt_py_fx]; % fixation
tgt_lo  = [s1.tgt_px_lo, s1.tgt_py_lo]; % low reward
tgt_hi  = [s1.tgt_px_hi, s1.tgt_py_hi]; % high reward

s2 = s1;
s2.eye_px_filt_trl = pos_new(1); 
s2.eye_py_filt_trl = pos_new(2);

if norm(pos_new - tgt_fx) <= reward_radius
    % moved to fixation target
    % expected: fixation goes away; low and high appear at origin
    s2.tgt_px_fx = inf; s2.tgt_py_fx = inf;
    s2.tgt_px_lo = 0;   s2.tgt_py_lo = 0;
    s2.tgt_px_hi = 0;   s2.tgt_py_hi = 0;

elseif norm(pos_new - tgt_lo) <= reward_radius
    % moved to low reward target 
    % expected: low goes away; fixation appears at [.5, -1]; high stays same
    s2.tgt_px_lo = inf; s2.tgt_py_lo = inf; 
    s2.tgt_px_fx = .5;  s2.tgt_py_fx = -1;

elseif norm(pos_new - tgt_hi) <= reward_radius
    % moved to high reward target 
    % expected: high goes away; fixation appears at [.5, -1]; low stays same
    s2.tgt_px_hi = inf; s2.tgt_py_hi = inf; 
    s2.tgt_px_fx = .5;  s2.tgt_py_fx = -1;

else
    % no movement to any target
    % expected: enviro stays the same (do nothing else)

end

s2 = centerCoords(s2);

end