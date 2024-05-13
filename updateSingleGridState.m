function s2 = updateGridState(s1, a1, posOpts, L)
% updateState, but both input state and action and output state are in grid world.  
% there may be a better way to do this

%s2 = gridifyState(updateState(ungridState(s1), ungridState(a1)));

if nargin < 3
    [~,posOpts] = gridifyState();
    posOpts = [-inf,posOpts,inf];
end
if nargin < 4
    L = length(posOpts);
end

tgt_fx = [s1.tgt_px_fx, tgt_py_fx]; % in coord sys of s1
tgt_lo = [s1.tgt_px_lo, tgt_py_lo]; % in coord sys of s1
tgt_hi = [s1.tgt_px_hi, tgt_py_hi]; % in coord sys of s1

s2_eye = a1{:,:}; % in coord sys of s1

[~,i0] = min(abs(posOpts)); % origin 
[~,ifx] = min(abs(posOpts-.5)); % fixation x home 
[~,ify] = min(abs(posOpts+1)); % fixation y home
% center to new coord sys 


% init s2
s2 = s1; 
% zero eye coords 
s2.eye_px_filt_trl = i0-1;
s2.eye_py_filt_trl = i0-1;

if isequal(s2_eye, tgt_fx)
    % moved to fixation target
    % expected: fixation goes away; low and high appear at origin
    s2.tgt_px_fx = L-1; s2.tgt_py_fx = L-1; 
    s2.tgt_px_lo = i0-1; s2.tgt_py_lo = i0-1;
    s2.tgt_px_hi = i0-1; s2.tgt_py_hi = i0-1;

elseif isequal(s2_eye, tgt_lo)
    % moved to low reward target 
    % expected: low goes away; fixation appears at [.5, -1]; high stays same
    s2.tgt_px_lo = L-1; s2.tgt_py_lo = L-1;
    s2.tgt_px_fx = ifx-1; s2.tgt_py_fx = ify-1;

elseif isequal(s2_eye, tgt_hi)
    % moved to high reward target 
    % expected: high goes away; fixation appears at [.5, -1]; low stays same
    s2.tgt_px_hi = L-1; s2.tgt_py_hi = L-1;
    s2.tgt_px_fx = ifx-1; s2.tgt_py_fx = ify-1;

else
    % no movement to any target
    % expected: enviro stays the same (do nothing else)

end

end