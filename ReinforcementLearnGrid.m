function [Qfun, Qtable] = ReinforcementLearnGrid(R)
% use Q learning to determine the optimal value function Q given the reward
% function R(s)
% Return Q as a table with possible states as rows and actions as columns 
% and as a function of [state, action] 

%% initialize s and Q  
L = 39; % position takes L possible values from 0 to L-1

A0 = timetable(seconds(0), 0, 0); % null action 
A0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl'};
A0 = gridifyState(A0); 

S0 = timetable(seconds(0), 0,0, 0,0, inf,inf, inf,inf); % default state
S0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl', ...
    'tgt_px_fx', 'tgt_py_fx', 'tgt_px_lo', 'tgt_py_lo', 'tgt_px_hi', 'tgt_py_hi'};

Sspace = repmat({S0}, L^8, 1);
Aspace = repmat({A0}, 1, L^2);

r = 1;
for eyeX = 0:(L-1)
    for eyeY = 0:(L-1)
        for fixX = 0:(L-1)
            for fixY = 0:(L-1)
                for loX = 0:(L-1)
                    for loY = 0:(L-1)
                        for hiX = 0:(L-1)
                            for hiY = 0:(L-1)
                                S1 = S0;
                                S1.eye_px_filt_trl = eyeX; 
                                S1.eye_py_filt_trl = eyeY;
                                S1.tgt_px_fx = fixX; 
                                S1.tgt_py_fx = fixY; 
                                S1.tgt_px_lo = loX; 
                                S1.tgt_py_lo = loY; 
                                S1.tgt_px_hi = hiX; 
                                S1.tgt_py_hi = hiY;
                                Sspace{r} = S1;
                                r = r+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

end