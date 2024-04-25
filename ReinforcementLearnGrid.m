function [Qfun, Qtable] = ReinforcementLearnGrid(R, tf)
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
S0 = gridifyState(S0);

Aspace = repmat({A0}, 1, L^2); % all possible actions 
c = 1;
for eyeX = 0:(L-1)
    for eyeY = 0:(L-1)
        A1 = A0; 
        A1.eye_px_filt_trl = eyeX; 
        A1.eye_py_filt_trl = eyeY; 
        Aspace{c} = A1;
        c = c+1;
    end
end

% start with a random state. Pick this by taking a random action from the
% default state. 
S0 = updateGridState(S0, Aspace{randi(L^2)});

% Q table will be filled when called 
Qtable = [nan, Aspace];
Qtable(2,1) = S0;

%% iteration 
for t = 0:tf

end

%% helper
    function q = getQ(s,a)
        [row,col] = findQ(s,a); 
        q = Qtable{row,col}; 
        if isempty(q)
            q = rand; 
            Qtable{row,col} = q;
        end
    end
    function setQ(s,a,q)
        [row,col] = fundQ(s,a);
        Qtable{row,col} = q;
    end

    function [r,c] = findQ(s,a)
        found = false;
        r = 2; 
        while ~found & (r <= size(Qtable,1))
            if isEqualState(s, Qtable{r,1})
                % found state (row) 
                c = findQc(a)
            else
                r = r+1;
            end
        end
        if ~found
            % state (row) not yet observed 
            Qtable{r,1} = s;
            c = findQc(a); 
        end
    end

    function c = findQc(a)
        found = false;
        c = 2; 
        while ~found & (c <= size(Qtable,2))
            if isEqualState(a, Qtable{1,c})
                % found action (col) 
                found = true;
            else
                c = c+1;
            end
        end
        if ~found
            % found row but not column 
            warning('unexpected action')
            Qtable{1,c} = a;
            % found = true;
        end
    end

end