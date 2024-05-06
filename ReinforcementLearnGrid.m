function [Qfun, Qtable, Sall, Aall] = ReinforcementLearnGrid(R, gamma, alpha, tf)
% use Q learning to determine the optimal value function Q given the reward
% function R(s)
% Return Q as a table with possible states as rows and actions as columns 
% and as a function of [state, action] 

%% initialize s and Q  
[~,posOpts] = gridifyState();
L = length(posOpts); % position takes L possible values from 0 to L-1

A0 = timetable(seconds(0), 0, 0); % null action 
A0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl'};
A0 = gridifyState(A0); 

S0 = timetable(seconds(0), 0,0, 0,0, inf,inf, inf,inf); % default state
S0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl', ...
    'tgt_px_fx', 'tgt_py_fx', 'tgt_px_lo', 'tgt_py_lo', 'tgt_px_hi', 'tgt_py_hi'};
S0 = gridifyState(S0);

Aspace = repmat({A0}, 1, L^2); % all possible actions 
ind = 1;
for eyeX = 0:(L-1)
    for eyeY = 0:(L-1)
        A1 = A0; 
        A1.eye_px_filt_trl = eyeX; 
        A1.eye_py_filt_trl = eyeY; 
        Aspace{ind} = A1;
        ind = ind+1;
    end
end

Sspace = repmat({S0}, 1, L^6); % all possible states 
ind = 1;
%for eyeX = 0:(L-1)
%    for eyeY = 0:(L-1)
        for fixX = 0:(L-1)
            for fixY = 0:(L-1)
                for loX = 0:(L-1)
                    for loY = 0:(L-1)
                        for hiX = 0:(L-1)
                            for hiY = 0:(L-1)
                                S1 = S0; 
                                %S1.eye_px_filt_trl = eyeX;
                                %S1.eye_py_filt_trl = eyeY;
                                S1.tgt_px_fx = fixX; 
                                S1.tgt_py_fx = fixY; 
                                S1.tgt_px_lo = loX; 
                                S1.tgt_py_lo = loY;
                                S1.tgt_px_hi = hiX; 
                                S1.tgt_py_hi = hiY;
                                Sspace{ind} = S1;
                                ind = ind+1; 
                            end
                        end
                    end
                end
            end
        end
%    end
%end

% start with a random state. Pick this by taking a random action from the
% default state. 
S0 = updateGridState(S0, Aspace{randi(L^2)});

% init Q table 
Qvals = rand(length(Sspace), length(Aspace)); Qvals = num2cell(Qvals);
Qtable = [nan, Aspace; Sspace', Qvals];

%% iteration 
Snow = S0; 
Sall = repmat(S0, tf, 1); Aall = repmat(A0, tf, 1);
for t = 1:tf
    Qs = getQ(Snow, []); [q1,indMax] = max(Qs); At = Qtable{1,indMax+1};
    Snxt = updateGridState(Snow, At); rew = R(Snxt);
    Qs = getQ(Snxt, []); [maxQ,indMax] = max(Qs); 
    q2 = rew + gamma * maxQ;
    setQ(Snow,At, (1-alpha)*q1 + alpha*q2);
    Sall(t,:) = Snow; Aall(t,:) = At;
    Snow = Snxt;
end

Qfun = @(s,a) getQ(s,a);

%% helper
    function q = getQ(s,a)
        if isempty(a)
            row = findQr(s);
            q = [Qtable{row,2:end}]; 
        elseif isempty(s)
            col = findQc(a);
            q = [Qtable{2:end,col}];
        else
            row = findQr(s);
            col = findQc(a);
            q = Qtable{row,col};
        end
    end

    function setQ(s,a,q)
        row = findQr(s);
        col = findQc(a);
        Qtable{row,col} = q;
    end

    function r = findQr(s)
        r = 0;
        stateFound = false;
        while ~stateFound & (r < length(Sspace))
            r = r+1;
            stateFound = isEqualState(s, Sspace{r});
        end
        if ~stateFound 
            warning('State not found in table')
            r = 1; % S0
        end
        r = r+1;
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
            c = 2; % A0
        end
    end

end