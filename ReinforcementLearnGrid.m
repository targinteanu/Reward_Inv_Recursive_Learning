function [Qfun, Qtable, Sall, Aall] = ReinforcementLearnGrid(wT, gamma, alpha, tf, L, Sspace, Aspace, Qtable)
% use Q learning to determine the optimal value function Q given the reward
% function R(s) = wT*Phi(s)
% Return Q as a table with possible states as rows and actions as columns 
% and as a function of [state, action] 

%% initialize s and Q  
A0 = timetable(seconds(0), 0, 0); % null action 
A0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl'};
A0 = gridifyState(A0); 

S0 = timetable(seconds(0), 0,0, 0,0, inf,inf, inf,inf); % default state
S0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl', ...
    'tgt_px_fx', 'tgt_py_fx', 'tgt_px_lo', 'tgt_py_lo', 'tgt_px_hi', 'tgt_py_hi'};
S0 = gridifyState(S0);

% start with a random state. Pick this by taking a random action from the
% default state. 
S0 = updateGridState(S0, Aspace{randi(length(Aspace))});

%% iteration 
Snow = S0; 
Sall = repmat(S0, tf, 1); Aall = repmat(A0, tf, 1);
for t = 1:tf
    Qs = getQ(Snow, []); [q1,indMax] = max(Qs); At = Aspace{indMax};
    Snxt = updateGridState(Snow, At); 
    unwrappedPhi = phiGrid(Snxt); unwrappedPhi = unwrappedPhi{:,:}'; rew = wT*unwrappedPhi;
    Qs = getQ(Snxt, []); [maxQ,indMax] = max(Qs); 
    q2 = rew + gamma * maxQ;
    setQ(Snow,At, (1-alpha)*q1 + alpha*q2);
    Sall(t,:) = Snow; Aall(t,:) = At;
    Snow = Snxt;
end

Qfun = @(s,a) getQ(s,a);

%% helper

    function [q, row, col] = getQ(s,a)
        row = []; col = [];
        if isempty(a)
            row = findQr(s);
            q = [Qtable(row,:)]; 
        elseif isempty(s)
            col = findQc(a);
            q = [Qtable(:,col)];
        else
            row = findQr(s);
            col = findQc(a);
            q = Qtable(row,col);
        end
    end

    function setQ(s,a,q)
        row = findQr(s);
        col = findQc(a);
        Qtable(row,col) = single(q);
    end

    function r = findQr(s)
        r = state2ind(s,L); 
        if ~isEqualState(s, Sspace{r})
            warning('State not at expected place in table.')
            r = 0;
            stateFound = false;
            while ~stateFound & (r < length(Sspace))
                r = r+1;
                stateFound = isEqualState(s, Sspace{r});
            end
            if ~stateFound
                warning('State not found in table.')
                r = 1; % S0
            end
        end
    end

    function c = findQc(a)
        c = act2ind(a,L);
        if ~isEqualState(a, Aspace{c})
            warning('Action not at expected place in table.')
            found = false;
            c = 1;
            while ~found & (c <= length(Aspace))
                if isEqualState(a, Aspace{c})
                    % found action (col)
                    found = true;
                else
                    c = c+1;
                end
            end
            if ~found
                % found row but not column
                warning('unexpected action')
                c = 1; % A0
            end
        end
    end

end