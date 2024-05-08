function [Qfun, Qtable, Sall, Aall] = ReinforcementLearnGrid3(wT, gamma, alpha, tf, L, Sspace, Aspace, Qtable)
% solve all possibilities of Bellman equations given the reward
% function R(s) = wT*Phi(s)
% Return Q as a table with possible states as rows and actions as columns 
% and as a function of [state, action] 

%% initialize 
Qtable = single(zeros(size(Qtable))); % ?
Qtable2 = Qtable; 

alpha = single(alpha); gamma = single(gamma);

T = uint8(Qtable2); % state transitions 

%% get all state transitions 
disp('Determining state transitions...')
for r = 1:length(Sspace)
    Snow = Sspace{r}; 
    for c = 1:length(Aspace)
        At = Aspace{c};
        Snxt = updateGridState(Snow,At);
        r2 = state2ind(Snxt,L);
        T(r,c) = r2;
    end
end

%% get all state rewards 
disp('Determining state rewards...')
R = single(zeros(size(Sspace))); 
for r = 1:length(R)
    Snow = Sspace{r}; 
    unwrappedPhi = phiGrid(Snow); unwrappedPhi = unwrappedPhi{:,:}';
    R(r) = single(wT*unwrappedPhi);
end

%% reverse iteration 
disp('Reverse Iterating...')
for t = fliplr(1:tf)
    for r = 1:length(Sspace)
        for c = 1:length(Aspace)
            r2 = T(r,c);
            Qs = Qtable(r2,:); [maxQ,indMax] = max(Qs);
            rew = R(r2);
            Qtable2(r,c) = rew + gamma*maxQ; 
        end
    end
    Qtable = alpha*Qtable2 + (single(1)-alpha)*Qtable;
end
Qfun = @(s,a) getQ(s,a);

%% forward iteration 
disp('Forward Iterating...')
S0 = Sspace{randi(length(Sspace))}; 
Sall = repmat(S0, tf, 1); Aall = repmat(A0, tf, 1);
Snow = S0; r = state2ind(Snow,L);
for t = 1:tf
    Qs = Qtable(r,:); [q1,indMax] = max(Qs); At = Aspace{indMax};
    r2 = T(r,indMax); Snxt = Sspace{r2};
    Sall(t,:) = Snow; Aall(t,:) = At;
    Snow = Snxt; r = r2;
end

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