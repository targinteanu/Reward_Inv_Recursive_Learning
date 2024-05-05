function [A, Q] = pickBestAct(S, Qtbl)

stateSpace = Qtbl(:,1); actSpace = Qtbl(1,:); 
stateSpace = stateSpace(2:end); actSpace = actSpace(2:end); 
Q = cell2mat(Qtbl(2:end,2:end));

r = 0;
stateFound = false; 
while ~stateFound & (r < length(stateSpace))
    r = r+1;
    stateFound = isEqualState(S, stateSpace{r});
end

if ~stateFound 
    warning('State not found in table')
    A = actSpace{1}; % move by [0, 0] 
    Q = [];
else
    Q = Q(r,:);
    [Q,c] = max(Q); 
    A = actSpace{c};
end

end