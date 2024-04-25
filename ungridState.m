function contState = ungridState(gridState)
% take a discrete grid state and estimate the continuous state space 
% works for both states and actions 

contState = timetable(gridState.Time);

[~,posOpts] = gridifyState([]); 
posOpts = [-inf, posOpts, inf]; % off-grid maps to inf

for c = 1:width(gridState)
    varname = gridState.Properties.VariableNames{c};
    gridind = gridState{:,c}+1; % ind 1 -> -inf
    contvals = posOpts(gridind)';
    contvals = table(contvals); contvals.Properties.VariableNames = {varname};
    contState = [contState, contvals];
end

end