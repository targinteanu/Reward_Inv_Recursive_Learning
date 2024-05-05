function [gridState, posOpts] = gridifyState(contState)
% take a continuous state and turn it into a discrete grid 
% works for actions as well as states 

if nargin == 0
    contState = [];
end

if ~isempty(contState)
    gridState = timetable(contState.Time);
else
    gridState = [];
end

% each position coordinate has 39 possible states: 
% off-grid, [-9:.5:9], off-grid
% mapped to integer indexes 

minpos = -9; maxpos = 9; steppos = 3; 
posOpts = minpos:steppos:maxpos; 
L = length(posOpts);

for c = 1:width(contState)
    varname = contState.Properties.VariableNames{c};
    gridvals = uint8(zeros(height(contState),1));
    for r = 1:height(contState)
        contval = contState{r,c};
        if contval > maxpos + steppos/2
            % off-grid + 
            gridvals(r) = L+1;
        elseif contval < minpos - steppos/2
            % off-grid -
            gridvals(r) = 0;
        else
            % on grid 
            [~,gridvals(r)] = min(abs(contval - posOpts));
        end
    end
    gridvals = table(gridvals); gridvals.Properties.VariableNames = {varname};
    gridState = [gridState, gridvals];
end

end