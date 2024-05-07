function sNorm = phiGrid(sGrid)
% normalize a grid state (or action) to the space [0, 1]

sNorm = timetable(sGrid.Time);

[~,posOpts] = gridifyState();
L = length(posOpts); % sGrid entries will be in range [0, L] 

for c = 1:width(sGrid)
    varname = sGrid.Properties.VariableNames{c};
    gridind = sGrid{:,c};
    normvals = double(gridind)/L;
    normvals = table(normvals); normvals.Properties.VariableNames = {varname};
    sNorm = [sNorm, normvals];
end

end