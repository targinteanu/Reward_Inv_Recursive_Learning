function sGrid = phiGridInv(sNorm)
% recover a state (or action) from normalization by phiGrid 

sGrid = timetable(sNorm.Time); 

[~,posOpts] = gridifyState();
L = length(posOpts); % sGrid entries will be in range [0, L] 

for c = 1:width(sNorm)
    varname = sNorm.Properties.VariableNames{c};
    normvals = sNorm{:,c};
    gridind = uint8(normvals*L);
    gridind = table(gridind); gridind.Properties.VariableNames = {varname};
    sGrid = [sGrid, gridind];
end

end