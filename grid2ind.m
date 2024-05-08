function ind = grid2ind(S, n, removeEye)
% produce an index (>= 1) that is unique to the grid-world coords
% S is a single state, not a list
% n is the base (number of possible values of each variable)
S = S{1,:};
S = num2str(S(:))';
if removeEye
    S = S(3:end); % remove eye position
end
ind = base2dec(S,n);
ind = ind+1;
end