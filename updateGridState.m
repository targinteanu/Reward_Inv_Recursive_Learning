function s2 = updateGridState(s1, a1)
% updateState, but both input state and action and output state are in grid world.  
% there may be a better way to do this

s2 = gridifyState(updateState(ungridState(s1), ungridState(a1)));

end