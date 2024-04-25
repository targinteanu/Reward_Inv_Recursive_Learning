function result = isEqualState(s1, s2)

result = isequal(s1{:,:}, s2{:,:}) & ...
         ... isequal(s1.Time, s2.Time) & ...
         isequal(s1.Properties.VariableNames, s2.Properties.VariableNames);

end