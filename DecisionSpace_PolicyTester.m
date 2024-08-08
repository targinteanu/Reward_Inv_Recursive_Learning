for qi = fliplr(1:length(Qtbls))

qi
Qtbl2 = Qtbls{qi};

%% pick best action at each state under this policy 
Abest = zeros(size(StateNames)); 
for ind3 = 1:length(Abest)
    [~,Abest(ind3)] = max(Qtbl2(ind3,:));
end

%% test accuracy of this policy against training data 
disp('Computing training accuracy...')
train_result = false(height(Srec),1);
for trl = 1:width(Srec)
    Strl = Srec(:,trl); Atrl = Arec(:,trl);
    StInd = find(Strl); 
    Apred = Abest(StInd);
    train_result(trl) = Apred == Atrl;
end
training_accuracy = mean(train_result)

%% test accuracy of this policy against other data
%{
disp('Computing testing accuracy...')
indTest = randi(length(data_all_trials)); 
dtasTest = data_all_trials{ind1};
test_result = logical([]);
for ind2 = 1:(length(dtasTest))
    dta = dtasTest(ind2);
    disp(['Testing ',dta.name,' ...']);
    for trl = 1:(length(dta.task_cond))
        [Strl, Atrl] = getStateSpace(dta,trl,false);
        Strl = gridifyState(Strl); Atrl = gridifyState(Atrl);
        if ~isempty(Strl)
            for trl = 1:height(Strl)
                Strl = Strl(trl,:); Atrl = Atrl(trl,:);
                StInd = state2ind(Strl,L);
                if ~isEqualState(Strl, Qstate{StInd})
                    error('State does not match index of table.')
                end
                Apred = Qaction{Abest(StInd)};
                test_result = [test_result; isEqualState(Apred, Atrl)];
            end
        end
    end
end
testing_accuracy = mean(test_result)
clear dta 
%}
end