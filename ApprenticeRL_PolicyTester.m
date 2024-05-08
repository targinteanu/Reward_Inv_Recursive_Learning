for qi = fliplr(1:size(Qtbls,3))

qi
Qtbl2 = Qtbls(:,:,qi);

%% pick best action at each state under this policy 
Abest = zeros(size(Qstate)); 
for ind3 = 1:length(Abest)
    [~,Abest(ind3)] = max(Qtbl2(ind3,:));
end

%% test accuracy of this policy against training data 
disp('Computing training accuracy...')
train_result = false(height(Srec),1);
for t = 1:height(Srec)
    St = Srec(t,:); At = Arec(t,:);
    StInd = state2ind(St,L); 
    if ~isEqualState(St, Qstate{StInd})
        error('State does not match index of table.')
    end
    Apred = Qaction{Abest(StInd)};
    train_result(t) = isEqualState(Apred, At);
end
training_accuracy = mean(train_result)

%% test accuracy of this policy against other data
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
            for t = 1:height(Strl)
                St = Strl(t,:); At = Atrl(t,:);
                StInd = state2ind(St,L);
                if ~isEqualState(St, Qstate{StInd})
                    error('State does not match index of table.')
                end
                Apred = Qaction{Abest(StInd)};
                test_result = [test_result; isEqualState(Apred, At)];
            end
        end
    end
end
testing_accuracy = mean(test_result)
clear dta 

end