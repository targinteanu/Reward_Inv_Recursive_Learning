%% test random 
ind1 = randi(length(data_all_trials));
dtas = data_all_trials{ind1};
ind2 = randi(length(dtas)); 
dta = dtas(ind2);
trl = randi(length(dta.task_cond));
%trl = trl+1;

[S,A] = getStateSpace(dta, trl, true);
dta.reward_area(trl)

%% gather all possibilities 
reward_areas = [];
Sall = []; Aall = [];
for ind1 = 1:(length(data_all_trials))
    dtas = data_all_trials{ind1};
    Ssess = []; Asess = [];
    for ind2 = 1:(length(dtas))
        dta = dtas(ind2);
        Srec = []; Arec = [];
        for trl = 1:(length(dta.task_cond))
            [Strl, Atrl] = getStateSpace(dta,trl,false);
            Srec = [Srec; Strl]; 
            Arec = [Arec; Atrl];
            reward_areas = [reward_areas; dta.reward_area(trl)];
        end
        Ssess = [Ssess; Srec];
        Asess = [Asess; Arec];
    end
    Sall = [Sall; Ssess]; 
    Aall = [Aall; Asess];
end

%% summarize results 
Senv = Sall(:,3:end); % environment only 
figure('Units','normalized', 'Position',[.1,.1,.5,.5]);
disp('reward_area [mean, sd]: ');
[mean(reward_areas), std(reward_areas)]
histogram(reward_areas); hold on; 
for col = 1:width(Senv)
    colname = Senv.Properties.VariableNames{col};
    coldata = Senv{:,col};
    coldata = coldata(~isinf(coldata));
    histogram(coldata);
    disp([colname,' [mean, sd]: ']);
    [mean(coldata), std(coldata)]
end
legend(['reward_area', Senv.Properties.VariableNames]);