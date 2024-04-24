gamma = .5; % discount factor 

% loop through all 
ind1 = randi(length(data_all_trials)); 
%for ind1 = 1:(length(data_all_trials))
    dtas = data_all_trials{ind1};
    ind2 = randi(length(dtas));
%    for ind2 = 1:(length(dtas))
        dta = dtas(ind2);
        dist_from_prim = zeros(length(dta.task_cond),1); 
        Srec = []; Arec = [];
        for trl = 1:(length(dta.task_cond))
            [Strl, Atrl] = getStateSpace(dta,trl,false);
            Srec = [Srec; Strl]; 
            Arec = [Arec; Atrl];

            Phi = Strl{1:height(Strl), 1:width(Strl)};
            Gamma = gamma.^(0:(width(Strl)-1));
            muEi = Phi*Gamma';
        end

        Srec.Time = Srec.Time - Srec.Time(1);
        Arec.Time = Arec.Time - Arec.Time(1);
%    end
%end