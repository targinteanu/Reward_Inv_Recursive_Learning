%% init rows/cols of Q table 

[~,posOpts] = gridifyState();
posOpts = [-inf,posOpts,inf];
L = length(posOpts); % position takes L possible values from 0 to L-1 (base L digit) 

A0 = timetable(seconds(0), 0, 0); % null action 
A0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl'};
A0 = gridifyState(A0); 

S0 = timetable(seconds(0), 0,0, 0,0, inf,inf, inf,inf); % default state
S0.Properties.VariableNames = {'eye_px_filt_trl', 'eye_py_filt_trl', ...
    'tgt_px_fx', 'tgt_py_fx', 'tgt_px_lo', 'tgt_py_lo', 'tgt_px_hi', 'tgt_py_hi'};
S0 = gridifyState(S0);

disp('Getting action space...')

Qaction = repmat({A0}, 1, L^2); % all possible actions 
for eyeX = 0:(L-1)
    for eyeY = 0:(L-1)
        A1 = A0; 
        A1.eye_px_filt_trl = eyeX; 
        A1.eye_py_filt_trl = eyeY; 
        ind = act2ind(A1,L);
        if ind > length(Qaction)
            error('Action exceeds table size.')
        end
        Qaction{ind} = A1;
    end
end

disp('Getting state space...')

Qstate = repmat({S0}, 1, L^6); % all possible states 
%for eyeX = 0:(L-1)
%    for eyeY = 0:(L-1)
        for fixX = 0:(L-1)
            for fixY = 0:(L-1)
                for loX = 0:(L-1)
                    for loY = 0:(L-1)
                        for hiX = 0:(L-1)
                            for hiY = 0:(L-1)
                                S1 = S0; 
                                %S1.eye_px_filt_trl = eyeX;
                                %S1.eye_py_filt_trl = eyeY;
                                S1.tgt_px_fx = fixX; 
                                S1.tgt_py_fx = fixY; 
                                S1.tgt_px_lo = loX; 
                                S1.tgt_py_lo = loY;
                                S1.tgt_px_hi = hiX; 
                                S1.tgt_py_hi = hiY;
                                ind = state2ind(S1,L); 
                                if ind > length(Qstate)
                                    error('State exceeds table size.');
                                end
                                Qstate{ind} = S1;
                            end
                        end
                    end
                end
            end
        end
%    end
%end

clear eyeX eyeY fixX fixY loX loY hiX hiY A1 S1

disp('State-Action space initialized.')

%% get all state transitions 
disp('Determining state transitions...')
LS = length(Qstate); LA = length(Qaction);
[C,R] = meshgrid(1:LA, 1:LS);
iter2c = C(:); iter2r = R(:); % flatten 
LRC = length(iter2c);
T = uint8(zeros(LS, LA)); % state transitions 
T = T(:); % flatten 
parfor iter = 1:LRC
    r = iter2r(iter); % unflatten
    c = iter2c(iter); % unflatten
    Snow = Qstate{r}; 
    At = Qaction{c};
    Snxt = updateGridState(Snow,At);
    r2 = state2ind(Snxt,L);
    T(iter) = r2;
end
T = reshape(T,[LS,LA]);
clear r c r2 Snow Snxt At LS LA LSA C R iter2r iter2c

%% save 
save('StateQT')