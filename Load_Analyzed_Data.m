function data_all_trials = Load_Analyzed_Data(basedirs)
% load data to array 

f_name_end = 'ANALYZED.mat';
f_name_str = 'reward_prediction';

if nargin == 0
    basedirs = [];
end

if isempty(basedirs)
    basedirs = {'Y:\Ephys\data_132F\2024-03\', 'Y:\Ephys\data_132F\2024-04\'}; 
end

data_all_trials = {};

for bd = basedirs
    data_this_dir = {};
    basedir = bd{:};
    sessions = dir(basedir);
    sessions = sessions(3:end); % remove '.', '..' 
    for session = sessions'
        if ~isempty(session) & session.isdir
            data_this_session = [];
            recs = dir([session.folder,filesep,session.name]);
            recs = recs(3:end); % remove '.', '..'
            for rec = recs'
                if ~isempty(rec) & rec.isdir
                    filepath = [rec.folder,filesep,rec.name,filesep,'analyzed_data',filesep,'behavior_data',filesep,'eye'];
                    if isfolder(filepath)
                        files = dir(filepath);
                        for f = files'
                            if ~isempty(f) & ~f.isdir
                                if strcmpi(f_name_end, f.name((end-length(f_name_end)+1):end) ) & ...
                                   strcmpi(f_name_str, f.name(1:length(f_name_str)) )
                                    load([f.folder,filesep,f.name], 'trials_data');
                                    tr_d = struct('task_cond', trials_data.task_cond, ...
                                                  'tgt_cond', trials_data.tgt_cond, ...
                                                  'rew_cond', trials_data.rew_cond, ...
                                                  'jump_cond', trials_data.jump_cond, ...
                                                  'choice', trials_data.choice, ...
                                                  'cue_x', trials_data.cue_x, ...
                                                  'cue_y', trials_data.cue_y, ...
                                                  'tgt_px', trials_data.tgt_px, ...
                                                  'tgt_py', trials_data.tgt_py, ...
                                                  'eye_px_filt', trials_data.eye_px_filt, ...
                                                  'eye_py_filt', trials_data.eye_py_filt, ...
                                                  'eye_vx_filt', trials_data.eye_vx_filt, ...
                                                  'eye_vy_filt', trials_data.eye_vy_filt ...
                                                  );
                                    data_this_session = [data_this_session, tr_d];
                                    clear trials_data sac_data meta_data tr_d
                                end
                            end
                        end
                    end
                end
            end
            data_this_dir = [data_this_dir, data_this_session];
            clear data_this_session
        end
    end
    data_all_trials = [data_all_trials, data_this_dir];
    clear data_this_dir
end

end