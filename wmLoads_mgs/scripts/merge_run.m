function ii_perrun=merge_run(ii_trial,run,taskinfo,same_field) 
% function for merge each run's different epoch data
% data in each epoch under the same filename are merged into 1 column in 
% table, with that specific filename as column name
% ii_perrun.fieldname could index the whole column

% Qingqing Yang, 6/16/2022

if nargin <4
% fields that are same in diff resp_epoch
same_field = {'X','Y','Pupil','XDAT','Velocity','i_saccall','f_saccall',...
    'i_sacc_allerr','f_sacc_allerr','n_sacc_epoch','calib_amt','calib_adj',...
    'calib_err','drift_amt','alltarg'};
% same_field = [];
end

ii_run = ii_trial{run};
fields = [];
for i = 1:size(ii_run,2)
    fields{i} = fieldnames(ii_run{i}); % 1 by resp_epochnum cell, each as fields by 1 cell
end

tf = 1==1;
for i = 1:size(ii_run,2)-1
    if ~strcmp(fields{i},fields{i+1})
        tf = 0==1;
    end
end

if tf == 0
    error("fieldnames in structs to combine are different, check input");
else
    ii_perrun = table;
    height = taskinfo.task.ntrials;
    ii_perrun.run = ones(height,run);
    ii_perrun.trial = [1:height]';
    ii_perrun.is_eyetracker = repmat(taskinfo.task.is_eyetracker,height,1);
    ii_perrun.is_TMS = repmat(taskinfo.task.is_TMS,height,1);
    ii_perrun.subject = repmat(taskinfo.task.subjectID,height,1);

    epoch_num = size(ii_run,2);
    for fn = 1:size(fields{1},1)
        field=fields{1}{fn};
        
        if strcmp(field,'params')
            ii_perrun.params = repmat(ii_run{1,1}.params,height,1);
            for ep = 2:epoch_num
                ii_perrun.params = [ii_perrun.params repmat(ii_run{1,ep}.params,height,1)];
            end
        elseif sum(strcmp(field,same_field))==0
            ii_perrun.(field) = [ii_run{1,1}.(field)];
            for ep = 2:epoch_num
                ii_perrun.(field) = [ii_perrun.(field) ii_run{1,ep}.(field)];
            end
        elseif sum(strcmp(field,same_field))>0 % only needa merge once
            ii_perrun.(field) = [ii_run{1,1}.(field)];  
        end
        
    end
    clear fn;
    
    
end

end