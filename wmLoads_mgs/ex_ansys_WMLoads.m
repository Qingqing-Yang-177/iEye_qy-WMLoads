% ex_anasys_WMLoads.m

% example analysis script for example 1/2/4 target eyetracking dataset
% (eyedata*.edf) that go through preprocessing, saccade extraction,
% MGS scoring, QC plotting, and some basic analyses/plots
% =======================
% Revised from Tommy Sprague, 6/12/2018
% which could be found at https://github.com/tommysprague/iEye_ts

% Adapted the package for multiple target items eye-tracking data processes
% The analysis is applicable to 1 target, 2 target, to 4 target task.
% For multiple target task settings, the target message should be sent as
% TarX1, TarY1,TarX2, TarY2, TarX3, TarY3,TarX4 and TarY4.

% The order of reporting whole 4 targets are freely chosen by participants,

% The code is applicable when feedbacks for all targets are shown either
% simultaneously after all responses are finished in each trial, or
% separately after each response

% For Separate Feedback, select ii_preproc_WMLoads_SepFb during the
% pre-processing. The calibration loops over all the feedbacks.

% ========================================
% There are sample data from three control experiment:

% task = 'targ1_control' is the control exp where all trials have only 1
% target at 45 degree, 1st quard.During resp_epoch, target would not be
% offset

% task = 'targ2_control' is the control exp where 1/2 trials have 1 target,
% which is at 1st or 2nd quard, and 1/2 trails have 2 target at 1st & 2nd
% quards

% task = 'targ4_control' is the control exp where 1/3 trials have 1 targt
% which is at 1st quard, 1/3 trials have 2 targets in 1st & 2nd quards, and
% 1/3 trials have 4 targets one in each quard.

% Before running, to do list:
% change task type, set epoch and params for specific task
% change path in ex_ansys_WMLoads.m
% change behavior filename
% change edf2asc path in init.m
% change edf2asc path in import_edf.m
% change exp pc settings
% change ii_loadparams.m params

% Qingqing Yang, 6/23/2022

%%
addpath(genpath('D:\A_Clay_Lab\iEye_qy-WMLoads'));

file_dir ='D:\A_Clay_Lab\iEye_qy-WMLoads\wmLoads_mgs';

ii_params = ii_loadparams; % load default set of analysis parameters, only change what we have to

% % Define the task
% task = 'targ1_control';
% task = 'targ2_control';
task = 'targ4_control'; % SeperateFb, back to fixation after each resp
% task = 'wmloads_SimuFb'; % For formal test edf

ii_params.task = task;
if strcmp(task,'targ1_control')
    ifg_fn = 'D:\A_Clay_Lab\iEye_qy-WMLoads\wmLoads_mgs\p_500hz.ifg';
    ii_params.valid_epochs =[1 2 3 4 5];
    ii_params.trial_end_value = 5;   % XDAT value for trial end
    ii_params.drift_epoch = [1 2]; % XDAT values for drift correction
    ii_params.calibrate_epoch = 4;
    ii_params.response_epoch = 3;
    ii_params.plot_epoch = [2 3 4];
    edf_prefix ='eyedata_OTOYFix';
    ii_params.calibrate_target = {'TarX','TarY'};
    
elseif strcmp(task,'targ2_control')
    ifg_fn = 'D:\A_Clay_Lab\iEye_qy-WMLoads\wmLoads_mgs\p_500hz_TarX12.ifg';
    ii_params.valid_epochs =[1 2 3 4 5 6 7 8];
    ii_params.trial_end_value = [5 8];   % XDAT value for trial end
    ii_params.drift_epoch = [1 2]; % XDAT values for drift correction
    ii_params.calibrate_epoch = [4 7];
    ii_params.response_epoch = [3 6];
    ii_params.plot_epoch = [2 3 4 5 6 7 8];
    edf_prefix ='eyedata_OTTTDLYFix';
    ii_params.calibrate_target = {'TarX1','TarY1','TarX2','TarY2'};
    
elseif strcmp(task,'targ4_control')
    ifg_fn = 'D:\A_Clay_Lab\iEye_qy-WMLoads\wmLoads_mgs\p_500hz_TarX1234.ifg';
    ii_params.valid_epochs =[1 2 3 4 5 6 7 8 9 10 11];
    ii_params.trial_end_value = 11;   % XDAT value for trial end
    ii_params.drift_epoch = [1 2];
    ii_params.calibrate_epoch = [3 5 7 9];   % XDAT value for when we calibrate (feedback stim)
    ii_params.response_epoch = [3 5 7 9];
    ii_params.plot_epoch = [2 3 4 5 6 7 8 9 10];  % what epochs do we plot for preprocessing?
    edf_prefix ='eyedata_CTTF';
    ii_params.calibrate_target = {'TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'};
 
else
    error('Unknown task type, please define specific epochs & params');
end

ii_params.calibrate_select_mode = 'last'; % how do we select fixation with which to calibrate?
ii_params.calibrate_mode = 'scale'; % scale: trial-by-trial, rescale each trial; 'run' - run-wise polynomial fit
ii_params.blink_window = [200 200]; % how long before/after blink (ms) to drop?
ii_params.calibrate_limits = [2.5]; % when amount of adj exceeds this, don't actually calibrate (trial-wise); ignore trial for polynomial fitting (run)
ii_params.calibrate_window =200; % 100 for 1000Hz, 200 for 500 Hz
ii_params.resolution = [1920 1080]; % Curtis lab, behavior room
ii_params.ppd = 32.2029; % % Curtis lab, behavior room
ii_params.aperture_size = 12; % 12 dva as the aperture_size

% all files we want to load are eyedata_*.mat - so let's list them (we
% know they're in the same directory as this script)
tmp = file_dir; tmp2 = strfind(tmp,filesep);
root = tmp;
edf_files = dir(fullfile(root,task,sprintf('%s*.edf',edf_prefix)));

% Load behavior data, needs change if more than 1 participant
taskinfo_prefix = 'MGS_WMLoads_Simul_Control_CTTF'; % targ4_control
task_files = dir(fullfile(root,task,sprintf('%s*.mat',taskinfo_prefix)));
taskinfo = cell(1,length(task_files));
for i = 1:length(task_files)
    taskinfo{i}=load(task_files(i).name);
end


% create empty cell array of all our trial data for combination later on
ii_trial = cell(1,length(edf_files));
ii_alltrial = cell(1,length(edf_files));

%% Delete extraneous .asc files in this directory
file_list = ls(); file_list = string(file_list);
files_to_delete = contains(file_list,'.asc');
if ~isempty(files_to_delete)
    files_to_delete = file_list(files_to_delete,:);
    for ii = 1:length(files_to_delete)
        f = files_to_delete{ii};
        delete(f);
    end
end

%% Preprocess and score the error
for ff = 1:length(edf_files)
    
    % what is the output filename?
    preproc_fn = sprintf('%s/%s_preproc.mat',edf_files(ff).folder,edf_files(ff).name(1:(end-4)));
    
    % run preprocessing!
    if strcmp(task,'targ2_control')|| strcmp(task,'targ4_control')
        [ii_data, ii_cfg, ii_sacc] = ii_preproc_WMLoads_SepFb(fullfile(edf_files(ff).folder,edf_files(ff).name),ifg_fn,preproc_fn,ii_params);
    else
        [ii_data, ii_cfg, ii_sacc] = ii_preproc_WMLoads_SimulFb(fullfile(edf_files(ff).folder,edf_files(ff).name),ifg_fn,preproc_fn,ii_params);
    end
    
    if ff == 1
        % plot some features of the data
        % (check out the docs for each of these; lots of options...)
        
        if contains(ii_cfg.vis,'TarX4')
            plot_coords = {'X','Y','TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'};
        elseif contains(ii_cfg.vis,'TarX2')
            plot_coords = {'X','Y','TarX1','TarY1','TarX2','TarY2'};
        else
            plot_coords = {'X','Y','TarX','TarY'};
        end
        
        ii_plottimeseries_WMLoads(ii_data,ii_cfg,plot_coords); % pltos the full timeseries
        
        ii_plotalltrials_WMLoads(ii_data,ii_cfg,plot_coords);
        
        % fixations from beginning to the end sample of each trial
        ii_plotalltrials2d_WMLoads(ii_data,ii_cfg,{'X','Y'}); % plots all trials, in 2d, overlaid on one another w/ fixations
        
    end
    
    % score trials
    % default parameters should work fine - but see docs for other
    % arguments you can/should give when possible
    pp = 1; % the id of the target
    for tt= 1:length(ii_params.response_epoch)
        resp_epoch = ii_params.response_epoch(tt);
        
        if contains(ii_cfg.vis,'TarX4')
            targ_coords = {'TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'};
        elseif contains(ii_cfg.vis,'TarX2')
            targ_coords = {'TarX1','TarY1','TarX2','TarY2'};
        else
            targ_coords = {'TarX','TarY'};
        end
        [ii_trial{ff}{pp},ii_cfg] = ii_scoreMGS_WMLoads(ii_data,ii_cfg,...
            ii_sacc,targ_coords,resp_epoch,ii_params.drift_epoch); % {1 by run}{1 by epoch}
        
        % save each trial's exp setup
        ii_trial{ff}{pp}.is_TMS = ones(size(ii_trial{ff}{pp}.X,1),1)==taskinfo{ff}.task.is_TMS;
        ii_trial{ff}{pp}.is_eyetracker =  ones(size(ii_trial{ff}{pp}.X,1),1)==taskinfo{ff}.task.is_eyetracker;
        
        % plot fixations in each resp epoch separately
        % [ii_data2,ii_cfg2] = ii_definetrial_respepoch_WMLoads(ii_data,ii_cfg,...
        %   ii_params.trial_start_chan, ii_params.trial_start_value,...
        %   ii_params.trial_end_chan,   ii_params.trial_end_value,...
        %   ii_params.drift_epoch, resp_epoch);
        %
        %   ii_plotalltrials2d_WMLoads(ii_data2,ii_cfg2,{'X','Y'});
        
        pp=pp+1;
    end
    clear tt pp;
    
end

%% combine the data based on each run (concat the resp_epoches)
for ff=1:length(edf_files)
    ii_alltrial{ff} = ii_combineresp(ii_trial{ff});
    ii_alltrial{ff}.params.resp_epoch = ii_params.response_epoch;
end

%% Combine the data in each run (concat the runs)
clear ii;
ii_sess = ii_combineruns_WMLoads(ii_alltrial);

% save setsize
ii_sess.setsize = [];
setsize = [];
MaxLoad = taskinfo{ff}.task.maxLoad; % 4 for my task
TrialNum = taskinfo{ff}.task.trialNum; %36 for my task
for r = 1:max(unique(ii_sess.r_num))
    for ii = 1:size(ii_sess.alltarg(ii_sess.r_num==r,:),1)/MaxLoad ;
        add = (r-1)*TrialNum*MaxLoad; 
        setsize(ii,1) = sum(ii_sess.alltarg(add+ii,:)~=0)/2; % since alltarg is for x and y
    end
    ii_sess.setsize=[ii_sess.setsize; repmat(setsize,[MaxLoad,1])];
    clear setsize
end

% save ii_sess as sumarized data
datafile = sprintf('%s/%s_%.fruns_ii_sess.mat',task_files(1).folder,taskinfo{1, 1}.task.subjectID,length(taskinfo));
save(datafile,'ii_sess');
    
%% look at the processed data, make sure it looks ok

% based on what criteria shoudl we exclude trials?
which_excl = [13 20 23]; % Added 23 when i_sacc and f_sacc is for diff target, Qingqing Yang, 6/9/2022

% TRIAL EXCLUSION CODES:
% first-digit:
% - 1 - trial-level exclusion (resp_epoch not exist in that trial [10];
%       bad drift correction [11], calibration [12], or delay-
%       fixation break [13]
% - 2 - primary saccade exclusion (no primary sacc [20]; too small/short [21], large error [22], inconsistent i&f saccade [23])

% first, plot an exclusion report over all runs
% - this will open a few figures. first one will be a 'dot plot', which
%   shows all exclusion criteria exceeded for a given trial by plotting a dot
%   above that criterion. if any dot for a trial is red, we're going to throw that
%   trial away. if black, we'll ignore that criterion for now.
% - second, plots a summary graph, showing % of trials excluded due to each
%   criterion, and how many overall will be excluded. If plotting
%   run-concatenated data, also show how this varies across runs

% To do:
% exclusion_plotQC: whether to assign 23 in all_excl based on setting

fh_excl = ii_plotQC_exclusions_WMLoads(ii_sess,ii_cfg,which_excl);

% second, plot every trial, scored - this will overlay the primary saccade
% (purple), final saccade (if one; green) on the raw channel traces for
% each trial. the trial # indicator will indicate whether we shoudl
% consider the trial carefully: if italics, at least one exclusion
% criterion reached. If red, we exclude trial based on a criterion. After
% the trial number, there are several letters indicating which, if any,
% exclusion criteria were passed, and an asterisk says whether that one was
% used to exclude.
% - d: drift correction
% - c: calibration
% - fb: fixation break
% - ei: iniital saccade error
% - xi: no initial saccade detected within response window
% - bi: bad initial saccade (duration/amplitude outside range)
% - ixf: initial & final saccade enpoints are closest to diff target
% Also plotted is the response window (vertical gray lines) and the target
% position (horizontal dashed lines)

fh_trial = ii_plotQC_alltrials_WMLoads(ii_sess,ii_cfg,which_excl);

% often, it's useful to save these out as pngs:
for ff = 1:length(fh_excl)
    fn2s = sprintf('%s/%s_QC_excl_%02.f.png',edf_files(1).folder,edf_prefix,ff);
    saveas(fh_excl(ff),fn2s);
end

for ff = 1:length(fh_trial)
    fn2s = sprintf('%s/%s_QC_trial_%02.f.png',edf_files(1).folder,edf_prefix,ff);
    saveas(fh_trial(ff),fn2s);
end

%% now some very basic analysis recipes
%  (RT histogram, all traces, and distribution of endpoints, etc)
%  below can be considered examples for how to plot a few typical
%  single-subj analyses of eye position data. Notice that the bulk of the
%  work has been done already (mostly in ii_score, but also in preproc,
%  etc). This makes our job very easy - just sorting, subplotting, etc. We
%  don't have different conditions in this experiment, so I sorted one of
%  the analyses by 'run' - this could easily be a condition variable that's
%  embedded within ii_sess.trialinfo instead. All the same logic applies.
%  Feel free to borrow from these recipies when conducting your own
%  analyses.


% decide which trials we want to include for analyses:
% - drop any trials that fail the which_excl test
% - make sure we're not using trials with faster than 100 ms RT

% check trial exclusions
% 1 in use_trial means no which_excl, use this trial
use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, which_excl), ii_sess.excl_trial, 'UniformOutput',false));
trial_exist = zeros(length(ii_sess.excl_trial),1)==0;
trial_exist=~cellfun(@any,cellfun(@(a) ismember(10,a),ii_sess.excl_trial,'UniformOutput',false));
use_trial = trial_exist&use_trial;
keep_rate = mean(use_trial)/mean(trial_exist);
num=length(ii_sess.excl_trial);

% check RT
use_trial(ii_sess.i_sacc_rt < 0.1) = 0;

% see how many trials we're keeping

fprintf('Analyzing %0.03f%% of trials\n%i out of %i\n',keep_rate*100, mean(use_trial)*num, mean(trial_exist)*num);

%% FIRST: plot RT histogram
use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, which_excl), ii_sess.excl_trial, 'UniformOutput',false));
trial_exist = zeros(length(ii_sess.excl_trial),1)==0;
trial_exist=~cellfun(@any,cellfun(@(a) ismember(10,a),ii_sess.excl_trial,'UniformOutput',false));
use_trial = trial_exist&use_trial;
use_trial(ii_sess.i_sacc_rt < 0.1) = 0;

f_ts = figure;
histogram(ii_sess.i_sacc_rt(use_trial==1),10);
xlabel('Response time (s)');
xlim([0 1]); % 800 ms is longest allowed RT
saveas(f_ts,sprintf('%s/%s_RT_his.png',edf_files(1).folder,taskinfo{1}.task.subjectID));
clear f_ts;

%% SECOND: plot all traces for primary, final saccade for each run
f_ts =figure;
ru = unique(ii_sess.r_num); % get run #'s

use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, which_excl), ii_sess.excl_trial, 'UniformOutput',false));
trial_exist = zeros(length(ii_sess.excl_trial),1)==0;
trial_exist=~cellfun(@any,cellfun(@(a) ismember(10,a),ii_sess.excl_trial,'UniformOutput',false));
use_trial = trial_exist&use_trial;
use_trial(ii_sess.i_sacc_rt < 0.1) = 0;

for rr = 1:length(ru)
    subplot(1,length(ru),rr); hold on;
    
    % only grab trials we'll use
    if strcmp(task,'wmloads_SimuFb')% only frist resp_epoch (SimuFb)
        thisidx = find(ii_sess.r_num==ru(rr) & ii_sess.resp_num==1 & use_trial==1);
    else
        thisidx = find(ii_sess.r_num==ru(rr) &  use_trial==1);
    end
    
    mycolors = lines(length(thisidx)); % color for each trial
    
    for tt = 1:length(thisidx)
        
        % plot trace(s) & endpoint(s)
        plot(ii_sess.i_sacc_trace{thisidx(tt)}(:,1),ii_sess.i_sacc_trace{thisidx(tt)}(:,2),'-','LineWidth',1.5,'Color',mycolors(tt,:));
        plot(ii_sess.i_sacc_raw(thisidx(tt),1),ii_sess.i_sacc_raw(thisidx(tt),2),'ko','MarkerFaceColor',mycolors(tt,:),'MarkerSize',5);
        if ~isempty(ii_sess.f_sacc_trace{thisidx(tt)})
            plot(ii_sess.f_sacc_trace{thisidx(tt)}(:,1),ii_sess.f_sacc_trace{thisidx(tt)}(:,2),'-','LineWidth',1.5,'Color',mycolors(tt,:));
            plot(ii_sess.f_sacc_raw(thisidx(tt),1),ii_sess.f_sacc_raw(thisidx(tt),2),'kd','MarkerFaceColor',mycolors(tt,:),'MarkerSize',5);
        end
        
        % plot target
        plot(ii_sess.targ(thisidx(tt),1),ii_sess.targ(thisidx(tt),2),'x','LineWidth',2,'MarkerSize',8,'Color',mycolors(tt,:));
    end
    plot(0,0,'ks','MarkerFaceColor','k','MarkerSize',5); % fixation point
    xlim([-15 15]); ylim([-15 15]);
    axis square off;
    
    title(sprintf('Run %i',ru(rr)));
end
saveas(f_ts,sprintf('%s/%s_i&f_traces.png',edf_files(1).folder,taskinfo{1}.task.subjectID));
clear f_ts;
% btw, a *very* fast way to do a quick approximation of the above, at least
% for one set of traces, is:
% cellfun(@(c) plot(c(:,1),c(:,2), ii_sess.i_sacc_trace{use_trial==1});

%% THIRD: plot aligned initial and final saccade endpoints
f_ts = figure;
to_plot = {'i_sacc','f_sacc'}; % what fields do we want to plot?

use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, which_excl), ii_sess.excl_trial, 'UniformOutput',false));
trial_exist = zeros(length(ii_sess.excl_trial),1)==0;
trial_exist=~cellfun(@any,cellfun(@(a) ismember(10,a),ii_sess.excl_trial,'UniformOutput',false));
use_trial = trial_exist&use_trial;
use_trial(ii_sess.i_sacc_rt < 0.1) = 0;

for pp = 1:length(to_plot)
    subplot(length(to_plot),1,pp);
    hold on;
    
    % use scatter because it gives us translucent datapoints (revert to
    % plot if necessary)
    
    % flip y for every other quadrant so that up is towards vertical
    % meridian, down is towards horizontal
    tmpy = ii_sess.(to_plot{pp})(use_trial==1,2);
    thisflipidx = sign(ii_sess.targ(use_trial==1,1))~=sign(ii_sess.targ(use_trial==1,2));
    tmpy(thisflipidx) = -1*tmpy(thisflipidx);
    
    %scatter(ii_sess.(to_plot{pp})(use_trial==1,1),ii_sess.(to_plot{pp})(use_trial==1,2),20,mycolors(pp,:),'filled','MarkerFaceAlpha',0.5);
    s = scatter(ii_sess.(to_plot{pp})(use_trial==1,1),tmpy,20,mycolors(pp,:),'filled','MarkerFaceAlpha',0.5);
    s.MarkerFaceAlpha = 0.15;
    s.MarkerFaceColor = [0 0 1];
    
    plot(0,0,'k+','MarkerSize',8,'MarkerFaceColor','k');
    plot(9,0,'kx','MarkerSize',8,'MarkerFaceColor','k');
    xlim([-2 15]);
    axis equal off;
    
    title(to_plot{pp},'Interpreter','none');
end
saveas(f_ts,sprintf('%s/%s_i&f_endpoints.png',edf_files(1).folder,taskinfo{1}.task.subjectID));
clear f_ts;
%% FOURTH: align all trial timecouress and saccade traces and plot
% sometimes it's useful to plot a stacked set of saccade timecourses
% aligned towards the direction of the target. Here, I'll do that for the
% extracted traces and the full response timecourse, as two different
% subplots.

use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, which_excl), ii_sess.excl_trial, 'UniformOutput',false));
trial_exist = zeros(length(ii_sess.excl_trial),1)==0;
trial_exist=~cellfun(@any,cellfun(@(a) ismember(10,a),ii_sess.excl_trial,'UniformOutput',false));
use_trial = trial_exist&use_trial;
use_trial(ii_sess.i_sacc_rt < 0.1) = 0;

% start with just i_sacc_trace
f_ts = figure;
subplot(3,1,[1 2]); hold on;
if strcmp(task,'wmloads_SimuFb')% only frist resp_epoch (SimuFb)
    thisidx = find(ii_sess.resp_num ==1 & use_trial==1);
else % SepFb
    thisidx = find(use_trial==1);
end

for tt = 1:length(thisidx)
    [tmpth,tmpr] = cart2pol(ii_sess.i_sacc_trace{thisidx(tt)}(:,1),ii_sess.i_sacc_trace{thisidx(tt)}(:,2));
    
    % change th, keeping r the same, based on th of ii_sess.targ
    [adjth,~] = cart2pol(ii_sess.targ(thisidx(tt),1),ii_sess.targ(thisidx(tt),2));
    
    [aligned_x,aligned_y] = pol2cart(tmpth-adjth,tmpr);
    
    % if top left or bottom right quadrant, flip y
    if sign(ii_sess.targ(thisidx(tt),1))~=sign(ii_sess.targ(thisidx(tt),2))
        aligned_y = -1*aligned_y;
    end
    
    plot(aligned_x,aligned_y,'-','LineWidth',1,'Color',[0.2 0.2 0.2]);
    clear tmpth tmpr adjth aligned_x aligned_y;
end
% fixation point
plot(0,0,'+','MarkerSize',8,'MarkerFaceColor',mycolors(2,:),'Color',mycolors(2,:),'LineWidth',2);

% target location
plot(9,0,'d','Color',mycolors(1,:),'MarkerFaceColor',mycolors(1,:),'MarkerSize',8);
xlim([-2 15]); axis equal;
set(gca,'TickDir','out','XTick',0:6:12);
title('Aligned primary saccade trajectories');

% now, the aligned trace
if strcmp(task,'targ4_control')
    xdat_to_plot = [3 4; 5 6; 7 8; 9 10]; % plot response & feedback epoch
elseif strcmp(task,'targ2_control')
    xdat_to_plot = [3 4; 6 7]; % plot response & feedback epoch
elseif strcmp(task,'targ1_control')
    xdat_to_plot = [3 4]; % plot response & feedback epoch
else % regular task
    xdat_to_plot = [4 5];
end

subplot(3,1,3); hold on;

for resp_id = 1:length(ii_params.response_epoch)
    xdat = xdat_to_plot(resp_id,:);
    thisidx = find(ii_sess.resp_num ==1 & use_trial==resp_id);
    for tt = 1:length(thisidx)
        % similar to before, but now we're going to rotate X,Y channels
        [tmpth,tmpr] = cart2pol(ii_sess.X{thisidx(tt)},ii_sess.Y{thisidx(tt)});
        
        % change th, keeping r the same, based on th of ii_sess.targ
        [adjth,~] = cart2pol(ii_sess.targ(thisidx(tt),1),ii_sess.targ(thisidx(tt),2));
        
        [aligned_x,~] = pol2cart(tmpth-adjth,tmpr); % we're not using y here, just x
        
        aligned_x = aligned_x(ismember(ii_sess.XDAT{thisidx(tt)},xdat));
        
        this_t = (1:length(aligned_x))/ii_cfg.hz;
        plot(this_t,aligned_x,'-','LineWidth',1,'Color',[0.2 0.2 0.2]);
        
    end
end
xlabel('Time (s) after GO');
ylabel('Eye position (towards target)');
xlim([0 1.5]); ylim([-2 15]);
set(gca,'XTick',[0:0.7:1.4],'YTick',[0:6:12],'TickDir','out');

saveas(f_ts,sprintf('%s/%s_timecourses_traces.png',edf_files(1).folder,taskinfo{1}.task.subjectID));
clear f_ts;

%% Further Basic Analysis
use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, which_excl), ii_sess.excl_trial, 'UniformOutput',false));
trial_exist = zeros(length(ii_sess.excl_trial),1)==0;
trial_exist=~cellfun(@any,cellfun(@(a) ismember(10,a),ii_sess.excl_trial,'UniformOutput',false));
use_trial = trial_exist&use_trial;
use_trial(ii_sess.i_sacc_rt < 0.1) = 0;

setsizes=unique(ii_sess.setsize);
for i = 1:length(setsizes)
thisidx = find(ii_sess.setsize==setsizes(i));
i_errors = [i_errors nanmean(ii_sess.i_sacc_err(thisidx))];
f_errors = [f_errors nanmean(ii_sess.f_sacc_err(thisidx))];
end
