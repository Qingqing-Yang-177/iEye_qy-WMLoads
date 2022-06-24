function [ii_data,ii_cfg,ii_sacc] = ii_preproc_WMLoads_SimulFb(edf_fn,cfg_fn,preproc_fn,ii_params,trialinfo,skip_steps)
% ii_preproc Performs default pre-processing stream
%
% For preprocessing saccade data (specifically, memory-guided saccade
% data), we do:
% 01. setup iEye and import our file
% 02. rescale the data given our known screen resolution, viewing distance,
% etc (pixels per degree is required)
% 03. invert the y-channel, which is measured in screen coords (+ is down)
% 04. identify and remove blinks
% 04b clean gaze channels if bad values found
% 05. use the 'XDAT' (epoch) channel to define trials
% 06. compute smoothed versions of gaze channels using a Gaussian kernel
% 07. use smoothed gaze channels to compute velocity
% 08. use velocity channel to find saccades, and filter those saccades based
% on distance, velocity, and duration thresholds
% 09. identify periods of stable fixation between saccades
% 10. select stable fixations during pre-target epoch (used to drift
% correct)
% 11. use stable fixations to drift correct to known fixation coord (0,0)
% to account for signal drift over runs in eyetracking data
% 12. add behavioral data to iEye data struct (ii_cfg)
% 13. select stable fixation at end of feedback stimulus
% 14. use this fixation to further calibrate gaze data to known target
% position. 
% 15. plot all individual trials for QC purposes, save a figure
% 16. save the data (will contain ii_data, ii_cfg - the important variables
% containing full preprocessed time-series with appropriate labels, etc)
%
% At the end of this script, ii_data and ii_cfg will be ready to use for
% scoring saccades. User will still need to identify which saccades are
% relevant (using other iEye functions), and compare those to known
% positions (note that, with ii_calibratebytrial, this is much easier).
% Several new channels are created: X_smooth, Y_smooth, X_fix, Y_fix, as
% well as a bunch of 'selection' channels in ii_cfg (trialvec, tsel,
% saccades, blinks). This pair of data structures can/should be used with
% the new ii_plot* functions: ii_plottrial, ii_plotalltrials, and
% ii_plotalltrials2d (among others on their way). 
%
% skip_steps: for now, only optional steps are drift-correction ('drift') and
% calibration ('calibration') - when it comes time to do those steps, don't
% execute them if ismember(skip_steps,those)
%
% Tommy Sprague, 8/20/2017

% =================
% Adapted for Multi-target all report, simultaneous feedback experiment
% settings.

% Revised for invert y chan, after control exp, I found out if we directly
% send dva=pixel/ppd for target y loc, we dont needa invert y (or here we
% also inverted the TarY etc. make sure to check the plot after invert for
% ur own purpose.

% Added another function for importting edf files. The alternative function
% ii_import_noedf2asc_edf.m is more compatible than ii_import_edf.m, and it
% gives out an output whatever edf input is, without the help of edf2asc.
% However, for more standard process, we still use ii_import_edf.m here.

% Added another ii_censorchans.m to remove the extreme-valued X,Y outside
% of the stimuli.aperture_size

% Revised ii_calibratebytrial.m into ii_calibratebytrial_WMLoads.m for
% calibrate during the feedback, by the first fixation and select the
% appropriate target

% For debugging, minus 1 trail within ii_selectfixationsbytrial_WMLoads
% which is from ii_selectfixationsbytrial.m (no diff, only trialnums) and 
% ii_calibratebytrial_WMLoads.m

% Qingqing Yang, 6/14/2022

% ===============================
% TODO: params contains resolution, ppd (or screen params), saccade
% detection params, etc.
% ii_params.aperture_size
% ii_selectfixationsbytrial whether is in debugging or not
% ii_params.calibrate_limits debugging or not

% initialize iEye - make sure paths are correct, etc


ii_init_qy; % use if ii_import_edf.m, remember to revise the path inside.

if nargin < 4 || isempty(ii_params)
    ii_params = ii_loadparams;
end

% in case a params file is given
if ischar(ii_params)
    ii_params = ii_loadparams(ii_params);
end

if nargin < 6 || isempty(skip_steps)
    skip_steps = {}; % don't skip anything!
end


% import data
[ii_data,ii_cfg] = ii_import_edf(edf_fn,cfg_fn,[edf_fn(1:end-4) '_iEye.mat']);
 %[ii_data,ii_cfg] = ii_import_noedf2asc_edf(edf_fn,cfg_fn,[edf_fn(1:end-4) '_iEye.mat']);

% imported_plot = plot_data(ii_data,{'X','Y','TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'})


% truncate data to relevant XDATs
[ii_data,ii_cfg] = ii_trim(ii_data,ii_cfg,ii_params.valid_epochs,ii_params.epoch_chan);
% trim_plot=plot_data(ii_data,{'X','Y','TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'})


% rescale X, Y based on screen info (pixel to dva)
[ii_data,ii_cfg] = ii_rescale(ii_data,ii_cfg,{'X','Y'},ii_params.resolution,ii_params.ppd);
% rescale_plot = plot_data(ii_data,{'X','Y','TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'})


% % Invert Y channel (the eye-tracker spits out flipped Y values)
if strcmp(ii_params.task,'targ4_control1')
    % this exp has all y already inverted, so do nothing here
else
[ii_data,ii_cfg] = ii_invert(ii_data,ii_cfg,'Y');
end
% invert_plot=plot_data(ii_data,{'X','Y','TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'})


% remove extreme-valued X,Y channels (further than the screen edges)
[ii_data,ii_cfg] = ii_censorchans(ii_data,ii_cfg,{'X','Y'},...
    {[-0.5 0.5]*ii_params.resolution(1)/ii_params.ppd,...
     [-0.5 0.5]*ii_params.resolution(2)/ii_params.ppd});
% censor_plot1 = plot_data(ii_data,{'X','Y','TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'})


% remove the X,Y which are outside of the stimuli.aperture_size
[ii_data,ii_cfg] = ii_censorchans(ii_data,ii_cfg,{'X','Y'},...
     {[-1 1]*ii_params.aperture_size,[-1 1]*ii_params.aperture_size});
% censor_plot2 = plot_data(ii_data,{'X','Y','TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'})


% Correct for blinks
[ii_data, ii_cfg] = ii_blinkcorrect(ii_data,ii_cfg,{'X','Y'},'Pupil', ...
      ii_params.blink_thresh,ii_params.blink_window(1),ii_params.blink_window(2),'prctile'); 
% blink_plot = plot_data(ii_data,{'X','Y','TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'})


% split into individual trials (so that individual-trial corrections can be
% applied)
[ii_data,ii_cfg] = ii_definetrial(ii_data,ii_cfg,...
    ii_params.trial_start_chan, ii_params.trial_start_value,...
    ii_params.trial_end_chan,   ii_params.trial_end_value); 


% Smooth data
[ii_data,ii_cfg] = ii_smooth(ii_data,ii_cfg,{'X','Y'},ii_params.smooth_type,...
    ii_params.smooth_amt);

% smooth_plot = plot_data(ii_data,{'X','Y','X_smooth','Y_smooth'})


% compute velocity using the smoothed data
[ii_data,ii_cfg] = ii_velocity(ii_data,ii_cfg,'X_smooth','Y_smooth');


% below - need to document still....

% look for saccades [[NOTE: potentially do this twice? once for macro, once
% for micro?]]
[ii_data,ii_cfg] = ii_findsaccades(ii_data,ii_cfg,'X_smooth','Y_smooth',...
    ii_params.sacc_velocity_thresh,ii_params.sacc_duration_thresh,...
    ii_params.sacc_amplitude_thresh); 

% sacc_plot = ii_plottimeseries(ii_data,ii_cfg,{'X_smooth','Y_smooth'})


% look for MICROsaccades %update to the above, 7/11/2019, geh
[ii_data,ii_cfg] = ii_findmicrosaccades(ii_data,ii_cfg,'X_smooth','Y_smooth',...
    ii_params.sacc_velocity_thresh,ii_params.sacc_duration_thresh,...
    ii_params.microsacc_amplitude_limit); %keep dur and vel thresh, change amplitude thresh 0.5-1.5 


% find fixation epochs (between saccades and blinks)
% [create X_fix, Y_fix channels? these could be overlaid with 'raw' data as
% 'stable' eye positions
[ii_data,ii_cfg] = ii_findfixations(ii_data,ii_cfg,{'X','Y'},ii_params.fixation_mode);

% fix_plot = ii_plottimeseries(ii_data,ii_cfg,{'X','Y','X_fix','Y_fix','TarX','TarY'},'noselections')
% legend({'X','Y','X fix','Y fix'})


% select the fixation used for drift correction
if ~ismember('drift',skip_steps)

    [ii_data,ii_cfg] = ii_selectfixationsbytrial( ii_data, ii_cfg, ii_params.epoch_chan,...
        ii_params.drift_epoch, ii_params.drift_fixation_mode );


% use those selections to drift-correct each trial (either using the
% fixation value, which may include timepoints past end of epoch, or using
% 'raw' or 'smoothed' data on each channel)
    [ii_data,ii_cfg] = ii_driftcorrect(ii_data,ii_cfg,{'X','Y'},ii_params.drift_select_mode,...
        ii_params.drift_target);
end


% add trial info [not explicitly necessary if simple experiment, TarX &
% TarY fields used as before]
%
% trial_info should be n_trials x n_params/features - can be indexed in
% some data processing commands below. can be cell or array. 


% if trialinfo is defined, otherwise skip
if nargin>=5 && ~isempty(trialinfo)
    [ii_data,ii_cfg] = ii_addtrialinfo(ii_data,ii_cfg,trialinfo);
end


% 'target correct' or 'calibrate' (which name is better?)
% adjust timeseries on a trial so that gaze at specified epoch (quantified
% w/ specified method, like driftcorrect) is at known coords (either given
% by a channel or a pair of cols in ii_cfg.trialinfo

% first, select relevant epochs (should be one selection per trial)

if ~ismember('calibration',skip_steps)
    % FOR SPECIFIED TARGETS/NEAREST FIXATION SELECTION
    % For debugging, keep ii_selectfixationsbytrial (line79) 
    % tu = tu(1:length(tu)-1), since the participant doesnt finish the
    % trial. For formal test, comment the line down
    % we only need to have the feedback epoch has TarX1,...,TarX4, and
    % TarY1,...,TarY4
    [ii_data,ii_cfg] = ii_selectfixationsbytrial_WMLoads(ii_data,ii_cfg,ii_params.epoch_chan,...
        ii_params.calibrate_epoch,ii_params.calibrate_select_mode,ii_params.calibrate_window,...
        ii_params.calibrate_target,{'X_fix','Y_fix'}); %
    
    if ismember(ii_params.calibrate_mode,{'scale','rotate'})
        % then, calibrate by trial, update {X,Y}
        % for debugging make ii_params.calibrate_limits bigger (10 for now,
        % 2.5 for formal test)
        % line 156: tu = tu(1:length(tu)-1); line 158 numtrials-1 % For debugging
        [ii_data,ii_cfg] = ii_calibratebytrial_WMLoads(ii_data,ii_cfg,{'X','Y'},...
            ii_params.calibrate_target,ii_params.calibrate_mode, ii_params.calibrate_limits);
    elseif ismember(ii_params.calibrate_mode,'run')
        % calibrate by run
        [ii_data,ii_cfg] = ii_calibratebyrun(ii_data,ii_cfg,{'X','Y'},...
            ii_params.calibrate_target,ii_params.calibrate_deg,ii_params.calibrate_limits);
    end
end

% plot all these - make sure they're good
f_all = ii_plotalltrials_WMLoads(ii_data,ii_cfg,{'X','Y'},[],[],ii_params.epoch_chan,ii_params.plot_epoch,ii_params.show_plots);

% save the figures for our records
if length(f_all)>1
    for ff = 1:lenght(f_all)
        saveas(f_all,sprintf('%s_%02.f.png',preproc_fn(1:end-4),ff),'png');
    end
else
    saveas(f_all,sprintf('%s.png',preproc_fn(1:end-4)),'png');
end

% and plot the final full timeseries
if ii_params.show_plots == 1
    f_ts = ii_plottimeseries(ii_data,ii_cfg,{'X','Y','PupilZ'});
else
    f_ts = ii_plottimeseries(ii_data,ii_cfg,{'X','Y','PupilZ'},'nofigure');
end
saveas(f_ts,sprintf('%s_timeseries.png',preproc_fn(1:end-4)),'png');


% save the preprocessing params
ii_cfg.params = ii_params;

% save this script?

% save ii_data,ii_cfg in _preproc.mat file
ii_savedata(ii_data,ii_cfg,preproc_fn);


% get the saccades
[ii_data,ii_cfg,ii_sacc] = ii_extractsaccades(ii_data,ii_cfg,{'X','Y'},...
    ii_params.extract_sacc_mode_start,...
    ii_params.extract_sacc_mode_end,...
    ii_params.epoch_chan);

% save the saccades
ii_savesacc(ii_cfg,ii_sacc,[preproc_fn(1:end-4) '_sacc.mat']);

end
