function [ii_trial,ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,targ_coords,resp_epoch,fix_epoch,excl_criteria,save_chans,score_mode,align_to)
%ii_scoreMGS Extract typical parameters from each trial of an MGS dataset:
% - primary saccade: first saccade, exceeding some amplitude threshold,
%   after a response cue in the direction of one or more target positions
% - final saccade: endpoint of final eye position before a feedback or
%   return-to-fixation stimulus appears
% - RT: time from beginning of go cue until primary saccade is detected
%
% Then, knowing these paramters & features of the trial like target
% position(s), we also save out an aligned version of all saccade endpoints
% and traces. 
%
% ii_trial contains one entry per trial of each of several fields:
% - i_sacc: n_trials x 2, ALIGNED X Y of primary saccade endpoint
% - f_sacc: n_trials x 2, ALIGNED X Y of final saccade endpoint (note:
%   final saccade & primary saccade endpoint can be same value when only a
%   single saccade is detected)
% - n_sacc: n_trials x 1, number of saccades detected after and including  
%    primary saccade (if no primary saccade, this will be 0)
% - excl_trial: structure listing several exclusion criterion applying to
%   drift correction, epochs that are considered 'fixation', calibration
%   criteria, etc. 
% - i/f_sacc_trace: cell array (n_trials x 1) with trace of initial/final
%   saccade through time, beginning at i_sacc_rt
% - i/f_sacc_err: euclidean distance between first target and endpoint of
%   initial/final saccade
% - targ: n_trials x 2, extracted (or input) target coordinates for X,Y on 
%   each trial
%
% Because ii_data, ii_sacc are not modified, we won't return them here.
%
% We extract the same fields as were extracted in ii_extractsaccades
%
% Usage:
%  [ii_trial, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc) uses default
%  saccade scoring parameters to extract primary and final eye positions
%  for a standard MGS-style task. Be default, looks for target coordinate
%  in TarX, TarY channels, extracts X,Y,Pupil and XDAT, and identifies
%  resp_epoch as first one following a the longest non-final epoch in a
%  trial (fixation epochs defined as all those previous). [TODO]
%
%  [ii_trial, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,targ_coords)
%  uses coordinates specified by targ_coords instead of default (TarX,TarY
%  channels of ii_data). targ_coords can be a cell array of 2 strings, in
%  which case those channels are looked up in ii_data, a set of n_trials x
%  2n coordiantes, in which case pairs of columns are considered X,Y pairs
%  per trial (and we determine which of the n coords was chosen), or an array
%  of 2n numbers, in which case ii_cfg.trialinfo(:,targ_coords) will act as
%  X,Y pairs)
%
%  [ii_trial, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],resp_epoch) only
%  looks at saccades beginning and ending with XDAT within resp_epoch 
%  (integer or array of integers)
%
%  [ii_trial, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],[],fix_epoch)
%  only considers fixation breaks when XDAT within fix_epoch (integer or
%  array of integers)
%
%  [ii_trial, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],[],[],excl_criteria)
%  uses values stored within excl_criteria to determine primary saccade
%  amplitude/duration threshold, error threshold, drift correction
%  threshold, and delay-period fixation break threshold (see below for
%  exact field names [TODO])
%
%  [ii_trial, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],[],[],[],save_chans)
%  saves trial timecourse for each trial from each channel listed in
%  save_chans rather than the default set ('X','Y',Pupil and XDAT; note
%  that X/Y are derived from those used when extracting saccades)
%
% [ii_trial, ii_cfg] =
% ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],[],[],[],[],score_mode) determines
% whether a scored saccade must begin AND end during indicated epoch(s)
% ('strict'), or just begin ('lenient'). Default is 'strict'.
%
% [ii_trial, ii_cfg] =
% ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],[],[],[],[],[],align_to) aligns
% i_sacc_raw and f_sacc_raw to coordinate in align_to (1x2). if absent,
% aligns only based on polar angle, DOES NOT adjust ecc!!! If input for X 
% is nan, also does not adjust ecc. (default: [NaN 0])
%
% Examples:
% TODO
%
% NOTES:
% - at present, we require at least one target coordinate to be defined; in
%   principle, that's not necessary - we can add functionality for that later
% - also, only a single target coordinate for now: ii_nearestcoord could be
%   called from here, but let's assume it's already been called
% - asusming coordinate pairs: X,Y (no support for 1D or 3D, etc, coords) -
%   this is consistent with ii_extractsaccades
% - for aligning all targets, do we want to scale in polar coords, or
%   translate in cartesian? [discuss!]
%
% TODO:
% - add a i_sacc, f_sacc idx to ii_sacc - it's zero if not a scored
%   saccade; trial number if it's either of those (useful for plotting main
%   sequences later on)
%
% TRIAL EXCLUSION CODES:
% first-digit:
% - 1 - trial-level exclusion (bad drift correction [11], calibration [12], or delay-
%       fixation break [13]
% - 2 - primary saccade exclusion (no primary sacc [20]; too small/short [21], large error [22])
%
%
% Tommy Sprague, 6/11/2018 - adapted from individual study
% _extractSaccadeData scripts (e.g., MGSMap_extractSaccadeData)

% =====================================
% Adapted for Multi-target all-report experimental setting

% Added following fields in ii_trial:
% - i/f_sacc_allerr: n_trials x 2 x 4, euclidean distance between all target and endpoint of
%   initial/final saccade， trialnum by tar_num
% - i/f_sacc_errid: n_trials x 1, the min euclidean distance is calcu from which target
% - i/f_saccall: n_trials x 2 x 4, all the results for i/f_sacc aligned to all targ
% - alltarg: n_trials x 2 x 4, extracted (or input) target coordinates for X,Y on 
%   each trial

% Revised following fields in ii_trial
% - i/f_sacc_err: n_trials x 1, euclidean distance between first target and endpoint of
%   initial/final saccade
% - i/f_sacc: sacc raw aligned to the i/f_sacc target (i/f_sacc_errid)
% - targ: n_trials x 2, get the i_sacc closest one for each epoch

% Added ii_cfg.targ_exist in each resp_epoch is trial by targ logic array,
% where 1 represents that that resp_epoch exist, 0 if not exist, if not
% exist, ii_trial.excl_trial{tt} get [10] (epoch doesnt exist), in which
% case, I dont assig this trial with any other excl_trial code

% revised excl_default.i_dur_thresh from 150 to 0.150s for unified unit

% Added 2 ii_trial.excl_trial values, 10 if no targ exist in that epoch dur
% that trial, and 23 if i and f saccade is not closet to the same target

% Revised to get the target coordinate info from that exact resp_epoch, not
% resp_epoch(end)+1

% Revised from getting targ_coords from the the epoch immediately
% after the end of the response epoch (resp_epoch(end)+1), to getting it
% from the current resp_epoch(end)

% Before running:
% Change ii_trial.n_sacc_epoch = nan(ii_cfg.numtrials,5) (line 267) but i
% found that it's not be used...) into n_trials x n_epochs (exclude ITI)

% excl_default.i_err_thresh and delay_fix_thresh from line 184


% Qingqing Yang, 6/8/2022

%% sanitize inputs & set default values

% just not enough inputs
if nargin < 3
    error('iEye:ii_scoreMGS:insufficientInputs','At minimum, to score MGS data you must submit ii_data,ii_cfg, and ii_sacc (extracted saccade data)');
end

% if targ_coords empty, use TarX, TarY
if nargin < 4 || isempty(targ_coords)
    targ_coords = {'TarX1','TarY1','TarX2','TarY2','TarX3','TarY3','TarX4','TarY4'};
end

% if no response epoch defined:
% deduce which epoch to use for response...first, we look for delay,
% which we define as the longest epoch >= 1 and <= max(epoch(within-tiral))
% (this assumes that the ITI is the last epoch...)
if nargin < 5 || isempty(resp_epoch)
    tmp_iti = max(ii_data.XDAT(ii_cfg.trialvec~=0));
    tmp_delay = mode(ii_data.XDAT(ii_cfg.trialvec~=0 & ii_data.XDAT~=tmp_iti));
    resp_epoch = tmp_delay + 1;
    clear tmp_delay tmp_iti;
end

% by default, use all epochs up until resp_epoch
if nargin < 6 || isempty(fix_epoch)
    fix_epoch = 1:(resp_epoch-1);
end

% these are default exclusion criteria, excl_criteria input only needs a
% subset of these, if a field is missing it'll be overwritten w/ a default
% value

excl_default.i_dur_thresh = 0.150; % must be shorter than 150 ms
excl_default.i_amp_thresh = 5;   % must be longer than 5 dva [if FIRST saccade in denoted epoch is not at least this long and at most this duration, drop the trial]
excl_default.i_err_thresh = 5;   % i_sacc must be within this many DVA of target pos to consider the trial
% excl_default.i_err_thresh = 10; % for debugging
excl_default.drift_thresh = 2.5;     % if drift correction norm is this big or more, drop
excl_default.delay_fix_thresh = 2.5; % if any fixation is this far from 0,0 during delay (epoch 3)
% excl_default.delay_fix_thresh = 5; % for debugging 

% if excl_criteria is empty, make it a struct and insert default into
% missing values
if nargin < 7 || isempty(excl_criteria)
    excl_criteria = struct(); % empty struct
end
    
excl_fields = fields(excl_default);
for ee = 1:length(excl_fields)
    if ~ismember(excl_fields{ee},fields(excl_criteria))
        excl_criteria.(excl_fields{ee}) = excl_default.(excl_fields{ee});
    end
end
% now, we have excl_criteria w/ all the wanted fields, and default values
% filling in those not provided originally
clear excl_default; % make sure we don't accidentally use this guy again


if nargin < 9 || isempty(score_mode)
    score_mode = 'strict';
elseif ~ismember(score_mode,{'strict','lenient'})
    error('iEye:ii_scoreMGS:invalidMode','Only strict and lenient allowed as score_mode inputs');
end

if nargin < 10 || isempty(align_to)
    align_to = [NaN 0];
end


%% Figure out which channels we're looking at. {X,Y}
% Channels to extract is defined for ii_extractsaccades, so we need to
% figure out what those fields were - can pull it out of field names in
% ii_sacc (use $$_trace) - look for whatever comes before _trace (sorry, a
% bit ugly below, but cellfun's weren't working the way I wanted)
tmp_fields = fields(ii_sacc);
tmp_idx = find(cellfun(@any,strfind(tmp_fields,'_trace')));
which_chans = cell(length(tmp_idx),1);
for cc = 1:length(tmp_idx)
    which_chans{cc} = tmp_fields{tmp_idx(cc)}(1:(strfind(tmp_fields{tmp_idx(cc)},'_trace')-1));
end
clear tmp_fields tmp_idx;


if nargin < 8 || isempty(save_chans)
    save_chans = {which_chans{:},'Pupil','XDAT','Velocity'};
end


%% NOW: start extracting data!

%% Build ii_trial; fill with nothing % {X,Y,Pupil,XDAT,Velocity}
for chan_idx = 1:length(save_chans)
    ii_trial.(save_chans{chan_idx}) = cell(ii_cfg.numtrials,1);
end

% also want to store the raw coordinates
ii_trial.i_sacc_raw = nan(ii_cfg.numtrials,2);
ii_trial.f_sacc_raw = nan(ii_cfg.numtrials,2);

% aligned coordinates
ii_trial.i_sacc = nan(ii_cfg.numtrials,2);
ii_trial.f_sacc = nan(ii_cfg.numtrials,2);

ii_trial.i_saccall = nan(ii_cfg.numtrials,length(targ_coords));
ii_trial.f_saccall = nan(ii_cfg.numtrials,length(targ_coords));

ii_trial.i_sacc_err = nan(ii_cfg.numtrials,1);
ii_trial.f_sacc_err = nan(ii_cfg.numtrials,1);

ii_trial.i_sacc_errid = nan(ii_cfg.numtrials,1);
ii_trial.f_sacc_errid = nan(ii_cfg.numtrials,1);

ii_trial.i_sacc_allerr = nan(ii_cfg.numtrials,length(targ_coords)/2);
ii_trial.f_sacc_allerr = nan(ii_cfg.numtrials,length(targ_coords)/2);

ii_trial.n_sacc = nan(ii_cfg.numtrials,1); % how many saccades are there total? (in each epoch maybe?)
ii_trial.n_sacc_epoch = nan(ii_cfg.numtrials,5); % TODO: fill this w/ n_trials x n_epochs (exclude ITI)
% ii_trial.n_sacc_epoch = nan(ii_cfg.numtrials,8)

ii_trial.i_sacc_rt = nan(ii_cfg.numtrials,1); % latency from go cue to each of these
ii_trial.f_sacc_rt = nan(ii_cfg.numtrials,1);

ii_trial.i_sacc_trace = cell(ii_cfg.numtrials,1);
ii_trial.f_sacc_trace = cell(ii_cfg.numtrials,1);

ii_trial.i_sacc_peakvel = nan(ii_cfg.numtrials,1);
ii_trial.f_sacc_peakvel = nan(ii_cfg.numtrials,1);

% save some calibration, drift correction info for convenience
if isfield(ii_cfg,'calibrate')
    ii_trial.calib_amt = ii_cfg.calibrate.amt;
    ii_trial.calib_adj = ii_cfg.calibrate.adj;
    ii_trial.calib_err = ii_cfg.calibrate.err;
end

if isfield(ii_cfg,'drift')
    ii_trial.drift_amt = ii_cfg.drift.amt;
end

ii_trial.excl_trial = cell(ii_cfg.numtrials,1);  % why is this trial excluded? each cell includes several markers


% let's copy over ii_cfg.trialinfo, if it exists
if isfield(ii_cfg,'trialinfo')
    ii_trial.trialinfo = ii_cfg.trialinfo;
end

% add parameters used for extracting saccades: ii_trial.params (as they
% were input/sanitized - but, for e.g., targ_coord, not updated [?])
ii_trial.params.excl_criteria = excl_criteria;
ii_trial.params.resp_epoch  = resp_epoch;
ii_trial.params.fix_epoch   = fix_epoch;
ii_trial.params.targ_coords = targ_coords;
ii_trial.params.save_chans  = save_chans;
ii_trial.params.score_mode  = score_mode;
ii_trial.params.score_chans = which_chans; % the two channels that were used for scoring


%% extract targ_coords for each trial

targ_coords_extracted = nan(ii_cfg.numtrials,8);

% if targ_coords is n_trials x 2, then just extract from there
if size(targ_coords,1)==ii_cfg.numtrials && size(targ_coords,2) == 2
    targ_coords_extracted = targ_coords;

% if a cell, then we need to look in ii_data at these fields for each trial
% we assume that the target coordinate is given by the epoch immediately
% after the end of the response epoch (resp_epoch(end)+1)_TS
% Here we get the targ_coord from that exact resp_epoch
% always takes the mode - so uses whatever is present during resp_epoch of
% this channel most commonly (to account for mis-sampling, etc)
elseif iscell(targ_coords)
    for tt = 1:ii_cfg.numtrials % for each trial, get targ_coords
        thisidx = ii_cfg.trialvec==tt & ismember(ii_data.XDAT,resp_epoch);
        for cc = 1:length(targ_coords)
            targ_coords_extracted(tt,cc) = mode(ii_data.(targ_coords{cc})(thisidx));
        end
    end
    
% if 2 #'s, pull info out of those columns ii_cfg.trialinfo
elseif numel(targ_coords)==2
    if max(targ_coords) > size(ii_cfg.trialinfo,2)
        error('iEye:ii_scoreMGS:invalidInput','Column index submitted for targ_coords exceeds size of ii_cfg.trialinfo (%i columns)',size(ii_cfg.trialinfo,2));
    end
    targ_coords_extracted = ii_cfg.trialinfo(:,targ_coords);
else
    error('iEye:ii_scoreMGS:invalidInput','Invalid input for targ_coords: input should be a pair of channels (cell array of strings), list of coordinates (numtrials x 2), or pair of columns into ii_cfg.trialinfo');
end

% save these in ii_trial as ii_trial.alltarg
% targ_coords_extracted will still be used in this script
ii_trial.alltarg = targ_coords_extracted;
% the ii_trial.targ will be extracted by the errid from primary sacc


% get idx saccades to look at (those that start [end] in resp_epoch), and
% corresponding trial index
if strcmpi(score_mode,'strict')
    which_sacc = find(ismember(ii_sacc.epoch_start,resp_epoch) & ismember(ii_sacc.epoch_end,resp_epoch));
elseif strcmpi(score_mode,'lenient')
    % if lenient, just needs to start in this epoch
    which_sacc = find(ismember(ii_sacc.epoch_start,resp_epoch));
end

which_trials = ii_sacc.trial_start(which_sacc); % trial number for each saccade in which_sacc

%% loop over trials
for tt = 1:ii_cfg.numtrials
    
    % save the data from each channel from each trial
    % ii_trial.(save_chans){i,1} is ii_data.(save_chans) whole trial
    for chan_idx = 1:length(save_chans)
        ii_trial.(save_chans{chan_idx}){tt} = ii_data.(save_chans{chan_idx})(ii_cfg.trialvec==tt);
    end
    
    % time that relevant epoch of trial started 
    t_start = find(ii_cfg.trialvec==tt & ismember(ii_data.XDAT,resp_epoch) ,1,'first')/ii_cfg.hz; 
    
    % init the test if the epoch in this trial exist
    targ_id=find(resp_epoch==ii_cfg.params.response_epoch);
    ii_cfg.targ_exist=zeros(ii_cfg.numtrials,length(ii_cfg.params.response_epoch))==0;
    %% score saccades: find primary/initial; final saccade
    
    % all saccades that started (and ended?) in resp_epoch on this trial
    this_i_sacc = which_sacc(which_trials==tt);
    
    % extract amplitude of each of these
    this_i_amp = ii_sacc.amplitude(this_i_sacc);
    this_i_dur = ii_sacc.duration(this_i_sacc);
    
    
    % if amplitude too small or duration too short (for ALL
    % saccades in response epoch
    
    % first, if no saccade detected (if i_sacc is empty),
    % record as "20" - no identified saccades (whatosever!)
    
    if isempty(t_start)
        ii_cfg.targ_exist(tt,targ_id)= 0==1; % if there's no that targ in this trial, targ_exist=0;
        ii_trial.excl_trial{tt}(end+1) = 10; % if there's no targ
    else % so this trial has this epoch (target)
        ii_cfg.targ_exist(tt,targ_id)= 1==1;
        if isempty(this_i_sacc)
            ii_trial.excl_trial{tt}(end+1) = 20; % no saccades identified in response epoch
        elseif ~any(this_i_dur<=excl_criteria.i_dur_thresh & this_i_amp>=excl_criteria.i_amp_thresh)  % if none of the saccades pass the amplitude & duration criteria
            ii_trial.excl_trial{tt}(end+1) = 21; % none of the identified saccades (>0) passed exclusion criteria for primary saccade
        else
            % ok, now we can use the first of the saccades that do
            % pass exclusion criteria as the 'primary' saccade

            % find the first of this_i_sacc for which amp & dur
            % pass threshold
            this_i_sacc = this_i_sacc(this_i_amp>=excl_criteria.i_amp_thresh & this_i_dur <= excl_criteria.i_dur_thresh);

            % just index into the first element of this_i_sacc... [TODO: use
            % which_chan to compute this....]
            ii_trial.i_sacc_raw(tt,:) = [ii_sacc.X_end(this_i_sacc(1)) ii_sacc.Y_end(this_i_sacc(1))];
            ii_trial.i_sacc_rt(tt) = ii_sacc.t(this_i_sacc(1),1)-t_start; % ii_sacc.t take 1st column, so start time of sacc, and t_start is the time that this epoch starts in this trial
            ii_trial.i_sacc_trace{tt} = [ii_sacc.X_trace{this_i_sacc(1)} ii_sacc.Y_trace{this_i_sacc(1)}];

            % get peak velocity 

            ii_trial.i_sacc_peakvel(tt) = ii_sacc.peakvelocity(this_i_sacc(1)); 

        end
    clear this_i_amp this_i_dur;
    
    end
    % same for final sacc % why not testing dur and amp for final sacc?
    this_f_sacc = which_sacc(find(which_trials==tt,1,'last')); %last sacc start in this epoch, in this trial
    if ~isempty(this_f_sacc)
        ii_trial.f_sacc_raw(tt,:) = [ii_sacc.X_end(this_f_sacc) ii_sacc.Y_end(this_f_sacc)];
        ii_trial.f_sacc_rt(tt) = ii_sacc.t(this_f_sacc,1)-t_start; % start time for this sacc - epoch start
        ii_trial.f_sacc_trace{tt} = [ii_sacc.X_trace{this_f_sacc} ii_sacc.Y_trace{this_f_sacc}];
        ii_trial.f_sacc_peakvel(tt) = ii_sacc.peakvelocity(this_f_sacc(1)); 
    end
    
    % count saccades in this trial during relevant epoch AFTER
    % primary saccade, up to and including final saccade
    
    % first saccade: primary
    
    % last saccade: final
    
    % so, I think we just need the length of i_sacc:f_sacc
    if ~isempty(this_i_sacc) && ~isempty(this_f_sacc)
        % in most cases, we'll identify both i_sacc & f_sacc -
        % use this algo
        ii_trial.n_sacc(tt) = length(this_i_sacc:this_f_sacc);
    else
        % if either of them is undefined (or both), this will
        % give us the # of saccades: 1 if only i_sacc or only
        % f_sacc, 0 if neither
        ii_trial.n_sacc(tt) = sum([~isempty(this_i_sacc) ~isempty(this_f_sacc)]);
    end

    %% trial exclusions: find reasons we may want to exclude each trial
    
    % ~~~~~ FIRST: exclude based on trial-level features (see above)
    
    % note: only reject trials based on these criteria if those steps were
    % run!
    if isempty(t_start)
        ii_cfg.targ_exist(tt,targ_id)= 0==1; % if there's no that targ in this trial, targ_exist=0;
%         ii_trial.excl_trial{tt}(end+1) = 10; % if there's no targ
    else
        ii_cfg.targ_exist(tt,targ_id)= 1==1;
    
        % DRIFT CORRECTION TOO BIG % didnt do anything if not 
        if isfield(ii_cfg,'drift')
            if sqrt(sum(ii_cfg.drift.amt(tt,:).^2)) > excl_criteria.drift_thresh
                ii_trial.excl_trial{tt}(end+1) = 11;
            end
        end

        % CALIBRATION OUTSIDE OF RANGE
        if isfield(ii_cfg,'calibrate') % adj, whether adjustment was performed, good if adj = 1, adjusted, w/in distance
            if ii_cfg.calibrate.adj(tt)~=1
                ii_trial.excl_trial{tt}(end+1) = 12;
            end
        end

        % DURING DELAY, FIXATION OUTSIDE OF RANGE

        % find fixations in this trial; epoch [TODO: make sure I'm using the
        % right channels here!!]
        this_fix_idx = ii_cfg.trialvec==tt & ismember(ii_data.XDAT,fix_epoch); % that trail, before resp epoch
        if max(sqrt(ii_data.X_fix(this_fix_idx).^2+ii_data.Y_fix(this_fix_idx).^2)) > excl_criteria.delay_fix_thresh
            ii_trial.excl_trial{tt}(end+1) = 13;
        end
    end % end of empty t_start
    
    % in each trial, compute error for initial, final saccade based on extracted targ coords
        % error for initial (note: mean absolute error..)

        % I_SACC ERROR TOO HIGH (use targ_coords_extracted)
        % Get the min error for all targs, and compare it with i_err_thresh
            if ~isempty(this_i_sacc)
                [trial,targs]=size(targ_coords_extracted);
                for pp = 1: (targs/2)
                    if targ_coords_extracted(tt,pp*2-1) ~=0
                        ii_trial.i_sacc_allerr(tt,pp)=sqrt(sum((ii_trial.i_sacc_raw(tt,1)-targ_coords_extracted(tt,pp*2-1)).^2+(ii_trial.i_sacc_raw(tt,2)-targ_coords_extracted(tt,pp*2)).^2));
                    else
                        ii_trial.i_sacc_allerr(tt,pp)=nan;
                    end
                    [minerr,minerr_id]=min(ii_trial.i_sacc_allerr(tt,:));
                    ii_trial.i_sacc_err(tt,1)=minerr; % if nan, maybe no sacc pass the thershold
                    ii_trial.i_sacc_errid(tt,1)=minerr_id;
                end
                    clear minerr minerr_id pp targs trial;
                
                if  ~isempty(t_start)&&ii_trial.i_sacc_err(tt,1)> excl_criteria.i_err_thresh
                    ii_trial.excl_trial{tt}(end+1) = 22;
                end
            end

        % error for final
        % ii_trial.f_sacc_err = sqrt(sum((ii_trial.f_sacc_raw-targ_coords_extracted).^2,2));
        [trial,targs]=size(targ_coords_extracted);
        for pp = 1: (targs/2)
            if targ_coords_extracted(tt,pp*2-1) ~=0
                ii_trial.f_sacc_allerr(tt,pp)=sqrt(sum((ii_trial.f_sacc_raw(tt,1)-targ_coords_extracted(tt,pp*2-1)).^2+(ii_trial.f_sacc_raw(tt,2)-targ_coords_extracted(tt,pp*2)).^2));
            else
                ii_trial.f_sacc_allerr(tt,pp)=nan;
            end
            [minerr,minerr_id]=min(ii_trial.f_sacc_allerr(tt,:));
            ii_trial.f_sacc_err(tt,1)=minerr;
            ii_trial.f_sacc_errid(tt,1)=minerr_id;
        end
        clear minerr minerr_id pp targs trial;

        % Check if the targid for primary and final sacc is the same
        if ~isempty(t_start)&& ii_trial.i_sacc_errid(tt,1)~= ii_trial.f_sacc_errid(tt,1)
            ii_trial.excl_trial{tt}(end+1) = 23;
        end

        clear this_i_sacc this_f_sacc pp;

    
    
end % end of trial
clear tt;

%% compute error for initial, final saccade based on extracted targ coords
    % error for initial (note: mean absolute error..)
    % ii_trial.i_sacc_err = sqrt(sum((ii_trial.i_sacc_raw-targ_coords_extracted).^2,2));

    % error for final
    % ii_trial.f_sacc_err = sqrt(sum((ii_trial.f_sacc_raw-targ_coords_extracted).^2,2));

%% rotate all trials such that i_sacc, f_sacc aligned to targ

align_fields = {'i_sacc','f_sacc'}; % use these_raw, align them

for aa = 1:length(align_fields)
    % Cartesian coordinates X,Y to polar coordinates
    % rotate all so that i_sacc, f_sacc are at known position... (align_to) 
    [trial,targs]=size(targ_coords_extracted);
    adjth = nan(ii_cfg.numtrials,4);
    adjr=nan(ii_cfg.numtrials,4);
    for pp = 1:targs/2
        % angle, radius
        [tmpth,tmpr] = cart2pol(ii_trial.(sprintf('%s_raw',align_fields{aa}))(:,1),  ii_trial.(sprintf('%s_raw',align_fields{aa}))(:,2));
        [adjth(:,pp),adjr(:,pp)] = cart2pol(targ_coords_extracted(:,pp*2-1),targ_coords_extracted(:,pp*2)); % for now, we don't use radius here for anything...
        tmpth = tmpth - adjth(:,pp);  % TODO: use align_to?
        [ii_trial.(sprintf('%sall',align_fields{aa}))(:,pp*2-1),ii_trial.(sprintf('%sall',align_fields{aa}))(:,pp*2)] = pol2cart(tmpth,tmpr);
        % adjust x if align_to is not nan [NOTE: this is a translation, not a
        % scaling]
        if ~isnan(align_to(1))
            ii_trial.(sprintf('%sall',align_fields{aa}))(:,pp*2-1) = ii_trial.(sprintf('%sall',align_fields{aa}))(:,pp*2-1) - (adjr(:,pp)-align_to(1));
        end
    end    
    clear tmpth tmpr adjth adjr;
end

% Get the i/f_saccall into i/f_sacc chosen by target (errid)
align_fields = {'i_sacc','f_sacc'}; % use these_raw, align them
i_sacc_used = ii_trial.i_sacc_errid;
f_sacc_used = ii_trial.f_sacc_errid;
for tt = 1:ii_cfg.numtrials
    for aa = 1:length(align_fields)
        if aa ==1
            if ~isnan(i_sacc_used(tt))
                ii_trial.(align_fields{aa})(tt,1)= ii_trial.(sprintf('%sall',align_fields{aa}))(tt,i_sacc_used(tt)*2-1);
                ii_trial.(align_fields{aa})(tt,2)= ii_trial.(sprintf('%sall',align_fields{aa}))(tt,i_sacc_used(tt)*2);
            elseif isnan(i_sacc_used(tt))
                ii_trial.(align_fields{aa})(tt,1)= nan;
                ii_trial.(align_fields{aa})(tt,2)= nan;
            end
        elseif aa ==2
            if ~isnan(f_sacc_used(tt))
                ii_trial.(align_fields{aa})(tt,1)= ii_trial.(sprintf('%sall',align_fields{aa}))(tt,f_sacc_used(tt)*2-1);
                ii_trial.(align_fields{aa})(tt,2)= ii_trial.(sprintf('%sall',align_fields{aa}))(tt,f_sacc_used(tt)*2);
            else 
                ii_trial.(align_fields{aa})(tt,1)= nan;
                ii_trial.(align_fields{aa})(tt,2)= nan;
            end
        end 
    end
end

%% extract the targets for this epoch in each trial
% define the targ by primary saccade endpoint here
i_sacc_used = ii_trial.i_sacc_errid;
for tt = 1:ii_cfg.numtrials
    if ~isnan(i_sacc_used(tt))
        ii_trial.targ(tt,1:2)=ii_trial.alltarg(tt,i_sacc_used(tt)*2-1:i_sacc_used(tt)*2);
    else % if no avialable sacc, default send targ 1 (excl in the future anyway)
        ii_trial.targ(tt,1:2)=ii_trial.alltarg(tt,1:2);
    end
end

return