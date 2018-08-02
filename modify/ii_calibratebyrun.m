function [ii_data,ii_cfg] = ii_calibratebyrun(ii_data,ii_cfg, chan_names, calib_targets, poly_deg, fit_limit )
%II_CALIBRATEBYRUN Adjust gaze coords across a run so that gaze is best-matched,
%on average, to a known value at a known period of trial.
%   Differs from ii_calibratebytrial: fits a polynomial individually to
%   each channel to best match observed and desired coordiante and applies
%   transform to all timepoints across all trials. Typically trial-wise
%   drift correction is done before this step.
%
%   One selection expected per trial, and either channels, column indices
%   to ii_cfg.trialinfo, or n_trials length vectors, to which the full
%   trial timeseries is warped.
%   POLY_DEG is the degree of polynomial used for readjustment. By default,
%   cubic polynomial (degree = 3) is used. 
%   FIT_LIMIT is the maximum allowable error, in dva, between the selected
%   fixation and the calib_target on that trial, in euclidean distance, for
%   a trial to be included in the FIT step (estimating the recalibration
%   function). All trials will be adjusted based on the single run-wise
%   recalibration function. By default, inf, so all trials included.
%
% Note: it's helpful to have run driftcorrect before this so that
% 'fixation' values are set to 0,0, and so scaling doesn't substantially
%  alter mean fixation
%
% Purpose of this version of calibration is to deal with somewhat
% substantial number of trials without a corrective saccade to visual
% target re-presentation during feedback, which previously resulted in
% 0-dva MGS error metrics for some subj.
%
% NOTE: does not scale velocity!!! this may take some work...
%
% Updated 6/13/2018 - now uses absolute distance threshold for deciding
% whether to calibrate or not. If distance between input fixation and
% target is > adj_limits, do not calibrate.
% - saves out calibrate.amt: amount of scaling applied (factor used for
%   dividing/multiplying) to each channel
% - also saves out calibrate.err: euclidean distance between target point
%   and fixation (1 number)
% 

% Tommy Sprague, 8/2/2018 



% make sure trials have been defined, and there's a selection for every
% trial
if ~ismember('trialvec',fieldnames(ii_cfg))
    error('iEye:ii_calibratebyrun:noTrialsDefined', 'Trials not defined; run ii_definetrials.m first');
end

% TODO:
% as we loop over trials, if a selection is missing within a trial we won't adjust and will just output a warning

if sum(ii_cfg.sel)==0
    error('iEye:ii_calibratebyrun:noSelections', 'No timepoints selected! Try ii_selectfixationsbytrial');
end

% DEFAULT VALUES

% if no channels pre-defined, use X,Y
if nargin <3 || isempty(chan_names)
    chan_names = {'X','Y'};
end

% if no calib targets defined, assume TarX, TarY (will be checked below..)
if nargin < 4 || isempty(calib_targets)
    calib_targets = {'TarX','TarY'};
end

if nargin < 5 || isempty(poly_deg)
    poly_deg = 3;
end


% CLEAN INPUTS

if ~iscell(chan_names)
    chan_names = {chan_names};
end

% do channels exist?
for cc = 1:length(chan_names)
    if ~ismember(chan_names{cc},fieldnames(ii_data))
        error('iEye:ii_calibratebyrun:invalidChannel', 'Channel %s not found',chan_names{cc});
    end
end


if ~iscell(calib_targets) && ischar(calib_targets)
    calib_targets = {calib_targets};
end

if nargin < 6 || isempty(fit_limit)
    fit_limit = inf;
end

% make sure only a single number provided and it's positve (distance)
if length(adj_limit)~=1||adj_limit<0
    error('iEye:ii_calibratebyrun:invalidAdjLimit','Adjustment limits specified in dva distance, so one single positive number required');
end

% check calibration targets 
% numeric inputs can be either indices into ii_cfg.trialinfo, or n_trials
% x n_channels values
if isnumeric(calib_targets)
    
    
    if numel(calib_targets)~=length(chan_names) && numel(calib_targets)~=(length(chan_names)*ii_cfg.numtrials)
        error('iEye:ii_calibratebyrun:invalidCalibTargets', 'For numeric calib_targets input, need either indices to ii_cfg.trialinfo columns, or one value per channel per row (n_trials x n_chans)');
    end
    
    % if there are n_chans numbers, make sure they're integers AND that
    % they're in range of ii_cfg.trialinfo AND that ii_cfg.trialinfo
    % exists...
    
    if numel(calib_targets)==length(chan_names)
        % this one requires trialinfo be present
        if ~ismember('trialinfo',fieldnames(ii_cfg))
            error('iEye:ii_calibratebyrun:invalidCalibTargets', 'For numeric calib_targets input, must have trialinfo installed (run ii_addtrialinfo)');
        end
        if max(calib_targets) > size(ii_cfg.trialinfo,2)
            error('iEye:ii_calibratebyrun:invalidCalibTargets', 'calib_targets out of range for ii_cfg.trialinfo (calib_targets: %s, trialinfo dims: %s)',num2str(calib_targets),num2str(size(ii_cfg.trialinfo)));
        end
    end
    
    if numel(calib_targets)==(length(chan_names)*ii_cfg.numtrials)
        % these are harder to check... (I think have to make sure not
        % zero?)
    end
    
% cell inputs can either be all strings, or cell of vectors each n_trials
% long
elseif iscell(calib_targets) 
    if sum(cellfun(@ischar,calib_targets))==length(calib_targets)        
        if sum(ismember(calib_targets,fieldnames(ii_data)))~=length(calib_targets)
            error('iEye:ii_calibratebyrun:invalidCalibTargets', 'For cell array of strings as calib_targets input, all cells must match field names in ii_data');
        end
    elseif sum(cellfun(@length,calib_targets)==ii_cfg.numtrials)~=length(calib_targets)        
        % error: missing 
        error('iEye:ii_calibratebyrun:invalidCalibTargets', 'For cell array of numeric calib_targets input, all cells must have numtrials elements');        
    else
        % general error
        error('iEye:ii_calibratebyrun:invalidCalibTargets', 'Invalid cell array of calib_targets');
    end   
else
    % invalid input (not numeric, or a cell of chan names)
    error('iEye:ii_calibratebyrun:invalidCalibTargets', 'Invalid input for calibration targets. Must be cell array of strings matching length of chan_names, list of trialinfo column indices matching length of chan_names, matrix with number of columns matching length of chan_names, or cell with length(chan_names) entries, each with one row per trial');
    
end


% check if .calibrate field already present, warn that you may be
% double-calibrating
if ismember('calibrate',fieldnames(ii_cfg))
    warning('iEye:ii_calibratebyrun:alreadyCalibrated','CALIBRATE field already present in ii_cfg; will overwrite!');
end

% START PROCESSING!
tu = unique(ii_cfg.trialvec(ii_cfg.trialvec~=0));

% get coords to calibrate to
calibrate_to = nan(ii_cfg.numtrials,length(chan_names));

% if cell, could either be strings or numerics
if iscell(calib_targets)
    % if each element of calib_targets is a string:
    if sum(cellfun(@ischar,calib_targets))==length(calib_targets)
        
        for cc = 1:length(calib_targets)
            for tt = 1:length(tu)
                
                tr_idx = ii_cfg.trialvec==tu(tt);
                
                calibrate_to(tt,cc) = mean(ii_data.(calib_targets{cc})(ii_cfg.sel&tr_idx));
                
            end
        end
    % if each cell is numeric & a vector (checked above)
    elseif sum(cellfun(@isnumeric,calib_targets))==length(calib_targets)
       
        for cc = 1:length(calib_targets)
            calibrate_to(:,cc) = calib_targets{cc};
        end
    end

% if there are just indices into ii_cfg.trialinfo or a matrix of coords
elseif isnumeric(calib_targets)
    % if just indices
    if numel(calib_targets)==length(chan_names)
        
        for cc = 1:length(calib_targets)
            calibrate_to(:,cc) = ii_cfg.trialinfo(:,calib_targets(cc));
        end
        
    % if all trial calibration targets
    elseif size(calib_targets,1)==ii_cfg.numtrials && size(calib_targets,2)==length(chan_names)
        calibrate_to = calib_targets;
    end
end




% ok - now we have the calibration targets numerically defined, let's fit a
% function to the gaze coordinate on each trial to these positions
