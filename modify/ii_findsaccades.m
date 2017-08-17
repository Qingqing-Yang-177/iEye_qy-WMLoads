function [ii_data,ii_cfg] = ii_findsaccades(ii_data,ii_cfg,xchan,ychan,vel_thresh,dur_thresh,amp_thresh)
%Saccade detection
%   This function will detect and select saccades based on a particular set
%   of criteria (velocity and sample/time length).


% if nargin ~= 4
%     prompt = {'X Channel', 'Y Channel', 'Velocity Threshold', 'Length Threshold'};
%     dlg_title = 'Saccade Finder';
%     num_lines = 1;
%     answer = inputdlg(prompt,dlg_title,num_lines);
%
%     x = answer{1};
%     y = answer{2};
%     t = str2num(answer{3});
%     l = str2num(answer{4});
% end

% updated TCS 8/14/2017 - requires ii_velocity already been run, uses
% xchan/ychan only for ensuring sufficient distance between
% startpoint/endpoint
% for now, runs on entire time series - does not limit to particular
% epochs. at end of this function, saccades will be 'selected' by iEye -
% then, a separate function can be used to identify fixation periods
% between eye movements, and to extract information/metrics about the
% saccades themselves. [at least for now]. ii_cfg will have relevant info.

%%%%%%%%%% ADD DISTANCE BETWEEN!!!!

% make sure relevant channels exist
if ~ismember(xchan,fieldnames(ii_data))
    error('iEye:ii_findsaccades:channelNotFound', 'Channel %s does not exist in ii_data',xchan)
end

if ~ismember(ychan,fieldnames(ii_data))
    error('iEye:ii_findsaccades:channelNotFound', 'Channel %s does not exist in ii_data',ychan)
end

if ~ismember('velocity',fieldnames(ii_cfg))
    error('iEye:ii_findsaccades:velocityNotComputed', 'Velocity has not yet been computed. Use ii_velocity before running ii_findsaccades.')
end

% by this point, all variables accounted for...

[ii_data,ii_cfg] = ii_selectempty(ii_data,ii_cfg);

% FIND SACCADES >= T

[ii_data,ii_cfg] = ii_selectbyvalue(ii_data,ii_cfg,ii_cfg.velocity,'greaterthanequalto',vel_thresh);




% compute saccade duration
sacc_dur = (ii_cfg.cursel(:,2)-ii_cfg.cursel(:,1))/ii_cfg.hz;



% compute saccade amplitude
sacc_amp = sqrt((ii_data.(xchan)(ii_cfg.cursel(:,1))-ii_data.(xchan)(ii_cfg.cursel(:,2))).^2 + (ii_data.(ychan)(ii_cfg.cursel(:,1))-ii_data.(ychan)(ii_cfg.cursel(:,2))).^2);


sacc_keep = sacc_dur >= dur_thresh & sacc_amp >= amp_thresh;

ii_cfg.cursel = ii_cfg.cursel(sacc_keep,:);
ii_cfg.sel = 0*ii_cfg.sel;
for ii = 1:size(ii_cfg.cursel,1)
    ii_cfg.sel(ii_cfg.cursel(ii,1):ii_cfg.cursel(ii,2)) = 1;
end

ii_cfg.saccades = ii_cfg.cursel;


ii_cfg.history{end+1} = sprintf('ii_findsaccades - dur, vel, amp thresh: %d, %d, %d, chans %s, %s - %s',dur_thresh, vel_thresh, amp_thresh, xchan, ychan,datestr(now,30));



end

