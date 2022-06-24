function [ii_data2,ii_cfg2] = ii_definetrial_respepoch(ii_data,ii_cfg, c1, v1, c2, v2,center_epoch,resp_epoch)
%II_DEFINETRIAL Splits data into individual trials based on signal in the
%               ii_data.(chan) variable
%   [ii_data,ii_cfg] = ii_definetrial(ii_data,ii_cfg,c1,v1,c2,v2) splits
%   trials based on values in c1 and c2: trials are defined as beginning
%   when ii_data.(c1) is v1, and ending when ii_data.(c2) is v2.
%
% c1, c2 are strings, must be fields of ii_data
% v1, v2 are numeric, must be present an equal number of times (consecutive
% epochs) in c1 and c2, respectively
%
% ii_cfg.trialvec will contain the trial number associated with each sample
% ii_cfg.tcursel will contain beginning & end sample of each trial
%
% Example:
% load('exdata1.mat');
% [ii_data,ii_cfg] = ii_definetrial(ii_data,ii_cfg,'XDAT',1,'XDAT',8);
% figure; hold on;
% plot(ii_data.X); plot(ii_data.Y);
% myy = get(gca,'YLim');
% plot(repmat(ii_cfg.tcursel(:,1),1,2),myy,'k--');
% xlabel('Samples');
% title('Trial splits');

% Updated TCS 8/13/2017 to use ii_data, ii_cfg (no base vars)
% TODO: input cleaning, potentially add a sequential mode? (so that a
% channel is already the trial number?)

% =====================
% Adapted for only choose the epoch when participant focus on the central
% fixation, and 1st resp epoch,so the tcursel here would be 1 2 3 4 
% This function is adapted specific for plotting, so make sure it doesnt
% update the original ii_cfg

% Qingqing Yang, 6/12/2022

if nargin == 2
    prompt = {'Start when Channel', 'is at Value', 'Until Channel', 'is at Value','center_epoch','resp_epoch'};
    dlg_title = 'Trial Parameters';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines);
    c1 = answer{1};
    v1 = str2num(answer{2});
    c2 = answer{3};
    v2 = str2num(answer{4});
    center_epoch = str2num(answer{5});
    resp_epoch = str2num(answer{6});
    
end
ii_data2=ii_data;
ii_cfg2=ii_cfg;


% TODO: maybe relax case sensitivity?
if ismember(c1,fieldnames(ii_data2)) && ismember(c2,fieldnames(ii_data2))
    
    chan1 = ii_data2.(c1);
    chan2 = ii_data2.(c2);
    tsel     = chan1*0;
    trialvec = chan1*0;
    
    %% Get center_epoch and first resp_epoch
    swhere1 = find(diff(chan1 == center_epoch(1))== 1)+1;
    ewhere2 = find(diff(chan1 == resp_epoch)==-1);
    
%      in case where swhere(1) > ewhere(1), shift them relative to one
%      another (this can be true when c1 matches v1
     if (swhere1(1) > ewhere2(1)) && chan1(1)==v1
         swhere1 = [1; swhere1];
     end
     
     if length(ewhere2) == (length(swhere1)-1)
         ewhere2 = [ewhere2; length(chan2)];
     end
    
    tcursel = [swhere1 ewhere2];
    
    %% update ii_cfg2
    
    tindex = 1;
    
    ii_cfg2.tcursel = tcursel;
    ii_cfg2.trialvec = trialvec;
    
    ii_cfg2.history{end+1} = sprintf('ii_definetrial_respepoch %s %i to %s %i - %s',c1,v1,c2,v2,datestr(now,30));
    
    
end

end

