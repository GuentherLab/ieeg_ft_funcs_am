 % run early steps of intracranial trace cleaning: apply manual artifact masks, remove large amplitude spikes, high pass filter, notch filter
 %
 % outputs full-run trace, not trialwise
 %
 % bml repo must be in path
 % manual artifact mask should be created first - this will be loaded and applied in this function

function clean_mask_hpf_notch_filter(op)

vardefault('op',struct); % initialize options if not present
field_default('op','sub','DM1005')

op.art_crit = 'G'; % criteria G, implemented below, masks high amplitude jumps detected by the cleaning function

set_project_specific_variables(); % set paths etc. based on data collection site, load timing and electrode data

% filenamename that we will save filtered data into
ft_file_savepath = [FT_FILE_PREFIX 'raw-filt_ar-', op.art_crit]; 

load(FT_RAW_FILENAME)

% % % Adjusting length of sessions for notch filter to work
D_annot = bml_raw2annot(D);
Fl=[60 120 180 240 300 360 420 480];
D_annot.nSamples2 = round(floor(D_annot.nSamples .* Fl(1)./D_annot.Fs) .* D_annot.Fs./Fl(1));
D_annot.nSamples2 = round(floor(D_annot.nSamples2 .* Fl(1)./D_annot.Fs) .* D_annot.Fs./Fl(1));
D_annot.nSamples2 = D_annot.nSamples2(:,1);
% D_annot.nSamples2 = 938640; % hardcoded edit for DM1037
D_annot.ends = D_annot.starts + D_annot.nSamples2 ./ D_annot.Fs;

cfg=[];
cfg.epoch=D_annot;
D1= bml_redefinetrial(cfg,D);


% % % Selecting electrodes and remasking with zeros instead of NaNs
% % % Some of FieldTrip's functions don't work with NaNs, so we are going to temporarily replace NaNs with zeros to avoid issues. 

cfg=[];
cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'};
D_sel = ft_selectdata(cfg,D1);

cfg=[];
cfg.remask_nan = true;
cfg.value = 0;
D_sel = bml_mask(cfg, D_sel);

% ORIGINAL SIGNAL
spike_dur = 0.1; %100ms
og_sig = cat(2,D_sel.trial{:});                             % og_sig contains an electrodes*timepoints matrix of "original values"

% IDENTIFY POTENTIAL OUTLIERS (as samples where diff_sig > Q3 + n_thr * IQR, or diff_sig < Q1 - n_thr * IQR)
diff_sig = diff(og_sig,1,2);                                % diff_sig contains an electrodes*(timepoints-1) matrix of differences (temporal derivative)
m = 2*round(spike_dur*D_annot.Fs/2) + 1; % force m to be odd
diff_sig_smoothed = convn(diff_sig(:,max(1,min(size(diff_sig,2),1-(m-1)/2:size(diff_sig,2)+(m-1)/2))), hanning(m)', 'valid');

% [iqr_diff,qart_diff] = iqr(diff_sig_smoothed,2);   % iqr_diff: interquartile range of differences (IQR); qart_diff: first and third quartiles (Q1 & Q3) (REQUIRES >=R2024a)
iqr_diff = iqr(diff_sig_smoothed,2);
qart_diff = prctile(diff_sig_smoothed, [25; 75], 2);
iqr_thr = 3;                                                  % threshold to identify outliers


% RECONSTRUCTED SIGNAL
% Apply diff_sig_mask
% diff_sig = max(min(diff_sig, qart_diff(:,2)+n_thr*iqr_diff),qart_diff(:,1)-n_thr*iqr_diff); % crops derivatives beyond minimum/maximum values
diff_sig_mask = diff_sig_smoothed > qart_diff(:,2)+iqr_thr*iqr_diff | diff_sig_smoothed < qart_diff(:,1)-iqr_thr*iqr_diff;
diff_sig(diff_sig_mask) = 0;


% Apply Manual Artifact Mask
% load artifact mask
if exist('/Volumes/Nexus4/DBS/derivatives','dir') % if we're working in RD's local folder
    PATH_DER = '/Volumes/Nexus4/DBS/derivatives'; 
end
t = readtable([PATH_DER, filesep, 'sub-', op.sub, filesep 'annot',filesep, 'sub-',op.sub,'_ses-intraop_task-smsl_artifact-manual.tsv'], "FileType","text",'Delimiter', '\t');

% convert global time to samples
t.starts_idx = zeros(size(t.starts));
t.ends_idx = zeros(size(t.ends));
for i_t = 1:size(t,1)
    [~, t.starts_idx(i_t)] = min(abs(t.starts(i_t) - D_sel.time{1}));
    [~, t.ends_idx(i_t)] = min(abs(t.ends(i_t) - D_sel.time{1}));
end



%% AM proposal for changing this section: apply the manual mask, then use function hpf_and_instantaneous_artifact_mask rather than recreating that code here

% for each electrode, set value to zero from starts:end
for i_t = 1:size(t,1)
    diff_sig(strcmp(D_sel.label, t.label(i_t)), t.starts_idx(i_t):t.ends_idx(i_t)) = 0;
end

% apply
og_sig = cumsum([zeros(size(diff_sig,1),1), diff_sig],2);   % reconstructs original signal by cumulative sum (temporal integral)

% 
f_c = 2; % cutoff freq
k = 1/(1 + 2*pi*f_c/D_annot.Fs); % first order IIR
for n=2:size(og_sig,2)    
    og_sig(:,n) = k*og_sig(:,n-1) + k*diff_sig(:,n-1); % HPF
    % og_sig(:,n) = og_sig(:,n-1) + k*(diff_sig(:,n) - og_sig(:,n-1)); % LPF    
end
og_sig = og_sig(:,1:end-1); % AM edit 2026-6-26 - we need to remove a column to match original number of columns

D_sel.trial{:} = og_sig;
% clearvars diff_sig og_sig m_sig std_sig n_std
clearvars og_sig diff_sig diff_sig_mask diff_sig_smoothed iqr_diff qart_diff n_thr t




%% AM question - shouldn't we get rid of the highpass filter in this section? redundant with previous section... just turn off cfg.hpfilter options, leave dftfilter options
%% AM question - what's going on with the "if 0 elseif 1" section? 

% % % Applying high pass filter and line noise removal filter
% AM note: this step uses dft filter instead of bandstop filter (like artifact rejection uses).... ask Alan why the difference
cfg=[];
cfg.hpfilter=0; %%%%% this used to be true; AM changed to false
    cfg.hpfreq=1;
    cfg.hpfreq = 5;
cfg.hpfilttype='but';
cfg.hpfiltord=5;
cfg.hpfiltdir='twopass';
cfg.dftfilter='yes';
%%%% this interpolation option is causing an error with 3012, artifact criterion E
% cfg.dftreplace='neighbour';  %using spectrum interpolation method Mewett et al 2004
cfg.dftfreq           = [60 120 180 240 300 360 420 480];
cfg.dftbandwidth      = [1   1   1   1   1   1   1   1];
cfg.dftneighbourwidth = [2   2   2   2   2   2   2   2];

if 1 %% RD note: this should take <15mins
    F = cfg.dftfreq;
    N = D_annot.nSamples2;
    Fs = D_annot.Fs;
    BW = 2; % Bandwith of notch filter in Hz

    newF = F + reshape(find(abs(1/N*Fs*(0:N-1) - F(1))<BW/2)-1,[],1)/N*Fs - F(1);
    newF = reshape(newF,1,[]);
    % newF = [newF cfg.dftfreq(2:end)];
    cfg.dftfreq = newF;
elseif 0
    cfg.dftreplace='neighbour';
    cfg.dftbandwidth      = 0.05 * [ 1   1   1   1   1   1   1   1];
    cfg.dftneighbourwidth = 0.05 * [ 2   2   2   2   2   2   2   2];
end


%% AM question - first condition appears hardcoded - do we want to keep second condition? 

if 1 % && ~exist(['/Users/rohandeshpande/Documents/School/Research/Code/data/ft/sub-' op.sub '_ft_notch_' num2str(BW) 'Hz_cont.mat'],"file")
    D_sel_filt = ft_preprocessing(cfg,D_sel);
%     save(['/Users/rohandeshpande/Documents/School/Research/Code/data/ft/sub-' op.sub '_ft_notch_' num2str(BW) 'Hz_cont.mat'], 'D_sel_filt')
elseif exist(['/Users/rohandeshpande/Documents/School/Research/Code/data/ft/sub-' op.sub '_ft_notch_cont.mat'],"file")
%     load(['/Users/rohandeshpande/Documents/School/Research/Code/data/ft/sub-' op.sub '_ft_notch_' num2str(BW) 'Hz_cont.mat'])
end


%% save cleaned,HPF+notch-filtered traces on Y drive
save(ft_file_savepath,'D_sel_filt')









