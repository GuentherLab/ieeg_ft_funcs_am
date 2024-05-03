 %%%% plot the timecourse and SD of an electrode's response
 %%%% .... optionally sort trials by specific conditions 
 %%%% .... responses will be time-locked to a particular timepoint (e.g. speech onset)
 % 
 %%%% the main variable that must be present is a table 'timecourses_unaligned' with fields 'cond' (sort condition), 'timecourse', 'times', 'align_time'
 % this script is intended by be called by project-specific wrapper scripts, e.g. plot_resp_timecourse_triplet
 
%% params

vardefault('srt_row', 1);
vardefault('plot_raster',0); 

set_project_specific_variables()

channame = srt.chan{srt_row};
thissub = srt.sub{srt_row};
 
% remove trials with no response data
non_empty_trials = ~cellfun(@isempty, timecourses_unaligned.timecourse);
timecourses_unaligned = timecourses_unaligned(non_empty_trials,:);


%% align responses using alignment times provided within trials_tmp by the wrapper function
subind = string(subs.subject)==thissub; 

ntrials = height(timecourses_unaligned);
 nans_tr = nan(ntrials,1); 

 timecourses_unaligned = [timecourses_unaligned, table(nans_tr,              nans_tr,...
     'VariableNames',     {'tpoints_pre_onset', 'tpoints_post_onset'})];
 
 samp_period = 1e-5 * round(1e5 * diff(timecourses_unaligned.times{1}(1:2))); % sampling interval
 
 %%% find trial lengths pre- and post- the alignment time
for itrial = 1:ntrials
    % n timepoints before or at align_time
    timecourses_unaligned.tpoints_pre_onset(itrial) = nnz(timecourses_unaligned.times{itrial} <= timecourses_unaligned.align_time(itrial,1)); 
    % n timepoints after align_time
    timecourses_unaligned.tpoints_post_onset(itrial) = nnz(timecourses_unaligned.times{itrial} > timecourses_unaligned.align_time(itrial,1)); 
end
 
% pad or cut each trial to fit a specific size, so that we can align and average trials
switch trial_time_adj_method
    case 'median_plus_sd'
        n_tpoints_pre_fixed = round(median(timecourses_unaligned.tpoints_pre_onset) + std(timecourses_unaligned.tpoints_pre_onset)); 
        n_tpoints_post_fixed = round(median(timecourses_unaligned.tpoints_post_onset) + std(timecourses_unaligned.tpoints_post_onset)); 
    case 'median'
        n_tpoints_pre_fixed = median(timecourses_unaligned.tpoints_pre_onset); 
        n_tpoints_post_fixed = median(timecourses_unaligned.tpoints_post_onset); 
    case 'max'
        n_tpoints_pre_fixed = max(timecourses_unaligned.tpoints_pre_onset); 
        n_tpoints_post_fixed = max(timecourses_unaligned.tpoints_post_onset); 
end
tpoints_tot = n_tpoints_pre_fixed + n_tpoints_post_fixed; 

% nan-pad or cut trial windows so that they are all the same duration
%%% pad and cut values must be non-negative
resp_align = struct; 
resp_align.resp = NaN(ntrials, tpoints_tot); % aligned responses for this electrode; rows = trials, columns = timepoints
for itrial = 1:ntrials
   pre_pad = max([0, n_tpoints_pre_fixed - timecourses_unaligned.tpoints_pre_onset(itrial)]); 
   pre_cut = max([0, -n_tpoints_pre_fixed + timecourses_unaligned.tpoints_pre_onset(itrial)]); 
   pre_inds = 1+pre_cut:timecourses_unaligned.tpoints_pre_onset(itrial); % inds from timecourses_unaligned... if pre_cut > 0, some timepoints from this trial will not be used
   resp_align.resp(itrial, pre_pad+1 : n_tpoints_pre_fixed) = timecourses_unaligned.timecourse{itrial}(pre_inds); % fill in pre-onset data... fill in electrode responses starting after the padding epoch

   post_pad = max([0, n_tpoints_post_fixed - timecourses_unaligned.tpoints_post_onset(itrial)]);
   post_cut = max([0, -n_tpoints_post_fixed + timecourses_unaligned.tpoints_post_onset(itrial)]); 
   post_inds = timecourses_unaligned.tpoints_pre_onset(itrial) +  [1 : timecourses_unaligned.tpoints_post_onset(itrial)-post_cut]; % inds from timecourses_unaligned
   resp_align.resp(itrial, n_tpoints_pre_fixed+1:end-post_pad) = timecourses_unaligned.timecourse{itrial}(post_inds); % fill in post-onset data
   
   timecourses_unaligned.trial_onset_adjust(itrial) = samp_period * [pre_pad - pre_cut]; % number of timepoints to add to time landmarks
end
resp_align.mean = mean(resp_align.resp,'omitnan'); % mean response timecourse
resp_align.std = std(resp_align.resp, 'omitnan'); % stdev of response timecourses
resp_align.std_lims = [resp_align.mean + resp_align.std; resp_align.mean - resp_align.std]; 
resp_align.n_nonnan_trials = sum(~isnan(resp_align.resp)); % number of usable trials for this aligned timepoint
resp_align.sem = resp_align.std ./ sqrt(resp_align.n_nonnan_trials);
resp_align.sem_lims = [resp_align.mean + resp_align.sem; resp_align.mean - resp_align.sem]; 

xtime = 0.5 + [linspace(-n_tpoints_pre_fixed, -1, n_tpoints_pre_fixed), linspace(0, n_tpoints_post_fixed-1, n_tpoints_post_fixed)];
xtime = samp_period * xtime; 

%% plotting
if newfig
    hfig = figure('Color',[1 1 1]);
end

[unq_conds, ~, trial_cond_ind] = unique( timecourses_unaligned.cond );
if isnumeric(unq_conds)
    if  any(isnan(unq_conds)) % remove NaN condition labels
        [unq_conds, ~, trial_cond_ind] = unq_conds(~isnan(unq_conds));
    end
    timecourses_unaligned.cond = cellstr(num2str(timecourses_unaligned.cond));
    unq_conds = cellstr(num2str(unq_conds));
elseif iscell(unq_conds) % remove empty condition labels... probably occurs for unusable trials (e.g. no speech)
    unq_conds = unq_conds(~cellfun(@isempty,unq_conds)); 
end
nconds = length(unq_conds);

celcol = cell(nconds,1);
resp_grpd = table(unq_conds,celcol,celcol,'VariableNames',{'condval','resp','resp_mean'}); 

for icondval = 1:nconds
    these_trial_inds = trial_cond_ind == icondval;
    resp_grpd.resp{icondval} = resp_align.resp(these_trial_inds,:);
    resp_grpd.resp_mean{icondval} = mean(resp_grpd.resp{icondval},1,'omitnan');
end

% if grouping val indices not specified, plot them all
if isempty (condval_inds_to_plot)
    condval_inds_to_plot = 1:nconds;
end
nvals_to_plot = length(condval_inds_to_plot); 

% plot error bars
for icond = 1:nconds
    thiscond = unq_conds{icond};
    resp_rows_match = strcmp(timecourses_unaligned.cond, thiscond);
    this_cond_resp = resp_align.resp(resp_rows_match,:); % aligned trial timecourses for trials that match this condition label
    this_cond_mean = nanmean(this_cond_resp);
     this_cond_std = std(this_cond_resp, 'omitnan'); % stdev of response timecourses
    this_cond_n_nonnan_trials = sum(~isnan(this_cond_resp)); % number of usable trials for this aligned timepoint
    this_cond_sem = this_cond_std ./ sqrt(this_cond_n_nonnan_trials);
    this_cond_sem_lims = [this_cond_mean - this_cond_sem; this_cond_mean + this_cond_sem]; 
    plotinds = this_cond_n_nonnan_trials > 0; % timepoints with computable error bars

    if nnz(plotinds) > 0
        lowlims = this_cond_sem_lims(1,plotinds); 
        uplims = fliplr(this_cond_sem_lims(2,plotinds));
        if smooth_timecourses
            lowlims = smoothdata(lowlims, 2, smooth_method, smooth_windowsize); 
            uplims = smoothdata(uplims, 2, smooth_method, smooth_windowsize); 
        end

        hfill = fill([xtime(plotinds), fliplr(xtime(plotinds))], [lowlims,uplims], [0.8 0.8 0.8], 'HandleVisibility','off'); % standard error
        hfill.LineStyle = 'none'; % no border

        hfill.EdgeColor = [0.8 0.8 0.8]; 
   end
   hold on
end

timecourses_to_plot = cell2mat(resp_grpd.resp_mean(condval_inds_to_plot,:))';
if smooth_timecourses
    timecourses_to_plot = smoothdata(timecourses_to_plot, 1, smooth_method, smooth_windowsize); 
end
hplot = plot(xtime, timecourses_to_plot); 
%         hplot.LineWidth = 1;
hax = gca;
for ival = 1:nvals_to_plot
    hplot(ival).LineWidth = plotop.linewidth;
    colormap(cmapname)
    cmap = colormap;
    colormap_ind = round(size(cmap,1) * ival/nvals_to_plot);
    hplot(ival).Color = cmap(colormap_ind,:);
end

xlimdefault = xlim;
if isfield(plotop,'x_ax_hardlims') && ~isempty(plotop.x_ax_hardlims)
    xlim([max(plotop.x_ax_hardlims(1),xlimdefault(1)), min(plotop.x_ax_hardlims(2),xlimdefault(2))])
end

ylimdefault = ylim;
if isfield(plotop,'yx_ax_hardlims')  && ~isempty(plotop.y_ax_hardlims)
    ylim([max(plotop.y_ax_hardlims(1),ylimdefault(1)), min(plotop.y_ax_hardlims(2),ylimdefault(2))])
end

xlabel('Time (sec)')
ylabel([op.resp_signal, ' power (normed)'])

set(gcf,'Color',[1 1 1])

if ~isempty(sort_cond)
    hleg = legend(unq_conds,'Interpreter','none');
end

hold off



box off

%%
if plot_raster
    figure
    imagesc(resp_align.resp)
%     xlabel('Time (sec)')
    ylabel('Trial')


end