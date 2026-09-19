 
%% this is the non-function version of this script - just kept for compatibility with AM's dbs-triplet scripts
% ---- next step when working on triplet should be updating triplet scripts to work with the updation function version, then archiving this non-function version

%%%% plot the timecourse and SD of an electrode's response
 %%%% .... optionally sort trials by specific conditions 
 %%%% .... responses will be time-locked to a particular timepoint (e.g. speech onset)
 % 
 % this script is intended by be called by project-specific wrapper scripts, e.g. plot_resp_timecourse_triplet and plot_resp_timecourse_seq
 
%% params

vardefault('srt_row', 1);
vardefault('plot_raster',0); 
vardefault('plotops',struct);
    field_default('plotops','linewidth',1);
vardefault('cmapname','jet');
vardefault('y_ax_hardlims',[]);
vardefault('times_to_plot',table()); 

field_default('op,','yline_zero_width', 0.25); 
field_default('op,','yline_zero_color',  [0.8 0.8 0.8]); 
field_default('op,','yline_zero_style', '-');

set_project_specific_variables()
 
% remove trials with no response data
non_empty_trials = ~cellfun(@isempty, timecourses_unaligned);
trials_tmp = trials_tmp(non_empty_trials,:);
    trial_conds = trial_conds(non_empty_trials,:); 
timecourses_unaligned = timecourses_unaligned(non_empty_trials);


%% align responses using alignment times provided within trials_tmp by the wrapper function
% subind = string(subs.subject)==thissub; 

ntrials = height(trials_tmp);
 nans_tr = nan(ntrials,1); 

 trials_tmp = [trials_tmp, table(nans_tr,              nans_tr,...
     'VariableNames',     {'tpoints_pre_onset', 'tpoints_post_onset'})];
 
 % compute sampling interval... assumes a static sample rate
 samp_period = 1e-5 * round(1e5 * diff(trials_tmp.times{1}(1:2))); 
 
 %%% find trial lengths pre- and post- the alignment time
for itrial = 1:ntrials
    % n timepoints before or at align_time
    trials_tmp.tpoints_pre_onset(itrial) = nnz(trials_tmp.times{itrial} <= trials_tmp.align_time(itrial,1)); 
    % n timepoints after align_time
    trials_tmp.tpoints_post_onset(itrial) = nnz(trials_tmp.times{itrial} > trials_tmp.align_time(itrial,1)); 
end
 
% pad or cut each trial to fit a specific size, so that we can align and average trials
switch trial_time_adj_method
    case 'median_plus_sd'
        n_tpoints_pre_fixed = round(median(trials_tmp.tpoints_pre_onset) + std(trials_tmp.tpoints_pre_onset)); 
        n_tpoints_post_fixed = round(median(trials_tmp.tpoints_post_onset) + std(trials_tmp.tpoints_post_onset)); 
    case 'median'
        n_tpoints_pre_fixed = median(trials_tmp.tpoints_pre_onset); 
        n_tpoints_post_fixed = median(trials_tmp.tpoints_post_onset); 
    case 'max'
        n_tpoints_pre_fixed = max(trials_tmp.tpoints_pre_onset); 
        n_tpoints_post_fixed = max(trials_tmp.tpoints_post_onset); 
end
tpoints_tot = n_tpoints_pre_fixed + n_tpoints_post_fixed; 

% nan-pad or cut trial windows so that they are all the same duration
%%% pad and cut values must be non-negative
resp_align = struct; 
resp_align.resp = NaN(ntrials, tpoints_tot); % aligned responses for this electrode; rows = trials, columns = timepoints
for itrial = 1:ntrials
   pre_pad = max([0, n_tpoints_pre_fixed - trials_tmp.tpoints_pre_onset(itrial)]); 
   pre_cut = max([0, -n_tpoints_pre_fixed + trials_tmp.tpoints_pre_onset(itrial)]); 
   pre_inds = 1+pre_cut:trials_tmp.tpoints_pre_onset(itrial); % inds from timecourses_unaligned... if pre_cut > 0, some timepoints from this trial will not be used
   resp_align.resp(itrial, pre_pad+1 : n_tpoints_pre_fixed) = timecourses_unaligned{itrial}(pre_inds); % fill in pre-onset data... fill in electrode responses starting after the padding epoch

   post_pad = max([0, n_tpoints_post_fixed - trials_tmp.tpoints_post_onset(itrial)]);
   post_cut = max([0, -n_tpoints_post_fixed + trials_tmp.tpoints_post_onset(itrial)]); 
   post_inds = trials_tmp.tpoints_pre_onset(itrial) +  [1 : trials_tmp.tpoints_post_onset(itrial)-post_cut]; % inds from timecourses_unaligned
   resp_align.resp(itrial, n_tpoints_pre_fixed+1:end-post_pad) = timecourses_unaligned{itrial}(post_inds); % fill in post-onset data
   
   trials_tmp.trial_onset_adjust(itrial) = samp_period * [pre_pad - pre_cut]; % number of timepoints to add to time landmarks
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

if ~isempty(sort_cond)
    [unq_conds, ~, trial_cond_ind] = unique( trial_conds );
    if isnumeric(unq_conds) % remove NaN condition labels
        [unq_conds, ~, trial_cond_ind] = unq_conds(~isnan(unq_conds));
        unq_conds = cellstr(num2str(unq_conds));
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
        resp_rows_match = strcmp(trial_conds, thiscond);
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
        hplot(ival).LineWidth = plotops.linewidth;
        cmap = colormap;
        colormap_ind = round(size(cmap,1) * ival/nvals_to_plot);
        hplot(ival).Color = cmap(colormap_ind,:);
    end

    

    ylimdefault = ylim;
    if ~isempty(y_ax_hardlims)
        ylim([max(y_ax_hardlims(1),ylimdefault(1)), min(y_ax_hardlims(2),ylimdefault(2))])
    end

    legend_strs = [repmat({''},nconds,1); unq_conds]; % empty entries match error bars

elseif isempty(sort_cond)

    timecourses_to_plot = nanmean(resp_align.resp); 
    lowlims = resp_align.sem_lims(1,:); % standard error
    uplims = resp_align.sem_lims(2,:); % standard error
    if smooth_timecourses
        timecourses_to_plot = smoothdata(timecourses_to_plot, 2, smooth_method, smooth_windowsize); 
        lowlims = smoothdata(lowlims, 2, smooth_method, smooth_windowsize); 
        uplims = smoothdata(uplims, 2, smooth_method, smooth_windowsize); 
    end
    hold off
    hfill = fill([xtime, fliplr(xtime)], [lowlims, fliplr(uplims)], [0.8 0.8 0.8]);
                hfill.LineStyle = 'none'; % no border
                hfill.EdgeColor = [0.8 0.8 0.8]; 
    hold on  
    hplot = plot(xtime, timecourses_to_plot);
        hplot.LineWidth = 1;

    legend_strs = {''}; 
end

%%% set xlims before the next section, which queries the xlims of the axis; if they aren't set first, matlab sometimes gets confused while querying them
vardefault('xlimits', samp_period * [-n_tpoints_pre_fixed, n_tpoints_post_fixed]); 
xlim(xlimits)

%%%% add y=0 line.... to remove, set op.yline_zero_width=0
h_yline = yline(0.3,'LineWidth',op.yline_zero_width, 'Color',op.yline_zero_color ,'LineStyle',op.yline_zero_style);

% plot xlines for all specified timepoints
%%% only implemented for seq, not triplet caller scripts
yproportion = 0.9; % how high up the plot to put the label text
for ilabel = 1:height(times_to_plot)
    thislabel = times_to_plot.varname{ilabel};
    t_thislabel = trials_tmp{:,thislabel}; % absolute times for this event on each trial
    t_post_aligntime_thislabel = t_thislabel - trials_tmp.align_time; % how long after the align time this event occurs each trial (negative if this event occurs pre-aligntime)
    mean_post_aligntime_thislabel = mean(t_post_aligntime_thislabel); % average window size between this event and aligntime across trials 
    h_xline(ilabel) = xline(mean_post_aligntime_thislabel, 'Color', times_to_plot.color{ilabel}); 

    % add text label for this timepoint
    switch times_to_plot.line_side{ilabel}; case 'R'; alignside = 'Left'; case 'L'; alignside = 'Right'; end % 'side' is the reverse of alignment
    hax = gca; 
    ax_pos = hax.Position;          % Convert axis coords to normalized figure coords....  [left bottom width height] in figure norm units
    xl = hax.XLim;
    yl = hax.YLim;
    xproportion = [mean_post_aligntime_thislabel - xl(1)] / diff(xl);
    xcoord = ax_pos(1) + xproportion*ax_pos(3);
    ycoord = ax_pos(2) + yproportion*ax_pos(4);
    annotation('textbox', [xcoord, ycoord, 0, 0], ...
        'String', times_to_plot.plot_label{ilabel}, ...
        'HorizontalAlignment', alignside, ...
        'VerticalAlignment', 'middle', ...
        'FitBoxToText', 'on', ...
        'EdgeColor', 'none');

end


xlabel('Time (sec)');
%     ylabel('HG power (normed)')
ylabel('normed power');

set(gcf,'Color',[1 1 1]);

% hleg = legend(legend_strs{:},'Interpreter','none');
hold off



box off

%%
if plot_raster
    figure
    imagesc(resp_align.resp)
%     xlabel('Time (sec)')
    ylabel('Trial')


end