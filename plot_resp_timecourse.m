 %%%% plot the timecourse and SD of an electrode's response
 % load preprocessed trialwise data and run resp_types.m first
 %
 % align all trialwise responses to speech onset
 %
 
%% params

vardefault('srt_row', 1);

vardefault('sort_by_trial_cond',0)

% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

vardefault('plot_raster',0); 


set_project_specific_variables()

switch TASK
    case 'smsl'
        vardefault('sort_cond','learn_con');

    case 'triplet'

    otherwise 
        error('unknown task')
end

channame = srt.chan{srt_row};
thissub = srt.sub{srt_row};

resprow = strcmp(resp.chan,channame) & strcmp(resp.sub,thissub);
 



%% align responses
 % align responses to first syl onset

subind = string(subs.subject)==thissub; 

ntrials = height(trials_tmp);
 nans_tr = nan(ntrials,1); 

 trials_tmp = [trials_tmp, table(nans_tr,              nans_tr,...
     'VariableNames',     {'tpoints_pre_onset', 'tpoints_post_onset'})];
 
 samp_period = 1e-5 * round(1e5 * diff(trials_tmp.times{1}(1:2))); % sampling interval
 
 %%% find trial lengths pre- and post-onset
    %%%% trials.times{itrial} use global time coordinates
    %%%% ....... start at a fixed baseline window before stim onset
    %%%% ....... end at a fixed time buffer after speech offset
for itrial = 1:ntrials
    % n timepoints before or at speech onset [time from beginning of pre-stim baseline to speech onset]
    trials_tmp.tpoints_pre_onset(itrial) = nnz(trials_tmp.times{itrial} <= trials_tmp.t_prod_on(itrial,1)); 
    % n timepoints after speech onset [time from voice onset to end of post-speech time buffer]
    trials_tmp.tpoints_post_onset(itrial) = nnz(trials_tmp.times{itrial} > trials_tmp.t_prod_on(itrial,1)); 
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


%% plotting
xtime = 0.5 + [linspace(-n_tpoints_pre_fixed, -1, n_tpoints_pre_fixed), linspace(0, n_tpoints_post_fixed-1, n_tpoints_post_fixed)];
xtime = samp_period * xtime; 



% close all
    %     hfig = figure('Color',[1 1 1]);

if sort_by_trial_cond
    if isnumeric(unq_conds)
        [unq_conds, ~, trial_cond_ind] = unq_conds(~isnan(unq_conds));
        unq_conds = cellstr(num2str(unq_conds));
        trialconds  = strtrim(cellstr(num2str(trials_tmp{:,sort_cond})));
    else
        [unq_conds, ~, trial_cond_ind] = unique( trial_conds );
        trialconds = trials_tmp{:,sort_cond}; 
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
        condval_inds_to_plot = 1:ncondvals;
    end
    nvals_to_plot = length(condval_inds_to_plot); 

    % plot error bars
    %%% plotting error bars separately from means makes it simpler to create legend
    for icond = 1:nconds
        thiscond = unq_conds{icond};
        resp_rows_match = strcmp(trialconds, thiscond);
        this_cond_resp = resp_align.resp(resp_rows_match,:); % aligned trial timecourses for trials that match this condition label
        this_cond_mean = nanmean(this_cond_resp);
         this_cond_std = std(this_cond_resp, 'omitnan'); % stdev of response timecourses
        this_cond_n_nonnan_trials = sum(~isnan(this_cond_resp)); % number of usable trials for this aligned timepoint
        this_cond_sem = this_cond_std ./ sqrt(this_cond_n_nonnan_trials);
        this_cond_sem_lims = [this_cond_mean - this_cond_sem; this_cond_mean + this_cond_sem]; 
        plotinds = this_cond_n_nonnan_trials > 0; % timepoints with computable error bars

        if nnz(plotinds) > 0
            hfill = fill([xtime(plotinds), fliplr(xtime(plotinds))], [this_cond_sem_lims(1,plotinds), fliplr(this_cond_sem_lims(2,plotinds))], [0.8 0.8 0.8]); % standard error
            hfill.LineStyle = 'none'; % no border
    
            hfill.EdgeColor = [0.8 0.8 0.8]; 
       end
       hold on
    end


    % % % % % % % % plot means.... code originally for dbsseq version of function
    % % % % % % % for icond = 1:nconds
    % % % % % % %     thiscond = unq_conds{icond};
    % % % % % % %     resp_rows_match = strcmp(trialconds, thiscond);
    % % % % % % %     this_cond_resp = resp_align.resp(resp_rows_match,:); % aligned trial timecourses for trials that match this condition label
    % % % % % % %     this_cond_mean = nanmean(this_cond_resp);
    % % % % % % %     hplot(icond) = plot(xtime, this_cond_mean);
    % % % % % % % end

    timecourses_to_plot = cell2mat(resp_grpd.resp_mean(condval_inds_to_plot,:))';
    if smooth_timecourses
        timecourses_to_plot = smoothdata(timecourses_to_plot, 1, smooth_method, smooth_windowsize); 
    end
    hplot = plot(xtime, timecourses_to_plot); 
    %         hplot.LineWidth = 1;
    hax = gca;
    for ival = 1:nvals_to_plot
        hplot(ival).LineWidth = plotops.linewidth;
    end

    xlim(xlimits)

    ylimdefault = ylim;
    if ~isempty(y_ax_hardlims)
        ylim([max(y_ax_hardlims(1),ylimdefault(1)), min(y_ax_hardlims(2),ylimdefault(2))])
    end

    legend_strs = [repmat({''},nconds,1); unq_conds]; % empty entries match error bars

else


    hold off
    hfill = fill([xtime, fliplr(xtime)], [resp_align.sem_lims(1,:), fliplr(resp_align.sem_lims(2,:))], [0.8 0.8 0.8]); % standard error
%         hfill.LineStyle = 'none'; % no border
        hfill.EdgeColor = [0.8 0.8 0.8]; 
    hold on 
    hplot = plot(xtime, nanmean(resp_align.resp));
        hplot.LineWidth = 1;

    legend_strs = {''}; 
end

    xlabel('Time (sec)')
%     ylabel('HG power (normed)')
    ylabel('normed power')

    set(gcf,'Color',[1 1 1])

    hleg = legend(legend_strs{:},'Interpreter','none');
    hold off



box off

%%
if plot_raster
    figure
    imagesc(resp_align.resp)
%     xlabel('Time (sec)')
    ylabel('Trial')


end