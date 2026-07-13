 %%%% plot the timecourse and SD of an electrode's response
 %%%% .... optionally sort trials by specific conditions 
 %%%% .... responses will be time-locked to a particular timepoint (e.g. speech onset) specified by trials_tmp.align_time
%
%  op.time_align_var must be a variable in trials table
%  op.sort_cond can be any of:
%       '' = plot all trials together, unsorted
%       string = name of table variable with only one value per row
%       {string, index} = name of table variable with muliple values per row, index of column to use.... e.g. {consonant,2}
%
% outputs: 
%       1. trials = original trials table appended with resp_aligned (responses aligned to intratrial event of inerest)
%       2. resp_grpd = table with row for each trial condition, and the responses occurring in each trial of that condition
%       3. align_stats = struct with fields containing simple analyses of aligned timecourses, including timecourse mean, sem, sem bar lims (for plotting), timepoints on each side of sync point
%               .... this contains align_stats.xtime added - match this with trials.resp_aligned for plotting
%       4. op = original op struct plus defaults that were filled in
 %
 % this script is intended by be called by project-specific wrapper scripts, e.g. plot_resp_timecourse_triplet and plot_resp_timecourse_seq
 
 function [trials,resp_grpd, align_stats, op_out] = plot_resp_timecourse(trials, op)

%% params

field_default('op','sort_cond','')
field_default('op','plot_raster',0); 
field_default('op','trace_width',1);
field_default('op','cmapname','jet');
field_default('op','y_ax_hardlims',[]);
field_default('op','times_to_plot',table()); 

field_default('op','yline_zero_width', 0.25); 
field_default('op','yline_zero_color',  [0.8 0.8 0.8]); 
field_default('op','yline_zero_style', '-');
field_default('op','y_timelabel_height',0.8); % how high up the plot to put the timepoint label text


 

if strcmp(op.sort_cond,'') % plot all trials in a single trace
    trials.sort_cond = ones(height(trials),1); 
elseif ~isempty(op.sort_cond)
    if iscell(op.sort_cond) % if we need to select only 1 column from the table variable
        trials.sort_cond = trials{:,op.sort_cond{1}}(:,op.sort_cond{2}); 
    else
        trials.sort_cond = trials{:,op.sort_cond}; 
    end



end











%% plotting
if op.newfig
    hfig = figure('Color',[1 1 1]); box off
end


if ~strcmp(op.sort_cond,'')


     [trials, align_stats, resp_grpd, op] = sort_responses_by_condition(trials,op); 

    % plot error bars
    nconds = height(resp_grpd); 
    for icond = 1:nconds
        this_cond_sem_lims = [resp_grpd.resp_mean{icond} - resp_grpd.sem{icond}; resp_grpd.resp_mean{icond} + resp_grpd.sem{icond}]; 
        plotinds = resp_grpd.n_good_trials{icond} > 0; % timepoints with computable error bars

        if nnz(plotinds) > 0
            lowlims = this_cond_sem_lims(1,plotinds); 
            uplims = fliplr(this_cond_sem_lims(2,plotinds));
            if op.smooth_timecourses
                lowlims = smoothdata(lowlims, 2, op.smooth_method, op.smooth_windowsize); 
                uplims = smoothdata(uplims, 2, op.smooth_method, op.smooth_windowsize); 
            end

            hfill = fill([align_stats.xtime(plotinds), fliplr(align_stats.xtime(plotinds))], [lowlims,uplims], [0.8 0.8 0.8], 'HandleVisibility','off'); % standard error
                hfill.LineStyle = 'none'; % no border
                hfill.EdgeColor = [0.8 0.8 0.8]; 
       end
       hold on
    end

    % if grouping val indices not specified, plot them all
    if isempty (op.condval_inds_to_plot)
        op.condval_inds_to_plot = 1:nconds;
    end
    nvals_to_plot = length(op.condval_inds_to_plot); 

    timecourses_to_plot = cell2mat(resp_grpd.resp_mean(op.condval_inds_to_plot,:))';
    if op.smooth_timecourses
        timecourses_to_plot = smoothdata(timecourses_to_plot, 1, op.smooth_method, op.smooth_windowsize); 
    end
    hplot = plot(align_stats.xtime, timecourses_to_plot); 
    %         hplot.LineWidth = 1;
    hax = gca;
    for ival = 1:nvals_to_plot
        hplot(ival).LineWidth = op.trace_width;
        cmap = colormap;
        colormap_ind = round(size(cmap,1) * ival/nvals_to_plot);
        hplot(ival).Color = cmap(colormap_ind,:);
    end

   

    legend_strs = [repmat({''},nconds,1); resp_grpd.condval]; % empty entries match error bars

elseif strcmp(op.sort_cond,'')

     [trials, align_stats, op] = align_timecourses(trials, op);

    timecourses_to_plot = nanmean(trials.resp_aligned); 
    lowlims = align_stats.sem_lims(1,:); % standard error
    uplims = align_stats.sem_lims(2,:); % standard error
    if op.smooth_timecourses
        timecourses_to_plot = smoothdata(timecourses_to_plot, 2, op.smooth_method, op.smooth_windowsize); 
        lowlims = smoothdata(lowlims, 2, op.smooth_method, op.smooth_windowsize); 
        uplims = smoothdata(uplims, 2, op.smooth_method, op.smooth_windowsize); 
    end
    hold off
    hfill = fill([align_stats.xtime, fliplr(align_stats.xtime)], [lowlims, fliplr(uplims)], [0.8 0.8 0.8]);
                hfill.LineStyle = 'none'; % no border
                hfill.EdgeColor = [0.8 0.8 0.8]; 
    hold on  
    hplot = plot(align_stats.xtime, timecourses_to_plot);
        hplot.LineWidth = 1;

    legend_strs = {''}; 

    resp_grpd = table; 
end

%%% set xlims before the next section, which queries the xlims of the axis; if they aren't set first, matlab sometimes gets confused while querying them
xlimits = op.samp_period * [-align_stats.n_tpoints_pre_fixed, align_stats.n_tpoints_post_fixed]; 
xlim(xlimits)

%%%% add y=0 line.... to remove, set op.yline_zero_width=0
h_yline = yline(0.3,'LineWidth',op.yline_zero_width, 'Color',op.yline_zero_color ,'LineStyle',op.yline_zero_style);

% plot xlines for all specified timepoints
%%% only implemented for seq, not triplet caller scripts
for ilabel = 1:height(op.times_to_plot)
    thislabel = op.times_to_plot.varname{ilabel};
    t_thislabel = trials{:,thislabel}; % absolute times for this event on each trial
    t_post_aligntime_thislabel = t_thislabel - trials.align_time; % how long after the align time this event occurs each trial (negative if this event occurs pre-aligntime)
    
    % need to use nanmean here because of e.g. go trials with no speech production, which would have t_prod_on==nan
    mean_post_aligntime_thislabel = nanmean(t_post_aligntime_thislabel); % average window size between this event and aligntime across trials 
    h_xline(ilabel) = xline(mean_post_aligntime_thislabel, 'Color', op.times_to_plot.color{ilabel}); 

    % add text label for this timepoint
    switch op.times_to_plot.line_side{ilabel}; case 'R'; alignside = 'Left'; case 'L'; alignside = 'Right'; end % 'side' is the reverse of alignment
    hax = gca; 
    ax_pos = hax.Position;          % Convert axis coords to normalized figure coords....  [left bottom width height] in figure norm units
    xl = hax.XLim;
    yl = hax.YLim;
    xproportion = [mean_post_aligntime_thislabel - xl(1)] / diff(xl);
    xcoord = ax_pos(1) + xproportion*ax_pos(3);
    ycoord = ax_pos(2) + op.y_timelabel_height*ax_pos(4);
    annotation('textbox', [xcoord, ycoord, 0, 0], ...   %%% need to debug for nans
        'String', op.times_to_plot.plot_label{ilabel}, ...
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
if op.plot_raster
    hfig_rater = figure('Color','w'); box off
    imagesc(resp_align.resp)
%     xlabel('Time (sec)')
    ylabel('Trial')


end

op_out = op; 
trials_out = trials; 

end