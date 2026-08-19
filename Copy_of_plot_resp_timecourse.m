function [trials, resp_grpd, align_stats, op_out] = plot_resp_timecourse(trials, op)
%% params
field_default('op','sort_cond','')
field_default('op','do_condition_sorting',1); 
field_default('op','plot_raster',0); 
field_default('op','trace_width',1);
field_default('op','cmapname','jet');
field_default('op','y_ax_hardlims',[]);
field_default('op','smooth_timecourses',1);
    field_default('op','smooth_method','gaussian');
    field_default('op','smooth_windowsize',30); 
field_default('op','condval_inds_to_plot',[]); 
field_default('op','newfig',1); 
field_default('op','yline_zero_width', 0.25); 
field_default('op','yline_zero_color',  [0.8 0.8 0.8]); 
field_default('op','yline_zero_style', '-');
% parameters for plotting event marker
field_default('op','xline_events',{}); 
field_default('op','include_xline_for_align_event',0); 
field_default('op','xline_event_line_width',2); 
field_default('op','xline_event_label_height',0.85); 
field_default('op','xline_event_label_border_width',1); 
field_default('op','leg_pos_adjust',0.22); 

%% plotting
if op.newfig
    hfig = figure('Color',[1 1 1]); box off
end

if ~op.do_condition_sorting
    align_stats = op.align_stats; 
elseif op.do_condition_sorting
    if strcmp(op.sort_cond,'') 
        trials.sort_cond = ones(height(trials),1); 
    elseif ~isempty(op.sort_cond)
        if iscell(op.sort_cond) 
            trials.sort_cond = trials{:,op.sort_cond{1}}(:,op.sort_cond{2}); 
        else
            trials.sort_cond = trials{:,op.sort_cond}; 
        end
    end
end

if ~strcmp(op.sort_cond,'')
    if op.do_condition_sorting
        [trials, align_stats, resp_grpd, op] = sort_responses_by_condition(trials,op); 
    elseif ~op.do_condition_sorting
        resp_grpd = op.resp_grpd; 
    end
    
    % plot error bars
    nconds = height(resp_grpd); 
    for icond = 1:nconds
        this_cond_sem_lims = [resp_grpd.resp_mean{icond} - resp_grpd.sem{icond}; resp_grpd.resp_mean{icond} + resp_grpd.sem{icond}]; 
        plotinds = resp_grpd.n_good_trials{icond} > 0; 
        if nnz(plotinds) > 0
            lowlims = this_cond_sem_lims(1,plotinds); 
            uplims = fliplr(this_cond_sem_lims(2,plotinds));
            if op.smooth_timecourses
                lowlims = smoothdata(lowlims, 2, op.smooth_method, op.smooth_windowsize); 
                uplims = smoothdata(uplims, 2, op.smooth_method, op.smooth_windowsize); 
            end
            hfill = fill([align_stats.times_aligned(plotinds), fliplr(align_stats.times_aligned(plotinds))], [lowlims,uplims], [0.8 0.8 0.8], 'HandleVisibility','off');
            hfill.LineStyle = 'none'; 
            hfill.EdgeColor = [0.8 0.8 0.8]; 
       end
       hold on
    end
    
    if isempty(op.condval_inds_to_plot)
        op.condval_inds_to_plot = 1:nconds;
    end
    nvals_to_plot = length(op.condval_inds_to_plot); 
    timecourses_to_plot = cell2mat(resp_grpd.resp_mean(op.condval_inds_to_plot,:))';
    if op.smooth_timecourses
        timecourses_to_plot = smoothdata(timecourses_to_plot, 1, op.smooth_method, op.smooth_windowsize); 
    end
    
    h_timecourse = plot(align_stats.times_aligned, timecourses_to_plot); 
    hax = gca;
    
    cmap = colormap;
    color_indices = round(linspace(1, size(cmap, 1), nvals_to_plot));
    set(h_timecourse, {'LineWidth'}, num2cell(repmat(op.trace_width, nvals_to_plot, 1)));
    set(h_timecourse, {'Color'}, num2cell(cmap(color_indices, :), 2));
    
    legend_strs = resp_grpd.condval; 

elseif strcmp(op.sort_cond,'')
     [trials, align_stats, op] = align_timecourses(trials, op);
    timecourses_to_plot = nanmean(trials.resp_aligned); 
    lowlims = align_stats.sem_lims(1,:); 
    uplims = align_stats.sem_lims(2,:); 
    if op.smooth_timecourses
        timecourses_to_plot = smoothdata(timecourses_to_plot, 2, op.smooth_method, op.smooth_windowsize); 
        lowlims = smoothdata(lowlims, 2, op.smooth_method, op.smooth_windowsize); 
        uplims = smoothdata(uplims, 2, op.smooth_method, op.smooth_windowsize); 
    end
    hold off
    hfill = fill([align_stats.xtime, fliplr(align_stats.xtime)], [lowlims, fliplr(uplims)], [0.8 0.8 0.8]);
    hfill.LineStyle = 'none'; 
    hfill.EdgeColor = [0.8 0.8 0.8]; 
    hold on  
    h_timecourse = plot(align_stats.xtime, timecourses_to_plot);
    h_timecourse.LineWidth = 1;
    resp_grpd = table; 
end

xlimits = op.samp_period * [-align_stats.n_tpoints_pre_fixed, align_stats.n_tpoints_post_fixed]; 
xlim(xlimits)
h_yline = yline(0,'LineWidth',op.yline_zero_width, 'Color',op.yline_zero_color ,'LineStyle',op.yline_zero_style);

%% plot xlines for specified timepoints
if ~isempty(op.xline_events)
    if size(op.xline_events, 1) == 1
        op.xline_events(2,:) = op.xline_events(1,:);
    end
    
    if ~op.include_xline_for_align_event
        op.xline_events(:, strcmp(op.xline_events(1,:), op.time_align_var)) = []; 
    end
    
    hax = gca;
    ax_pos = hax.Position;
    xlimits = xlim; 
    y_text_norm = ax_pos(2) + ax_pos(4) * op.xline_event_label_height;
    
    n_ev = size(op.xline_events, 2);
    all_x_values = NaN(n_ev, nconds); % Initialize with NaNs to cleanly ignore missing data
    all_colors = cell(n_ev, nconds); 
    mean_ev_times = NaN(n_ev, 1);
    
    sort_cond_data = trials.(op.sort_cond);
    align_var_data = trials{:, op.time_align_var};
    cond_colors = get(h_timecourse, {'Color'});
    if nvals_to_plot == 1, cond_colors = {cond_colors}; end 
    
    % FIX 1: Safely handle missing trials (NaNs) in the timing data
    safe_mean = @(x) mean(x, 'omitnan'); 
    
    for iev = 1:n_ev
        thislabel = op.xline_events{1,iev}; 
        t_diff = trials{:, thislabel} - align_var_data;
        
        [G, group_ids] = findgroups(sort_cond_data);
        mean_t = splitapply(safe_mean, t_diff, G);
        
        mean_ev_times(iev) = safe_mean(mean_t);
        
        for icond = 1:nconds
            thiscond = resp_grpd.condval{icond};
            if isnumeric(group_ids)
                idx = find(strcmp( cellstr(string(group_ids)), thiscond));
            else
                idx = find(strcmp(group_ids, thiscond));
            end
            
            if ~isempty(idx)
                all_x_values(iev, icond) = mean_t(idx);
            end
            all_colors{iev, icond} = cond_colors{icond};
        end
    end
    
    % Vectorized plotting of xlines (xline natively ignores NaNs)
    x_flat = all_x_values';
    colors_flat = all_colors';
    hxline = xline(x_flat(:), 'LineWidth', op.xline_event_line_width);
    set(hxline, {'Color'}, colors_flat(:));
    
    % Draw text boxes
    h_annno_event = cell(n_ev, 1); 
    for iev = 1:n_ev
        % FIX 2: Skip annotation drawing entirely if the mean time is NaN
        if isnan(mean_ev_times(iev))
            continue;
        end
        
        thislabel_to_plot = op.xline_events{2,iev}; 
        x_text_norm = ax_pos(1) + ax_pos(3) * (mean_ev_times(iev) - xlimits(1)) / (xlimits(2) - xlimits(1));
        
        h_annno_event{iev} = annotation('textbox', [x_text_norm, y_text_norm, 0.5, 0], ...
            'String', thislabel_to_plot, ...
            'Interpreter', 'none', ...
            'HorizontalAlignment', 'center', ... 
            'EdgeColor', 'k', ... 
            'LineWidth', op.xline_event_label_border_width, ... 
            'BackgroundColor', 'w', ... 
            'FitBoxToText', 'on');
    end
    
    drawnow; 
    for iev = 1:n_ev
        % FIX 3: Make sure we don't try to reposition a skipped (empty) annotation
        if isempty(h_annno_event{iev}) || ~isvalid(h_annno_event{iev})
            continue; 
        end
        
        x_text_norm = ax_pos(1) + ax_pos(3) * (mean_ev_times(iev) - xlimits(1)) / (xlimits(2) - xlimits(1));
        current_pos = h_annno_event{iev}.Position; 
        h_annno_event{iev}.Position(1) = x_text_norm - (current_pos(3) / 2); 
    end
end

xlabel('Time (sec)');
if isfield(op,'resp_signal')
    ylabel(['normed ', op.resp_signal, ' power']);
else
    ylabel('normed power');
end
set(gcf,'Color',[1 1 1]);
hleg = legend(legend_strs{:},'Interpreter','none');
hleg.Position(1) = hleg.Position(1)-op.leg_pos_adjust; 
hold off
box off

if op.plot_raster
    hfig_rater = figure('Color','w'); box off
    imagesc(resp_align.resp)
    ylabel('Trial')
end

op_out = op; 
trials_out = trials; 
end