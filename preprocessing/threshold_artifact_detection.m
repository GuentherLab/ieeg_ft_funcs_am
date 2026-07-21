function D = threshold_artifact_detection(D, cfg)
vardefault('cfg',struct);
field_default('cfg', 'threshold', 2)

y = [D.trial{:}];
y_mean = mean(y, 2, 'omitnan');
y_std = abs(std(y, 0, 2, 'omitnan'));

for i_trial = size(D.trial,1) % double check
    D.trial{i_trial}(D.trial{i_trial} > y_mean + cfg.threshold*y_std | D.trial{i_trial} < y_mean - cfg.threshold*y_std) = nan;
end

end