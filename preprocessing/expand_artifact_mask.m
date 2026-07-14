function D = expand_artifact_mask(D, cfg)
% Inputs:
% D_in: Input fieldtrip with nans indicating where artifacts exist
% win_buffer_samples: double value indicating how many samples to turn into
% nans on each end of the artifact border
%
% Outputs:
% D_out: Output fieldtrip, identical to D_in but with expanded nan mask
% values

vardefault('cfg',struct);
field_default('cfg', 'fs', 100)
field_default('cfg', 'win_buffer_samples', 100);

% Identify artifact
artifact_mask_in = isnan(D.trial{1});
artifact_mask_out = artifact_mask_in;
[am_in_row, am_in_col] = find(artifact_mask_in);

% expand artifact mask 
for i_nan = 1:length(am_in_row)
    % identify instance of nan in original mask
    r = am_in_row(i_nan);
    c = am_in_col(i_nan);   

    % change +- win_buffer_samples in output    
    artifact_mask_out(r, max(1, c-cfg.win_buffer_samples):min(size(artifact_mask_out,2), c+cfg.win_buffer_samples)) = 1;
end

% apply expanded artifact mask
D.trial{1}(artifact_mask_out) = nan;
end