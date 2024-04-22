%%% get average normalized power across multiple freq bands for an electrode type
%
%%% inputs: 
    % cfg.param = row of 'parameter' table used for BML artifact rejeection
    % cfg.out_freq = sampling freq of output data (argument to bml_envelope_wavpow.m)
    % cfg.suppress_output = whether to print bml_envelope_wavpow.m output to command line
    %
    % D_in = fieldtrip struct with electrode data to be processed.... usually highpass filtered first

function D_avg_pow_eltype = multifreq_avg_power(cfg, D_in)

field_default('cfg','out_freq',100)
field_default('cfg','suppress_output',1)

param = cfg.param; 
ntrials = numel(D_in.trial);

el_type = strip(param.electrode_type{1});
wav_width = param.wav_width(1);

%selecting channels for artifact rejection
cfg1=[];
cfg1.channel = [el_type,'_*'];
D_eltype = ft_selectdata(cfg1,D_in);

if isempty(D_eltype.label)
    %channel type not available
    D_avg_pow_eltype = []; 
elseif ~isempty(D_eltype.label)


    % compute log-spaced frequencies between wav_freq_min and wav_freq_max
    nfreqs = param.n_wav_freqs(1); 
    wav_freqs = round(logspace(log10(param.wav_freq_min(1)),log10(param.wav_freq_max(1)),nfreqs));
    D_multifreq_eltype = cell(nfreqs,1);
    
    normed_pow = cell(1,ntrials); 
    for ifreq = 1:nfreqs
      %calculating absolute value envelope at 1Hz (1s chunks)
      cfg1=[];
      cfg1.out_freq = cfg.out_freq;
      cfg1.wav_freq = wav_freqs(ifreq);
      cfg1.wav_width = wav_width;
      if cfg.suppress_output
        cmd  = 'D_multifreq_eltype{ifreq} = bml_envelope_wavpow(cfg1,D_eltype);';
        evalc(cmd); % use evalc to suppress console output
      elseif ~cfg.suppress_output
        D_multifreq_eltype{ifreq} = bml_envelope_wavpow(cfg1,D_eltype);
      end
      
      nchannels = length(D_multifreq_eltype{ifreq}.label);
      D_multifreq_eltype{ifreq}.med_pow_per_block = NaN(nchannels, ntrials); % initialize
      for iblock = 1:ntrials % for each block, normalize by median power
        % rows are channels, so take the median across columns (power at timepoints for each channel)
          D_multifreq_eltype{ifreq}.med_pow_per_block(:,iblock) = median(D_multifreq_eltype{ifreq}.trial{iblock},2);
          % normalize power by median values within each channel for this block
          %%% normed_pow will be filled with all normed powers across blocks and frequencies; we will average across the 3rd dimension (frequency)
          normed_pow{iblock}(:,:,ifreq) = D_multifreq_eltype{ifreq}.trial{iblock} ./ D_multifreq_eltype{ifreq}.med_pow_per_block(:,iblock);
      end
    end
    
    D_avg_pow_eltype = struct; % averaged high gamma
        D_avg_pow_eltype.hdr = D_multifreq_eltype{1}.hdr;
        D_avg_pow_eltype.trial = D_multifreq_eltype{1}.trial;
        D_avg_pow_eltype.trial = cell(1,ntrials); % to be filled
        D_avg_pow_eltype.time = D_multifreq_eltype{1}.time;
        D_avg_pow_eltype.label = D_multifreq_eltype{1}.label;
        if isfield(D_multifreq_eltype{1},'sampleinfo')
            D_avg_pow_eltype.sampleinfo = D_multifreq_eltype{1}.sampleinfo;
        end
        
    %%%%% get averaged wave power
    % put median power across into a single array, with dimorder chans-trials-freq 
    med_pow = cell2mat(permute(cellfun(@(x)x.med_pow_per_block,D_multifreq_eltype,'UniformOutput',false), [2 3 1])); 
    med_pow_mean = mean(med_pow,3); % channel/block median powers, averaged across frequencies of interest
    for iblock = 1:ntrials
        D_avg_pow_eltype.trial{iblock} = mean(normed_pow{iblock},3);
        % multiply the [channel/block]-specific median powers back, to differentiate between absolute power values of channels and trials
        D_avg_pow_eltype.trial{iblock} = D_avg_pow_eltype.trial{iblock}  .* med_pow_mean(:,iblock); 
    end

end
