%% first method
SUBJECT = 'sub-DM1005';
TRIAL = 2;
ELECTRODE = 5;

PATH_DATASET = 'Y:\DBS';
PATH_DER = [PATH_DATASET filesep 'derivatives'];
PATH_SUB = [PATH_DER filesep SUBJECT];
PATH_FT = [PATH_SUB filesep 'fieldtrip'];
PATH_ANNOT = [PATH_SUB filesep 'annot'];

load([PATH_FT filesep SUBJECT '_ses-intraop_task-smsl_ft-raw_freqbands']);

%C = [0 2 4 6; 8 10 12 14; 16 18 20 22];
%image(C);

% rows are frequency, columns are time, data value is ecog reading
freqmap = trialed{1,1}.trial{1,TRIAL}(ELECTRODE,:);
for i=2:(length(trialed))
    freqmap = vertcat(trialed{1,i}.trial{1,TRIAL}(ELECTRODE,:),freqmap);
end

image(freqmap,'CDataMapping','scaled');
xlabel('time (s)');
ylabel('frequency');

%% second method
%{
SUBJECT = 'sub-DM1005';
TRIAL = 2;
ELECTRODE = 5;

PATH_DATASET = 'Y:\DBS';
PATH_DER = [PATH_DATASET filesep 'derivatives'];
PATH_SUB = [PATH_DER filesep SUBJECT];
PATH_FT = [PATH_SUB filesep 'fieldtrip'];
PATH_ANNOT = [PATH_SUB filesep 'annot'];

load([PATH_FT filesep SUBJECT '_ses-intraop_task-smsl_ft-raw-filt-trial-ar-ref']);

spectrogram(D_trial_ref.trial{1,TRIAL}(ELECTRODE,:),'yaxis');
%}