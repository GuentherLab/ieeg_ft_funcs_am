% first attempts at ASR of fluent speech

%% get audio from video
% Define file paths
inputFile = 'C:\ieeg_stut\sub-sap022\ses-3\Video\GH010063_task-reading_trimmed_530-605.mp4';
outputFile = 'C:\ieeg_stut\sub-sap022\ses-3\Video\GH010063_task-reading_trimmed_530-605_mono.wav';

% Extract audio
[y, Fs] = audioread(inputFile);

% Convert to mono by averaging all channels
y_mono = mean(y, 2);

% Export mono audio to WAV
audiowrite(outputFile, y_mono, Fs);

%% transcription
% Load Audio
wavFile = 'C:\ieeg_stut\sub-sap022\ses-3\Video\GH010063_task-reading_trimmed_530-605_mono.wav';
[y, Fs] = audioread(wavFile);

% Initialize Whisper Client (ModelSize options: "tiny", "base", "small", "medium", "large")
client = speechClient("whisper", ModelSize="medium"); 

% Transcribe and extract word timestamps table
transcriptTable = speech2text(y, Fs, Client=client);

% View word list, onset times, and confidence
disp(transcriptTable);
