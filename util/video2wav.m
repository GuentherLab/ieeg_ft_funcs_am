function wavFilePath = video2wav(videoFilePath)
% video2wav - Extracts audio from a video file and saves it as a WAV file.
%
%   wavFilePath = extractAudioToWav() opens a file selection dialog.
%   wavFilePath = extractAudioToWav(videoFilePath) processes the specified file.
%
%   Supports MP4, MOV, MKV, WebM, AVI, and WMV. Preserves original sample rate
%   and channel count, saving as uncompressed 16-bit PCM WAV (or 24-bit if
%   source bit depth exceeds 16-bit).

    % 1. Prompt user if no file path is provided
    if nargin < 1 || isempty(videoFilePath)
        filterSpec = {'*.mp4;*.mov;*.mkv;*.webm;*.avi;*.wmv', ...
                      'Video Files (*.mp4, *.mov, *.mkv, *.webm, *.avi, *.wmv)'; ...
                      '*.*', 'All Files (*.*)'};
        [fileName, pathName] = uigetfile(filterSpec, 'Select a Video File');
        
        if isequal(fileName, 0)
            disp('File selection canceled.');
            wavFilePath = '';
            return;
        end
        videoFilePath = fullfile(pathName, fileName);
    end

    % Validate file existence
    if ~isfile(videoFilePath)
        error('extractAudioToWav:FileNotFound', 'File not found: %s', videoFilePath);
    end

    % 2. Build output WAV path in the same directory
    [folder, baseName, ~] = fileparts(videoFilePath);
    wavFilePath = fullfile(folder, [baseName, '.wav']);

    % 3. Check if FFmpeg is available on the system path
    [status, ~] = system('ffmpeg -version');
    hasFFmpeg = (status == 0);

    if hasFFmpeg
        % Probe original bit depth to preserve high-res audio (e.g., 24-bit) if present
        probeCmd = sprintf('ffprobe -v error -select_streams a:0 -show_entries stream=bits_per_raw_sample,bits_per_sample -of default=noprint_wrappers=1:nokey=1 "%s"', videoFilePath);
        [pStatus, pOut] = system(probeCmd);
        
        codec = 'pcm_s16le'; % Standard lossless CD-quality WAV (16-bit)
        if pStatus == 0
            bits = str2double(regexp(pOut, '\d+', 'match', 'once'));
            if ~isnan(bits) && bits > 16
                codec = 'pcm_s24le'; % Preserve 24-bit / 32-bit studio audio
            end
        end

        % Extract audio: -vn (no video), -acodec (uncompressed PCM), -y (overwrite)
        % Sample rate and channel count are preserved automatically by default.
        cmd = sprintf('ffmpeg -y -i "%s" -vn -acodec %s "%s"', ...
                      videoFilePath, codec, wavFilePath);
        [cmdStatus, cmdOut] = system(cmd);

        if cmdStatus ~= 0
            error('extractAudioToWav:FFmpegError', ...
                  'FFmpeg failed to extract audio:\n%s', cmdOut);
        end
    else
        % Fallback to native MATLAB audioread/audiowrite (works for MP4/MOV/WMV on most OS)
        warning('extractAudioToWav:NoFFmpeg', ...
                'FFmpeg not found on system PATH. Falling back to native audioread (may not support MKV/WebM).');
        try
            info = audioinfo(videoFilePath);
            [y, Fs] = audioread(videoFilePath);
            
            % Match source bit depth if available (8, 16, 24, 32)
            bitDepth = 16;
            if isfield(info, 'BitsPerSample') && ~isempty(info.BitsPerSample)
                bitDepth = info.BitsPerSample;
            end
            
            audiowrite(wavFilePath, y, Fs, 'BitsPerSample', bitDepth);
        catch ME
            error('extractAudioToWav:NativeReadFailed', ...
                  'Native extraction failed (install FFmpeg for full format support).\nDetails: %s', ME.message);
        end
    end

    fprintf('Saved WAV audio to: %s\n', wavFilePath);
end
