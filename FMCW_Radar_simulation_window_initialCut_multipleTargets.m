%% FMCW Radar Simulation with Multiple Targets + Range & Doppler FFT + Windows
clc; clear;

%% Radar Parameters
c = 3e8;                  % Speed of light
f0 = 12e9;                % Carrier frequency (12 GHz)
B = 30e6;                 % Sweep bandwidth
Tchirp = 20e-6;           % Chirp duration
Fs = 12.5e6;              % ADC sampling rate
Nsamples = 256;           % Range FFT size
Nchirps = 1024;           % Number of chirps
numCh = 4;                % Number of channels (I/Q pairs)

%% Derived parameters
lambda = c/f0;            % Wavelength
PRF = 1/Tchirp;           % Pulse repetition frequency

%% Initial cut settings
cutTime = 2e-6;           % 2 us initial cut
cutSamples = round(cutTime * Fs);

%% Multiple Targets Parameters
% Define as vectors
%%%Rtargets = [300, 600, 1200];        % meters
%%%Vtargets = [-50, 100, -200];        % m/s
Rtargets = [400, 700, 1000, 1100];       % meters
Vtargets = [-100, 50, -200, 100];        % m/s
numTargets = length(Rtargets);

%% Time vector for one chirp
t_full = (0:Nsamples-1)/Fs;
t = t_full(cutSamples+1:end);
Nsamples_cut = length(t);

%% Generate baseband signal for 4 channels
%txSig = zeros(Nsamples_cut, numCh);
%for ch = 1:numCh
%    sig_ch = zeros(Nsamples_cut,1);
%    for target = 1:numTargets
%        beatFreq = 2*B*Rtargets(target)/(c*Tchirp);  % Beat frequency for this target
%        fD = 2*Vtargets(target)*f0/c;               % Doppler frequency
%        phiDoppler = 2*pi*fD*(0:Nchirps-1)*Tchirp; % Doppler phase per chirp
%        sig_ch = sig_ch + exp(1j*2*pi*beatFreq*t'); % Sum targets in this channel
%    end
%    txSig(:,ch) = sig_ch;
%end

%% Stack signals for all chirps (simulate Doppler)
rxData = zeros(Nsamples_cut, Nchirps, numCh);
for ch = 1:numCh
    for k = 1:Nchirps
        % sum Doppler effect of multiple targets
        sig_ch = zeros(Nsamples_cut,1);
        for target = 1:numTargets
            beatFreq = 2*B*Rtargets(target)/(c*Tchirp);
            fD = 2*Vtargets(target)*f0/c;
            phiDoppler = 2*pi*fD*(k-1)*Tchirp;
            sig_ch = sig_ch + exp(1j*2*pi*beatFreq*t') * exp(1j*phiDoppler);
        end
        rxData(:,k,ch) = sig_ch;
    end
end

%% Range FFT window selection
rangeWindowType = 'blackman';  %%%'hanning'; % 'hanning' or 'blackman'
switch lower(rangeWindowType)
    case 'hanning'
        w_range = hann(Nsamples_cut);
    case 'blackman'
        w_range = blackman(Nsamples_cut);
end

%% Doppler FFT window selection
dopplerWindowType = 'blackmanharris'; % 'hanning' or 'blackmanharris'
switch lower(dopplerWindowType)
    case 'hanning'
        w_doppler = hann(Nchirps);
    case 'blackmanharris'
        w_doppler = blackmanharris(Nchirps);
end

%% Range FFT (window + zero-padding)
rangeFFT = zeros(Nsamples, Nchirps, numCh); % zero-pad to 256 points
for ch = 1:numCh
    for k = 1:Nchirps
        sig_win = rxData(:,k,ch) .* w_range;
        rangeFFT(:,k,ch) = fft(sig_win, Nsamples);
    end
end

%% Doppler FFT
dopplerFFT = zeros(Nsamples, Nchirps, numCh);
for ch = 1:numCh
    for r = 1:Nsamples
        %%%sig_dop = rangeFFT(r,:) .* w_doppler.';
        sig_dop = rangeFFT(r,:,ch) .* w_doppler.';
        dopplerFFT(r,:,ch) = fftshift(fft(sig_dop, Nchirps));
    end
end

%% Range & Velocity axes
rangeAxis = (0:Nsamples-1)*(Fs*Tchirp*c)/(2*B*Nsamples);        % meters
%%%dopplerAxis = linspace(-Fs/2, Fs/2, Nchirps)*lambda/2; % m/s
dopplerAxis = linspace(-PRF/2, PRF/2, Nchirps) * lambda/2;       % m/s

%% Plot Range-Doppler Map for channel 1
figure;
imagesc(dopplerAxis, rangeAxis, 20*log10(abs(dopplerFFT(:,:,1))));
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title(['Range-Doppler Map - Channel 1 - Multiple Targets - RangeWindow: ', rangeWindowType, ', DopplerWindow: ', dopplerWindowType]);
colorbar;
axis xy;
