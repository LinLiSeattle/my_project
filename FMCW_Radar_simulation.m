%% FMCW Radar Simulation with Range & Doppler FFT
clc;
close all;
clear;

%% Radar Parameters
c = 3e8;                  % Speed of light
f0 = 12e9;                % Carrier frequency (12 GHz)
B = 30e6;                 % Sweep bandwidth (30 MHz)
Tchirp = 20e-6;           % Chirp duration
Fs = 12.5e6;              % ADC sampling rate
Nsamples = 256;           % Range FFT size (samples per chirp)
Nchirps = 1024;           % Number of chirps (for Doppler FFT)
numCh = 4;                % Number of channels (I/Q pairs: I1/Q1 ...)

%% Target parameters (for simulation)
Rtarget = 800;   %%%1.5e3;          % Target range (meters)
Vtarget = -30;     %%%300;            % Target radial velocity (m/s)
fD = 2*Vtarget*f0/c;      % Doppler frequency

%% Time vector for one chirp
t = (0:Nsamples-1)/Fs;

%% Generate FMCW chirp baseband signal for 4 channels
txSig = zeros(Nsamples, numCh);
for ch = 1:numCh
    % Simple single target per channel, I/Q simulation
    beatFreq = 2*B*Rtarget/(c*Tchirp);      % Beat frequency
    phiDoppler = 2*pi*fD*(0:Nchirps-1)*Tchirp; % Doppler phase per chirp
    % replicate over chirps
    txSig(:,ch) = exp(1j*2*pi*beatFreq*t'); % baseband signal for one chirp
end

%% Stack signals for all chirps (simulate Doppler)
rxData = zeros(Nsamples, Nchirps, numCh);
for ch = 1:numCh
    for k = 1:Nchirps
        rxData(:,k,ch) = txSig(:,ch) .* exp(1j*phiDoppler(k));
    end
end

%% Range FFT
rangeFFT = fft(rxData, Nsamples, 1);  % FFT along fast-time (range)

%% Doppler FFT
dopplerFFT = fftshift(fft(rangeFFT, Nchirps, 2),2);  % FFT along slow-time (Doppler)

%% Calculate range & velocity axes
lambda = c/f0;  %%% moved here
rangeAxis = (0:Nsamples-1)*(Fs*Tchirp*c)/(2*B*Nsamples); % meters
%%%dopplerAxis = linspace(-Fs/2, Fs/2, Nchirps)*lambda/2;   % m/s approximation
fd_max = 1 / (2 * Tchirp);
dopplerAxis = linspace(-fd_max, fd_max, Nchirps) * lambda / 2;

%%%lambda = c/f0;

%% Plot Range-Doppler Map for first channel
figure;
imagesc(dopplerAxis, rangeAxis, 20*log10(abs(dopplerFFT(:,:,1))));
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Range-Doppler Map (Channel 1)');
colorbar;
axis xy;

