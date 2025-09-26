% Define the input signal parameters
Fs = 1000; % Sampling frequency (Hz)
T = 1/Fs;  % Sampling period (s)
L = 1000;  % Length of signal
t = (0:L-1)*T; % Time vector

% Create a real input signal (e.g., a sinusoid with noise)
f1 = 50; % Frequency of sinusoid 1
f2 = 120; % Frequency of sinusoid 2
X = 0.7*sin(2*pi*f1*t) + sin(2*pi*f2*t) + 0.5*randn(size(t));

% Perform the FFT
Y = fft(X);

% Compute the two-sided spectrum P2
P2 = abs(Y/L);

% Compute the single-sided spectrum P1
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); % Multiply by 2 for non-DC, non-Nyquist components

% Define the frequency vector
f = Fs*(0:(L/2))/L;

% Plot the single-sided amplitude spectrum
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('Frequency (Hz)')
ylabel('|P1(f)|')
grid on;