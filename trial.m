clc
clear all
close all

[y, Fs] = audioread('10-18-2022-7-59-19__stereo199.wav');
N = length(y);
dt = 1/Fs;
T = (N-1)*dt;
t = 0:dt:T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time series observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all data
figure(1)
plot(t,y)
xlabel('time (s)')
ylabel('amplitude')
title('all data')

% a pop is heard at first 0.05s
timestamp_pop = 0.05;
idx_pop = timestamp_pop/T * N;
y_pop = y(1:idx_pop);
t_pop = t(1:idx_pop);
figure(2)
hold on
scatter(t_pop,y_pop,'.')
plot(t_pop,y_pop)
xlabel('time (s)')
ylabel('amplitude')
title(' y_{pop}(t) first 0.05s')
legend('scatter plot', 'line plot')
hold off

% transient at 0.5s can hear humming
timestap_hum1 = 0.5;
timestap_hum2 = 6;
idx_hum1 = timestap_hum1/T * N;
idx_hum2 = timestap_hum2/T * N;
y_hum = y(idx_hum1:idx_hum2);
t_hum = t(idx_hum1:idx_hum2);
figure(3)
plot(t_hum,y_hum)
xlabel('time (s)')
ylabel('amplitude')
title('y_{hum}(t) 0.5s - 6s')

audiowrite('y_hum.wav', 10*y_hum, Fs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% freq domain observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% during humming window
x = y_pop;
[P1, f] = fftplotter(x, Fs);

figure(4)
plot(f,20*log10(P1))
title("Single-Sided Power Spectrum of Y_{pop}(f)")
xlabel("f (Hz)")
ylabel("|Y_{hum}(f)| dB")

figure(5)
plot(f,P1)
title("Single-Sided Amplitude Spectrum of Y_{pop}(f)")
xlabel("f (Hz)")
ylabel("|Y_{hum}(f)|")

% during humming window
x = y_hum;
[P1, f] = fftplotter(x, Fs);

figure(6)
plot(f,20*log10(P1))
title("Single-Sided Power Spectrum of Y_{hum}(f)")
xlabel("f (Hz)")
ylabel("|Y_{hum}(f)| dB")

figure(7)
plot(f,P1)
title("Single-Sided Amplitude Spectrum of Y_{hum}(f)")
xlabel("f (Hz)")
ylabel("|Y_{hum}(f)|")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter humming window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low frequency bank
y_hum_low = lowpass(y_hum,100,Fs);
y_hum_low = 10*y_hum_low;

% time domain analysis
figure(8)
plot(t_hum,y_hum_low)
xlabel('time (s)')
ylabel('amplitude')
title('y_{hum_{low}}')

% freq domain analysis
x = y_hum_low;
[P1, f] = fftplotter(x, Fs);

figure(9)
plot(f,20*log10(P1))
title("Single-Sided Power Spectrum of Y_{hum_low}(f)")
xlabel("f (Hz)")
ylabel("|Y_{hum_{low}}(f)| dB")

figure(10)
plot(f,P1)
title("Single-Sided Amplitude Spectrum of Y_{hum_{low}}(f)")
xlabel("f (Hz)")
ylabel("|Y_{hum_{low}}(f)|")

audiowrite('y_hum_low.wav', y_hum_low, Fs)


% high frequency bank (38kHz - 43kHz)
y_hum_high = bandpass(y_hum,[38e3 43e3],Fs);

% downconvert
f_im = 38e3;
f_range = 5e3;
y_hum_high_converted = y_hum_high .* 10.*cos(2*pi*f_im.*t_hum);
y_hum_high_converted = lowpass(y_hum_high_converted,5e3,Fs);

% time domain analysis
figure(12)
plot(t_hum,y_hum_high_converted)
xlabel('time (s)')
ylabel('amplitude')
title('y_{hum_{high}}')


% freq domain analysis
x = y_hum_high_converted;
[P1, f] = fftplotter(x, Fs);

figure(13)
plot(f,20*log10(P1))
title("Single-Sided Power Spectrum of Y_{hum_{high converted}}(f)")
xlabel("f (Hz)")
ylabel("|Y_{hum_{high}}(f)| dB")

figure(14)
plot(f,P1)
title("Single-Sided Amplitude Spectrum of Y_{hum_{high converted}}(f)")
xlabel("f (Hz)")
ylabel("|Y_{hum_{high}}(f)|")

audiowrite('y_hum_high.wav', y_hum_high, Fs)
audiowrite('y_hum_high_converted.wav', y_hum_high_converted, Fs)

function [X, f] = fftplotter(x,Fs)
    L = length(x);
    Y = fft(x);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    X = P1;
end

