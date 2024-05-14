%% project2_400101186_400101045_400101164

%% increase the sound step

clear; clc ;
% Load the voice file
[input_wave, Fs] = audioread('test5.mp3');
% Load the input voice file into the variable input_wave and get the sampling frequency Fs

% Define the pitch factor (e.g., to increase pitch by a factor of 1.5, set pitch_factor = 1.5)
pitch_factor = 1.5;
% Define the pitch factor by which the pitch will be changed (increased or decreased)

% Calculate the new length of the signal based on the pitch factor
new_length = round(numel(input_wave) / pitch_factor);
% Calculate the new length of the output signal based on the pitch factor and the length of the input signal

% Initialize the output signal
output_wave = zeros(new_length, 1);
% Create an empty array of zeros to store the modified output signal

% Implement the PSOLA method
for i = 1:new_length
    % Calculate the original index using the pitch factor
    original_index = min(round(i * pitch_factor), numel(input_wave));
    % Calculate the corresponding index in the original signal based on the pitch factor and the current index of the output signal
    
    % Copy the sample from the original signal to the modified signal
    output_wave(i) = input_wave(original_index);
    % Copy the sample from the original signal to the modified output signal at the current index
end

% Plot the FFT of the main voice
subplot(2,2,1);
Y = fft(input_wave);
L = length(input_wave);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1);
title('Single-Sided Amplitude Spectrum of Main Voice');
xlabel('f (Hz)');
ylabel(" main voice");

% Plot the FFT of the output voice
subplot(2,2,3);
Yout = fft(output_wave);
Lout = length(output_wave);
P2out = abs(Yout/Lout);
P1out = P2out(1:Lout/2+1);
P1out(2:end-1) = 2*P1out(2:end-1);
fout = Fs*(0:(Lout/2))/Lout;
plot(fout,P1out);
title('Single-Sided Amplitude Spectrum of Output Voice');
xlabel('f (Hz)');
ylabel("output voice");

% Plot the main voice
subplot(2,2,2);
plot(input_wave);
title('Main Voice');
xlabel('Sample');
ylabel('Amplitude');

% Plot the secondary (output) voice
subplot(2,2,4);
plot(output_wave);
title('Secondary (Output) Voice');
xlabel('Sample');
ylabel('Amplitude');
% Play the output signal
sound(output_wave, Fs);
% Play the modified output signal at the original sampling frequency

% Save the output signal to a file
audiowrite('test5_output_increase.wav', output_wave, Fs);
% Play the main voice
% sound(input_wave, Fs);
% Play the original input voice
% which('test1_output_increase.wav');

%%  decrease the sound step

clear; clc ;
% Load the voice file
[input_wave, Fs] = audioread('test5.mp3');
% Load the input voice file into the variable input_wave and get the sampling frequency Fs

% Define the pitch factor (e.g., to decrease pitch by a factor of 1.5, set pitch_factor = 0.5)
pitch_factor = 0.5;
% Define the pitch factor by which the pitch will be changed (increased or decreased)

% Calculate the new length of the signal based on the pitch factor
new_length = round(numel(input_wave) / pitch_factor);
% Calculate the new length of the output signal based on the pitch factor and the length of the input signal

% Initialize the output signal
output_wave = zeros(new_length, 1);
% Create an empty array of zeros to store the modified output signal

% Implement the PSOLA method
for i = 1:new_length
    % Calculate the original index using the pitch factor
    original_index = min(round(i * pitch_factor), numel(input_wave));
    % Calculate the corresponding index in the original signal based on the pitch factor and the current index of the output signal
    
    % Copy the sample from the original signal to the modified signal
    output_wave(i) = input_wave(original_index);
    % Copy the sample from the original signal to the modified output signal at the current index
end

% Plot the FFT of the main voice
subplot(2,2,1);
Y = fft(input_wave);
L = length(input_wave);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1);
title('Single-Sided Amplitude Spectrum of Main Voice');
xlabel('f (Hz)');
ylabel(" main voice");

% Plot the FFT of the output voice
subplot(2,2,3);
Yout = fft(output_wave);
Lout = length(output_wave);
P2out = abs(Yout/Lout);
P1out = P2out(1:Lout/2+1);
P1out(2:end-1) = 2*P1out(2:end-1);
fout = Fs*(0:(Lout/2))/Lout;
plot(fout,P1out);
title('Single-Sided Amplitude Spectrum of Output Voice');
xlabel('f (Hz)');
ylabel("output voice");

% Plot the main voice
subplot(2,2,2);
plot(input_wave);
title('Main Voice');
xlabel('Sample');
ylabel('Amplitude');

% Plot the secondary (output) voice
subplot(2,2,4);
plot(output_wave);
title('Secondary (Output) Voice');
xlabel('Sample');
ylabel('Amplitude');
% Play the output signal
sound(output_wave, Fs);
% Play the modified output signal at the original sampling frequency

% Save the output signal to a file
audiowrite('test5_output_decrease.wav', output_wave, Fs);

% Play the main voice
% sound(input_wave, Fs);
% Play the original input voice

%% change spectrum of the input voice by delta_shifting

clear; clc;

load music.mat
[Y,FS] = audioread("test3.mp3");
%sound(Y,FS);

leftchannel= Y(: , 1);
rightchannel= Y(: , 2);
y1=leftchannel+rightchannel;
sound1=(1/2).*y1;
sound1 = transpose(sound1);

%sound(sound1,FS);
t_Max = length(sound1)/FS;
t = linspace(0,t_Max,length(sound1));

figure
plot(t,sound1,LineWidth=2);
xlabel("t");
ylabel("input voice signal")

% spectrum of the input voice

X = fft(sound1);
X2 = fftshift(X);
l = length(sound1);
fshift = (-l/2:l/2-1)*(FS/l);

figure
plot(fshift,X2,LineWidth=2);
xlabel("frequency range (Hz)");
ylabel("F(f)");
title("fourier transform of input voice signal");

% spectrum of the result of multipication by cos(2*pi*f1*t)

Cos1 = cos(2*pi*9000.*t);
Multipication1 = sound1.*Cos1;

figure
plot(t,Multipication1,LineWidth=2);
xlabel("t");
ylabel("M(t)")
title("Multipication of input signal and cos(2*pi*f1*t) ")

X3 = fft(Multipication1);
X4 = fftshift(X3);
l = length(Multipication1);
fshift = (-l/2:l/2-1)*(FS/l);


figure
plot(fshift,X4,LineWidth=2);
xlabel("frequency range (Hz)");
ylabel("X1(f)");
title("fourier transform of Multipication of input signal and cos(2*pi*f1*t)");

% HPF applied to the result

f_HPF = (-l/2:l/2-1)*(FS/l);
HPF = zeros(1,l);


for i = 1:l

    if(i<0.296*l)
        HPF(i) = 1;
    end

    if( i>=0.704*l)
        HPF(i) = 1;
    end
end

figure 
plot(f_HPF,HPF,LineWidth=2);
xlabel("frequency range (Hz)")
ylabel("ideal HPF with fc = 9000 ")


Y1 = X4.*HPF;

figure
plot(f_HPF,Y1,LineWidth=2);
xlabel("frequency range (Hz)")
ylabel("Y1(f)")
title("output signal of ideal HPF")


reconstructedSignal = ifftshift(Y1);
reconstructedSignal2 = ifft(reconstructedSignal);

figure
plot(t,reconstructedSignal2,LineWidth=2);
xlabel('t')
ylabel(" r1(t)")
title("output signal of ideal HPF in time_scaling")

Cos2 = cos(2*pi*7000.*t);
Multipication2 = reconstructedSignal2.*Cos2;

X5 = fft(Multipication2);
X6 = fftshift(X5);
l = length(Multipication2);
fshift = (-l/2:l/2-1)*(FS/l);


figure
plot(fshift,X6,LineWidth=2);
xlabel("frequency range (Hz)");
ylabel("X2(f)");
title("fourier transform of Multipication of output signal of ideal HPF and cos(2*pi*f2*t)",FontSize=8);

% spectrum of multipication result by cos(2*pi*f2*t)

f_LPF = (-l/2:l/2-1)*(FS/l);
LPF = zeros(1,l);


for i = 1:l

    if(i>0.296*l && i<=0.704*l)
        LPF(i) = 1;
    end

end

figure 
plot(f_LPF,LPF,LineWidth=2);
xlabel("frequency range (Hz)")
ylabel("ideal LPF with fc = 9000 ")

% final LPF result

Y2 = X6.*LPF;

figure
plot(f_LPF,Y2,LineWidth=2);
xlabel("frequency range (Hz)")
ylabel("Y2(f)")
title("fourier transform of input signal with delta_shifting of spectrum",FontSize=8)


reconstructedSignal3 = ifftshift(Y2);
reconstructedsignal4 = ifft(reconstructedSignal3);
new_sound = abs(reconstructedsignal4);

figure
plot(t,new_sound,LineWidth=2);
xlabel('t')
ylabel(" r2 (t)")
title("input signal with delta_shifting of spectrum  in time_scaling",FontSize=8)

%sound(new_sound,FS);

%audiowrite("sample5_9000_7000.wav",new_sound,FS)