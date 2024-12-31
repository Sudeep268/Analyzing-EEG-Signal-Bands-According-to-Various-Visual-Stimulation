%Analyzing EEG signal bands according to various visual stimulation
%1915012 & 1915013
clc; clear all;

%importing data

yellow  = readtable('E:\Study docu\3-1\BSP Lab\Project\DATA\Final\eeg1_f_yellow_1-L03.txt');
yellow  = table2array(yellow);

red  = readtable('E:\Study docu\3-1\BSP Lab\Project\DATA\Final\eeg1_f_red_2-L03.txt');
red  = table2array(red);

blue  = readtable('E:\Study docu\3-1\BSP Lab\Project\DATA\Final\eeg1_f_blue_3-L03.txt');
blue  = table2array(blue);

green  = readtable('E:\Study docu\3-1\BSP Lab\Project\DATA\Final\eeg1_f_green_4-L03.txt');
green  = table2array(green);

black  = readtable('E:\Study docu\3-1\BSP Lab\Project\DATA\Final\eeg1_f_black_5-L03.txt');
black  = table2array(black);

text  = readtable('E:\Study docu\3-1\BSP Lab\Project\DATA\Final\eeg1_f_screenread_6-L03.txt');
text  = table2array(text);

%Segmentation of the data

%EEG
eeg(:,1)=yellow(:,1);
eeg(:,2)=red(:,1);
eeg(:,3)=blue(:,1);
eeg(:,4)=green(:,1);
eeg(:,5)=black(:,1);
eeg(:,6)=text(:,1);

%ALPHA
alpha(:,1)=yellow(:,2);
alpha(:,2)=red(:,2);
alpha(:,3)=blue(:,2);
alpha(:,4)=green(:,2);
alpha(:,5)=black(:,2);
alpha(:,6)=text(:,2);

%BETA
beta(:,1)=yellow(:,3);
beta(:,2)=red(:,3);
beta(:,3)=blue(:,3);
beta(:,4)=green(:,3);
beta(:,5)=black(:,3);
beta(:,6)=text(:,3);

%DELTA
delta(:,1)=yellow(:,4);
delta(:,2)=red(:,4);
delta(:,3)=blue(:,4);
delta(:,4)=green(:,4);
delta(:,5)=black(:,4);
delta(:,6)=text(:,4);

%THETA
theta(:,1)=yellow(:,5);
theta(:,2)=red(:,5);
theta(:,3)=blue(:,5);
theta(:,4)=green(:,5);
theta(:,5)=black(:,5);
theta(:,6)=text(:,5);

% Designing 50 Hz notch filter 
wo=(50/(400/2));
bw=wo/35;
[b,a]=iirnotch(wo,bw);

% Cancelling the 50Hz component from data
for i = 1:6
   eeg(:,i) = filter(b,a,eeg(:,i));
   alpha(:,i) = filter(b,a,alpha(:,i));
   beta(:,i) = filter(b,a,beta(:,i));
   delta(:,i) = filter(b,a,delta(:,i));
   theta(:,i) = filter(b,a,theta(:,i));
end
% Designing a lowpass filter of 30Hz
fc=30;
fs=400;
[d,c]=butter(10,fc/(fs/2));

%Plotting the frequency response of the filter
figure(1)
freqz(b,a);
figure(2) 
freqz(d,c);

% Applying the lowpass filter
for i = 1:6
   eeg(:,i) = filter(d,c,eeg(:,i));
   alpha(:,i) = filter(d,c,alpha(:,i));
   beta(:,i) = filter(d,c,beta(:,i));
   delta(:,i) = filter(d,c,delta(:,i));
   theta(:,i) = filter(d,c,theta(:,i));
end

%Finding the magnitude response 
for i=1:6
    fft_eeg(:,i)= fft(eeg(:,i));
    fft_eeg(:,i)= abs(fft_eeg(:,i));
    fft_alpha(:,i)= fft(alpha(:,i));
    fft_alpha(:,i)= abs(fft_alpha(:,i));
    fft_beta(:,i)= fft(beta(:,i));
    fft_beta(:,i)= abs(fft_beta(:,i));
    fft_delta(:,i)= fft(delta(:,i));
    fft_delta(:,i)= abs(fft_delta(:,i));
    fft_theta(:,i)= fft(theta(:,i));
    fft_theta(:,i)= abs(fft_theta(:,i));
end

%Trimming the values
rowstokeep = 1:5085;
for i=1:6
    fft_eeg_f(:,i) = fft_eeg(rowstokeep,i);  
    fft_alpha_f(:,i) = fft_alpha(rowstokeep,i);
    fft_beta_f(:,i) = fft_beta(rowstokeep,i);
    fft_delta_f(:,i) = fft_delta(rowstokeep,i);
    fft_theta_f(:,i) = fft_theta(rowstokeep,i);
end

%Creating the frequency axis
axis_eeg = linspace(.5,30, 5085);
axis_alpha = linspace(7.5,12,5085);
axis_beta = linspace(12,30,5085);
axis_delta = linspace(0.5,14,5085);
axis_theta = linspace(4,7.5,5085);
axis_eeg = axis_eeg';
axis_alpha = axis_alpha';
axis_beta = axis_beta';
axis_delta = axis_delta';
axis_theta = axis_theta';

% Plotting the freqency response
figure(3)
subplot(321)
plot(axis_eeg,fft_eeg_f(:,1));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the EEG signal for yellow color');

subplot(322)
plot(axis_eeg,fft_eeg_f(:,2));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the EEG signal for red color');

subplot(323)
plot(axis_eeg,fft_eeg_f(:,3));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the EEG signal for blue color');

subplot(324)
plot(axis_eeg,fft_eeg_f(:,4));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the EEG signal for green color');

subplot(325)
plot(axis_eeg,fft_eeg_f(:,5));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the EEG signal for black color');

subplot(326)
plot(axis_eeg,fft_eeg_f(:,6));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the EEG signal for text read');

figure(4)
subplot(321)
plot(axis_alpha,fft_alpha_f(:,1));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Alpha wave for yellow color');

subplot(322)
plot(axis_alpha,fft_alpha_f(:,2));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Alpha wave for red color');

subplot(323)
plot(axis_alpha,fft_alpha_f(:,3));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Alpha wave for blue color');

subplot(324)
plot(axis_alpha,fft_alpha_f(:,4));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Alpha wave for green color');

subplot(325)
plot(axis_alpha,fft_alpha_f(:,5));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Alpha wave for black color');

subplot(326)
plot(axis_alpha,fft_alpha_f(:,6));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Alpha wave for text read');

figure(5)
subplot(321)
plot(axis_beta,fft_beta_f(:,1));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Beta wave for yellow color');

subplot(322)
plot(axis_beta,fft_beta_f(:,2));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Beta wave for red color');

subplot(323)
plot(axis_beta,fft_beta_f(:,3));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Beta wave for blue color');

subplot(324)
plot(axis_beta,fft_beta_f(:,4));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Beta wave for green color');

subplot(325)
plot(axis_beta,fft_beta_f(:,5));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Beta wave for black color');

subplot(326)
plot(axis_beta,fft_beta_f(:,6));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Beta wave for text read');

figure(6)
subplot(321)
plot(axis_delta,fft_delta_f(:,1));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Delta wave for yellow color');

subplot(322)
plot(axis_delta,fft_delta_f(:,2));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Delta wave for red color');

subplot(323)
plot(axis_delta,fft_delta_f(:,3));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Delta wave for blue color');

subplot(324)
plot(axis_delta,fft_delta_f(:,4));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Delta wave for green color');

subplot(325)
plot(axis_delta,fft_delta_f(:,5));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Delta wave for black color');

subplot(326)
plot(axis_delta,fft_delta_f(:,6));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Delta wave for text read');

figure(7)
subplot(321)
plot(axis_theta,fft_theta_f(:,1));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Theta wave for yellow color');

subplot(322)
plot(axis_theta,fft_theta_f(:,2));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Theta wave for red color');

subplot(323)
plot(axis_theta,fft_theta_f(:,3));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Theta wave for blue color');

subplot(324)
plot(axis_theta,fft_theta_f(:,4));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Theta wave for green color');

subplot(325)
plot(axis_theta,fft_theta_f(:,5));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Theta wave for black color');

subplot(326)
plot(axis_theta,fft_theta_f(:,6));
xlabel('Frequency (Hz)'); grid on;
ylabel('magnitude');
title('Magnitude Responses of the Theta wave for text read');

%Creating the frequency axis
axis_eeg = linspace(.5,30, 501);
axis_alpha = linspace(7.5,12,501);
axis_beta = linspace(12,30,501);
axis_delta = linspace(0.5,14,501);
axis_theta = linspace(4,7.5,501);
axis_eeg = axis_eeg';
axis_alpha = axis_alpha';
axis_beta = axis_beta';
axis_delta = axis_delta';
axis_theta = axis_theta';

figure(8)
[pxx, f] = pmtm(eeg(:,1),15,1000);
subplot(3,2,1)
plot(axis_eeg, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of EEG for yellow color'); 

[pxx, f] = pmtm(eeg(:,2),15,1000);
subplot(3,2,2)
plot(axis_eeg, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of EEG for red color'); 

[pxx, f] = pmtm(eeg(:,3),15,1000);
subplot(3,2,3)
plot(axis_eeg, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of EEG for blue color');

[pxx, f] = pmtm(eeg(:,4),15,1000);
subplot(3,2,4)
plot(axis_eeg, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of EEG for green color'); 

[pxx, f] = pmtm(eeg(:,5),15,1000);
subplot(3,2,5)
plot(axis_eeg, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of EEG for black color'); 

[pxx, f] = pmtm(eeg(:,6),15,1000);
subplot(3,2,6)
plot(axis_eeg, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of EEG for text read'); 

figure(9)
[pxx, f] = pmtm(alpha(:,1),15,1000);
subplot(3,2,1)
plot(axis_alpha, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Alpha wave for yellow color'); 

[pxx, f] = pmtm(alpha(:,2),15,1000);
subplot(3,2,2)
plot(axis_alpha, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Alpha wave for red color'); 

[pxx, f] = pmtm(alpha(:,3),15,1000);
subplot(3,2,3)
plot(axis_alpha, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Alpha wave for blue color');

[pxx, f] = pmtm(alpha(:,4),15,1000);
subplot(3,2,4)
plot(axis_alpha, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Alpha wave for green color'); 

[pxx, f] = pmtm(alpha(:,5),15,1000);
subplot(3,2,5)
plot(axis_alpha, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Alpha wave for black color'); 

[pxx, f] = pmtm(alpha(:,6),15,1000);
subplot(3,2,6)
plot(axis_alpha, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Alpha wave for text read'); 


figure(10)
[pxx, f] = pmtm(beta(:,1),15,1000);
subplot(3,2,1)
plot(axis_beta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Beta wave for yellow color'); 

[pxx, f] = pmtm(beta(:,2),15,1000);
subplot(3,2,2)
plot(axis_beta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Beta wave for red color'); 

[pxx, f] = pmtm(beta(:,3),15,1000);
subplot(3,2,3)
plot(axis_beta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Beta wave for blue color');

[pxx, f] = pmtm(beta(:,4),15,1000);
subplot(3,2,4)
plot(axis_beta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Beta wave for green color'); 

[pxx, f] = pmtm(beta(:,5),15,1000);
subplot(3,2,5)
plot(axis_beta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Beta wave for black color'); 

[pxx, f] = pmtm(beta(:,6),15,1000);
subplot(3,2,6)
plot(axis_beta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Beta wave for text read'); 

figure(11)
[pxx, f] = pmtm(delta(:,1),15,1000);
subplot(3,2,1)
plot(axis_delta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Delta wave for yellow color'); 

[pxx, f] = pmtm(delta(:,2),15,1000);
subplot(3,2,2)
plot(axis_delta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Delta wave for red color'); 

[pxx, f] = pmtm(delta(:,3),15,1000);
subplot(3,2,3)
plot(axis_delta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Delta wave for blue color');

[pxx, f] = pmtm(delta(:,4),15,1000);
subplot(3,2,4)
plot(axis_delta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Delta wave for green color'); 

[pxx, f] = pmtm(delta(:,5),15,1000);
subplot(3,2,5)
plot(axis_delta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Delta wave for black color'); 

[pxx, f] = pmtm(delta(:,6),15,1000);
subplot(3,2,6)
plot(axis_delta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Delta wave for text read'); 

figure(12)
[pxx, f] = pmtm(theta(:,1),15,1000);
subplot(3,2,1)
plot(axis_theta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Theta wave for yellow color'); 

[pxx, f] = pmtm(theta(:,2),15,1000);
subplot(3,2,2)
plot(axis_theta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Theta wave for red color'); 

[pxx, f] = pmtm(theta(:,3),15,1000);
subplot(3,2,3)
plot(axis_theta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Theta wave for blue color');

[pxx, f] = pmtm(theta(:,4),15,1000);
subplot(3,2,4)
plot(axis_theta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Theta wave for green color'); 

[pxx, f] = pmtm(theta(:,5),15,1000);
subplot(3,2,5)
plot(axis_theta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Theta wave for black color'); 

[pxx, f] = pmtm(theta(:,6),15,1000);
subplot(3,2,6)
plot(axis_theta, 10*log10(pxx));
xlabel('Frequency (Hz)'); grid on;
ylabel('Power Spectral Density (dB/Hz)');
title('Estimation of PSD of Theta wave for text read'); 