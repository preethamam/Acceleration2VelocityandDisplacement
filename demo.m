%% Start parameters
%--------------------------------------------------------------------------
clear; close all; clc;
Start = tic;

%% Inputs
%--------------------------------------------------------------------------
%
% Data files
data = load('Output_Data_Final_Test1.lvm');

% Filter Cut-off Frequency and alpha
Fcut = [0.5];
alpha = 0.0001;
filtertype = 'highpass';      % lowpass | highpass | bandpass
filtermethod = 'fft';         % FFT-fft | FIR-fir
firorder = 2500;

% Sensor Constants
lvdtcons = 0.52725 /100;          % M/Volt - LVDT constant
accbiasV = 2.5;                   % Vdc -  accelerometer bias
accsensi = 1015/1000;             % V/mV/g -  accelerometer sensitivity
Lvdtmat = [];

%% Processed inputs
%--------------------------------------------------------------------------
time   = data(:,1);         % Time vector
accval_g = data(:,2);       % Acceleration in g
accval = data(:,2)*9.81;    % Acceleration in m/s^2
L = size(data,1);           % Length of signal
Fs = 1/(time(2)-time(1));   % Sampling frequency                    
Ts = 1/Fs;                  % Sampling period       

%% Displacement, Velocity and Acceleration
%--------------------------------------------------------------------------
[~, ~, ~, filtered_acc_g] ...
          = accelo2disp(time,Ts, Fs, Fcut,alpha, accval_g, Lvdtmat,...
                            lvdtcons, accbiasV, accsensi, filtertype...
                            ,filtermethod,firorder);
                        
[LVDTfilt, filtered_disp, filtered_vel, filtered_acc] ...
          = accelo2disp(time,Ts, Fs, Fcut,alpha, accval, Lvdtmat,...
                            lvdtcons, accbiasV, accsensi, filtertype...
                            ,filtermethod,firorder);

%% Compute the frequency
%--------------------------------------------------------------------------
% Compute the Fourier transform of the signal.
Y = fft(filtered_acc);

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 
% based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain f and plot the single-sided amplitude spectrum P1. 
% The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. 
% On average, longer signals produce better frequency approximations.
f = Fs*(0:(L/2))/L;

%% Plot Acceleration vs.time
%--------------------------------------------------------------------------
figure;
subplot(5,2,1)
ax1 = plot(time, accval_g, 'LineWidth', 1);
grid on
xlabel('Time ($sec.$)','Interpreter', 'latex');
ylabel('Acceleration ($g$)','Interpreter', 'latex');
title('Acceleration vs.time (with DC bias)')

ax2 = subplot(5,2,2);
plot(time, filtered_acc_g, 'LineWidth', 1)
grid on
xlabel('Time ($sec.$)','Interpreter', 'latex');
ylabel('Acceleration ($g$)','Interpreter', 'latex');
title('Acceleration vs.time (DC bias removed)')

ax3 = subplot(5,2,3:4);
plot(time,filtered_acc, 'LineWidth', 1)
grid on
xlabel('Time ($sec.$)','Interpreter', 'latex');
ylabel('Acceleration ($\frac{m}{s^2}$)','Interpreter', 'latex');
title('Acceleration vs.time (DC bias removed)')

ax4 = subplot(5,2,5:6);
plot(time,filtered_vel, 'LineWidth', 1)
grid on
xlabel('Time ($sec.$)','Interpreter', 'latex');
ylabel('Velocity ($\frac{m}{s}$)','Interpreter', 'latex');
title('Velocity vs.time')

ax5 = subplot(5,2,7:8);
plot(time,filtered_disp, 'LineWidth', 1)
grid on
xlabel('Time ($sec.$)','Interpreter', 'latex');
ylabel('Displacement ($m$)','Interpreter', 'latex');
title('Displacement vs.time')

ax6 = subplot(5,2,9:10);
plot(f,P1, 'LineWidth', 1)
grid on
title('Frequency Spectrum of Acceleration')
xlabel('f(Hz)')
ylabel('|P1(f)|')

%% Print the results
%--------------------------------------------------------------------------
% Extract the 
[M,I] = max(P1);

% Frequency of the signal
fprintf('--------------------------------------------------------------------------\n');
fprintf('Frequency of the signal\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('Frequncy is: %.3fHz \n', f(I));
fprintf('\n\n\n');

% Peaks of the signals
fprintf('--------------------------------------------------------------------------\n');
fprintf('Peaks of the signals\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('Peak acceleration is: %.3f (g) \n', max(abs(filtered_acc))/9.81)
fprintf('Peak acceleration is: %.3f (m/s^2) \n', max(abs(filtered_acc)))
fprintf('Peak velocity is: %.3f (m/s) \n', max(abs(filtered_vel)))
fprintf('Peak displactment is: %.3f (m) \n', max(abs(filtered_disp)))
fprintf('\n\n\n');

%% End parameters
%--------------------------------------------------------------------------
Runtime = toc(Start);
fprintf('--------------------------------------------------------------------------\n');
fprintf('Computational Run time\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('Run time is: %.3f (sec.) \n', Runtime)