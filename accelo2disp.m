function [LVDTfilt, filtered_disp, filtered_vel, filtered_acc ] ...
          = accelo2disp(t, Ts, Fs, Fcut,alpha, Accmat, Lvdtmat,...
                            lvdtcons, accbiasV, accsensi, filtertype...
                            ,filtermethod,firorder)

%%***********************************************************************%
%*               Acceleration to velocity and displacement              *%
%*             Filters and produces velocity and displacement           *% 
%*                         from acceleration data                       *%
%*                                                                      *%
%* Code author: Preetham Manjunatha                                     *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 08/03/2021                                                     *%
%************************************************************************%
% My apologies for too many input parameters. This function was created
% during my infancy period in learning MATLAB!!!
%************************************************************************%
%
% Usage: [LVDTfilt, filtered_disp, filtered_vel, filtered_acc ] ...
%           = accelo2disp(t, Ts, Fs, Fcut,alpha, Accmat, Lvdtmat,...
%                             lvdtcons, accbiasV, accsensi, filtertype,...
%                             filtermethod,firorder);
%
% Inputs: t             - time vector
%         Ts            - sampling period
%         Fs            - sampling frequency
%         Fcut          - filter Cut-off Frequency scalar for high and low
%                         pass filters. [minimum value maximum value] for bandpass filter.
%         alpha         - scaling constant
%         Accmat        - Accleration matrix for M accelerometers [N x M]
%         Lvdtmat       - LVDT matrix for K accelerometers [N x K] (usually K = 1)
%         lvdtcons      - LVDT calibration constant (M/Volt)
%         accbiasV      - accelerometer bias (Vdc)  
%         accsensi      - accelerometer sensitivity (V/mV/g)
%         filtertype    - lowpass | highpass | bandpass
%         filtermethod  - FFT-fft | FIR-fir
%         firorder      - order of FIR filter (1500 to 2500)
% 
% Outputs: LVDTfilt       - Filtered LVDT output
%          filtered_disp  - Filtered displacements
%          filtered_vel   - Filtered velocities
%          filtered_acc   - Filtered acceleration
% 
% Demo input data is available at https://github.com/preethamam/Acceleration2VelocityandDisplacement

    % Accelerometer to Displcment
    filtered_disp = zeros(size(Accmat));
    filtered_vel = zeros(size(Accmat));
    filtered_acc = zeros(size(Accmat));
    
    for i=1:size(Accmat,2)
        [x_Time, xd_Time, xdd_Time] = ...
            func_acc2disp(t, Accmat(:,i), Fcut, alpha, filtertype,...
                            filtermethod,firorder);
        filtered_disp(:,i) = x_Time;
        filtered_vel(:,i)  = xd_Time;
        filtered_acc(:,i)  = xdd_Time;
    end

    % LVDT Displacement DC Clean
    LVDTfilt = zeros(size(Lvdtmat));
    if ~(isempty(Lvdtmat))
        for i=1:size(Lvdtmat,2)    
        [sigflt, sigfft]...
                = func_DCclean(Lvdtmat(:,i)*lvdtcons(i), Fcut,alpha, Fs, filtertype,filtermethod,...
                                firorder);
        LVDTfilt(:,i)  = sigflt;
        end
    end

    function [x_Time, xd_Time, xdd_Time] = func_acc2disp(t, signalIN, cutoff, alpha, filtertype...
                                                         ,filtermethod,firorder)

    % Function acc2disp performs the filtering in frequency domain and then double integration to
    % obtain the displacment from acceleration

    %% Input Analyzer
    xdd = signalIN; 
    dt = abs(t(5)-t(4));
    Fs = 1/dt;

    %% Double Integration
    % Filter the Acceleration Signal and Set the First k values constant (suppress  frequency content)
    [xdd_Time,~] = ...
              func_DCclean(xdd, cutoff, alpha, Fs, filtertype,filtermethod,firorder);

    %//%------------------------------------------------------------------------%
    % Perform 1st Integration and Filter the Acceleration Signal
    Xd_int = dt * cumtrapz(xdd_Time);

    % Filter the Velocity Signal and Set the First k values constant (suppress  frequency content)
    [xd_Time, ~] = ...
              func_DCclean(Xd_int, cutoff, alpha, Fs, filtertype,filtermethod,firorder);

    %//%------------------------------------------------------------------------%
    % Perform 2st Integration and Filter the Velocity Signal
    X_int = dt * cumtrapz(xd_Time);

    % Filter the displacement Signal and Set the First k values constant (suppress  frequency content)
    % Displacement unfiltered
    [x_TimeNF, ~] = ...
              func_DCclean(X_int, cutoff, alpha, Fs, filtertype,filtermethod,firorder);

    %//%------------------------------------------------------------------------%
    % Final Filtering of Retrived Displacement Data and Set the First k values constant (suppress frequency content)

    %Displacement filtered
    [x_Time, ~] = ...
              func_DCclean(x_TimeNF, cutoff, alpha, Fs, filtertype,filtermethod,firorder);


    end


    function [sigflt, sigfrq] = ...
              func_DCclean(sig, cutoff, alpha, fs, passtype,filtermethod,firorder)

    % This function FUNC_DCCLEAN removes the low frequency content (for the
    % given cutoff value) and the DC offset

    x = length(sig); 
    k = cutoff*(x/fs);
    n = round(k);


    switch filtermethod
        case 'fft'
            SIG_fft = fft(sig, x);
            switch passtype 
                case 'lowpass'
                    SIG_fft(1) = complex(alpha*abs(SIG_fft(n(1))),0);
                    SIG_fft(x) = complex(alpha*abs(SIG_fft(n(1))),0);
                    if (length(n)>1)
                       error('Too many cutoff frequencies') 
                    end
                    for i=n:round((x/2)-1)
                        SIG_fft(i) = alpha*SIG_fft(n);
                        SIG_fft(x-i) = conj(SIG_fft(i));
                    end

                case 'highpass'
                    SIG_fft(1) = complex(alpha*abs(SIG_fft(n(1))),0);
                    if (length(n)>1)
                       error('Too many cutoff frequencies') 
                    end
                    for i=2:n-1
                        SIG_fft(i) = alpha*SIG_fft(n);
                        SIG_fft(x-(i-2)) = conj(SIG_fft(i));
                    end

                otherwise
                    if (length(n)>2)
                       error('Too many cutoff frequencies') 
                    elseif (length(n)<2)
                       error('Too few cutoff frequencies') 
                    end

                    SIG_fft(1) = complex(alpha*abs(SIG_fft(n(1))),0);
                    SIG_fft(x) = complex(alpha*abs(SIG_fft(n(1))),0);
                    for i = n(2):round((x/2)-1)
                        SIG_fft(i) = alpha*SIG_fft(n(2));
                        SIG_fft(x-i) = conj(SIG_fft(i));
                    end

                    for i = n(1):-1:2
                        SIG_fft(i) = alpha*SIG_fft(n(1));
                        SIG_fft(x-i) = conj(SIG_fft(i));
                    end        

            end 
            sigfrq = SIG_fft;
            sigflt = real(ifft(SIG_fft));

        case 'fir'

            switch passtype
                case 'lowpass'
                    if (length(cutoff)>1)
                       error('Too many cutoff frequencies') 
                    end
                    h=fdesign.lowpass('N,Fc',firorder,cutoff,fs);
                    d=design(h); %Lowpass FIR filter
                    sigflt=filtfilt(d.Numerator,1,sig); %zero-phase filtering 
                    sigfrq = abs(fft(sigflt,x));

                case 'highpass'
                    if (length(cutoff)>1)
                       error('Too many cutoff frequencies') 
                    end
                    h=fdesign.highpass('N,Fc',firorder,cutoff,fs);
                    d=design(h); %Highpass FIR filter
                    sigflt=filtfilt(d.Numerator,1,sig); %zero-phase filtering 
                    sigfrq = abs(fft(sigflt,x));              

                otherwise 
                    if (length(cutoff)>2)
                       error('Too many cutoff frequencies') 
                    elseif (length(cutoff)<2)
                       error('Too few cutoff frequencies') 
                    end
                    h=fdesign.bandpass('N,Fc1,Fc2',firorder,cutoff(1),cutoff(2),fs);
                    d=design(h); %Bandpass FIR filter
                    sigflt=filtfilt(d.Numerator,1,sig); %zero-phase filtering 
                    sigfrq = abs(fft(sigflt,x));     

            end      
    end

    end
end
