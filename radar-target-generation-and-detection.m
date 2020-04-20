clear all
clc;

%% Radar Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity
% remains contant
target_range = 50.0;
target_velocity = 5; %m/s

c = 3e8; %m/s
max_range = 200; % m/s
range_resolution = 1; %m
max_velocity = 100; %m/s


%% FMCW Waveform Generation:
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.


%Operating carrier frequency of Radar
fc= 77e9;             %carrier freq


%Bsweep
Bsweep = c / (2*range_resolution);
trip_time = (2 * max_range)/c;
Tchirp = 5.5 * trip_time;
slope = Bsweep/Tchirp;

%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation.
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp.
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));



%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.

for i=1:length(t)
    %For each time stamp update the Range of the Target for constant velocity.
    r_t(i) = target_range + (target_velocity*t(i));
    td(i) = 2*r_t(i)/c; %Time delay

    %For each time sample we need update the transmitted and
    %received signal.
    Tx(i) = cos(2*pi*(fc*t(i) + (slope/2)*t(i)*t(i)));
    time_difference = t(i) - td(i);
    Rx(i)  = cos(2*pi*(fc*(time_difference) + (slope/2)*(time_difference^2)));

    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i) .* Rx(i);
end

%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.

Mix = reshape(Mix, [Nr,Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.

fft_1d = fft(Mix,Nr);

fft_1d = fft_1d./Nr;

% Take the absolute value of FFT output

fft_1d = abs(fft_1d);

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
fft_1d = fft_1d(1:(Nr/2));

%plotting the range
%figure ('Name','Range from First FFT')
subplot(2,1,1)

% plot FFT output

plot(fft_1d);
axis ([0 100 0 0.5]);
xlabel('Frequency [Hz]');
ylabel('Range [m]');
title("Target Range calculated through 1D FFT");




%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.

Tr = 5;
Td = 5;

%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation

Gr = 2;
Gd = 2;

% offset the threshold by SNR value in dB

offset = 5;

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);
thresholded_cells = zeros(Nr,Nd);


%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
n_trn_guard_cells_range = Gr + Tr;
n_trn_guard_cells_doppler = Gd + Td;

n_training_cells = (2 * Tr + 2 * Gr +1) * (2 * Td + 2 * Gd +1) - (2*Gr+1)*(2*Gd+1);

thresholded_RDM = RDM;

for range_ix = (Gr+Tr+1) : (Nr/2 - (Gr+Tr))
    for doppler_ix = (Gd+Td+1) : (Nd - (Gd+Td))
        noise_level = zeros(1,1);
        for row = range_ix - (Gr+Tr): range_ix + (Gr + Tr)
            for column = doppler_ix -(Gd + Td): doppler_ix + (Gd + Td)
                if(abs(range_ix - row) > Gr || abs(doppler_ix - column) > Gd)
                    noise_level = noise_level + db2pow(RDM(row,column));
                end
            end
        end
        noise_level = noise_level/n_training_cells;
        noise_level = pow2db(noise_level);
        noise_level = noise_level + offset;

        if(RDM(range_ix,doppler_ix) > noise_level)
            thresholded_RDM(range_ix,doppler_ix) = 1;
        else
            thresholded_RDM(range_ix,doppler_ix) = 0;
        end
    end
end

% The process above will generate a thresholded block, which is smaller
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0.

non_zero_indexes = thresholded_RDM ~= 0
non_one_indexes = thresholded_RDM ~= 1

non_thresholded_mask = non_one_indexes & non_zero_indexes;

thresholded_RDM(non_thresholded_mask) = 0;

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,thresholded_RDM);
title("Target Velocity calculated through 2D FFT");
colorbar;
