
%Given Values
tx_power = 100;
tx_gain = 1000;
rx_gain = 1000;
freq = [(500 * 10^6) (1 * 10^9) (10 * 10^9) (28 * 10^9)];
wavelength = (3*10^8) ./ freq;
rcs = 10;
loss = 10^(5/10);
pulse_width = 3.333 * 10^-6;
rx_noise = 10^(5/10);
temperature = 290;
prf = 20 * 10^3;
cpi = 30 * 10^-3;
c = 3 * 10 ^ 8;
initial_pos = [2500 2500];
velocity = [10 0];

file_output = fopen('output.txt', 'w+');

%Calculate Max Unambigous Range
max_unambiguous_range = MaxUnambiguousRange(prf);
fprintf(file_output, 'Max Unambiguous Range: %f\n', max_unambiguous_range);
%Calculate Sampling Frequency
rx_bandwidth = 1/pulse_width;
fprintf(file_output, 'Sampling Frequency: %f\n', rx_bandwidth);
%Calculate number of samples per pulse
N = floor(((1/prf) / pulse_width) - 1);
%Calculate number of pulses per CPI
M = cpi * prf;
%Create Data Array
data = zeros(N, M, length(freq));
rdm = zeros(N, M, length(freq));

%Loop through all frequencies
for f = 1:length(freq)
    %Loop through all pulses
    current_wavelength = wavelength(f);
    for m = 1:M
        %Calculate Scatterer Range at a given sample
        scatterer_range = initial_pos + ((m / M) * cpi .* velocity);
        %Check if the signal will be recieved before the next pulse
        if(max(scatterer_range) > max_unambiguous_range)
           fprintf('Error: You have to develop logic to handle unambiguous ranges'); 
        end
        %Calculate Received Power from all scatterers
        signal_power = RadarRange(tx_power, tx_gain, rx_gain, current_wavelength, loss, rcs, scatterer_range);
        %Convert power to sampled information
        scat_sig = sqrt(signal_power) .* exp((-1j .* 4 .* pi .* scatterer_range) ./ current_wavelength);
        %Determine when each received sample occured
        scatterer_delay = (2 .* scatterer_range) ./ c;
        scatterer_sample = ceil(scatterer_delay ./ pulse_width);

        %Loop through the scatterers
        for s = 1:length(scatterer_sample)
            %Add received signal to the pre-determined spots
           data(scatterer_sample(s), m, f) = data(scatterer_sample(s), m, f) + scat_sig(s);
        end
        
        %Generate noise
        noise_power_w = NoisePower(rx_bandwidth, rx_noise);
        noise_power_db = 10*log10(noise_power_w);
        noise = wgn(N,1,noise_power_db,'complex');
        %Add noise to current column
        data(:, m, f) = data(:, m, f) + noise;
    end
    
    for n = 1:N
        rdm(n,:,f) = fftshift(abs(fft(data(n, :, f))));
    end
    rdm(:, :, f) = 10 .* log10(rdm(:, :, f));
    axis_range = ((1:N) ./ N) .* max_unambiguous_range;
    axis_velocity = ((1:M) ./ M) .* MaxUnambiguousVelocity(current_wavelength, prf);
    surf(axis_velocity, axis_range, rdm(:, :, f));
end