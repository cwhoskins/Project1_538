
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

%Loop through all frequencies
for f = 1:length(freq)
    %Loop through all pulses
    current_wavelength = wavelength(f);
    for m = 1:M
        %Calculate Scatterer Range at a given sample
        scatterer_range = initial_pos + ((m / M) * cpi .* velocity);
        %Check if the signal will be recieved before the next pulse
        if(max(scatterer_range) > max_unambiguous_range)
           fprintf('Error: You have to develop logic to handle unambiguous rnages'); 
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
           data(scatterer_sample(s), m, f) = data(scatterer_sample(s), m, f) + scat_sig;
        end
        
        %Generate noise
        noise_power_w = NoisePower(rx_bandwidth, rx_noise);
        noise_power_db = 10*log10(noise_power_w);
        noise = wgn(N,1,noise_power_6,'complex');
        %Add noise to current column
        data(:, m, f) = data(:, m, f) + noise;
    end
end

%Potentially relevant stuff from HW2 that hasn't been implemented
scatterer_delay = (2 * scatterer_range) / c;
scatterer_sample = ceil(scatterer_delay / pulse_width);
fprintf(file_output, 'Return Sample # %d\n', scatterer_sample);
fprintf(file_output, 'Signal Power: %e W\n', signal_power);
signal = zeros(N, 1);
signal(scatterer_sample) = sqrt(signal_power);

noise_power_6 = NoisePower(rx_bandwidth, rx_noise_6);
noise_power_6 = 10*log10(noise_power_6);
fprintf(file_output, 'Noise Power: %f dB\n', noise_power_6);
noise = wgn(N,1,noise_power_6,'complex');

range_profile = signal + noise;
sample_time = 1:N;
plot(sample_time, abs(range_profile));
xlabel('Sample #');
ylabel('Received Signal');
title('Range Profile');
saveas(gcf,'Problem6.png');