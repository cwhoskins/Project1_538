
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
fprintf(file_output, 'Max Unambiguous Range: %g m\n', max_unambiguous_range);
%Calculate Max Unambiguous Velocity
max_velocity = MaxUnambiguousVelocity(wavelength, prf);
%Calculate Sampling Frequency
rx_bandwidth = 1/pulse_width;
fprintf(file_output, 'Sampling Frequency: %g KHz\n', (rx_bandwidth / 1000));
%Calculate Range Resolution
range_resolution = c / (2 * rx_bandwidth);
fprintf(file_output, 'Range Resolution: %g m\n', range_resolution);
%Calculate number of samples per pulse
N = floor(((1/prf) / pulse_width));
%Calculate number of pulses per CPI
M = cpi * prf;
%Create Data Array
data = zeros(N, M, length(freq));
rdm = zeros(N, M, length(freq));
%Calculate Velocity Resolution
vel_resolution = max_velocity ./ M;

%Loop through all frequencies
for f = 1:length(freq)
    %Loop through all pulses
    current_wavelength = wavelength(f);
    
    fprintf(file_output, '\nFrequency: %d MHz\n', (freq(f) / (10^6)));
    fprintf(file_output, 'Max Velocity: %g m/s\n', max_velocity(f));
    fprintf(file_output, 'Velocity Resolution: %g m/s\n', vel_resolution(f));
    
    for m = 1:M
        %Calculate Scatterer Range at a given sample
        scatterer_range = initial_pos + ((m / M) * cpi .* velocity);
        %Check if the signal will be recieved before the next pulse
        
        if(max(scatterer_range) > max_unambiguous_range)
           fprintf('Error: You have to develop logic to handle ambiguous ranges'); 
        end
        %Calculate Received Power from all scatterers
        signal_power = RadarRange(tx_power, tx_gain, rx_gain, current_wavelength, loss, rcs, scatterer_range);
        %Convert power to sampled information
        scat_sig = sqrt(signal_power) .* exp((-1j .* 4 .* pi .* scatterer_range) ./ current_wavelength);
        %Determine when each received sample occured
        scatterer_sample = round(scatterer_range./range_resolution)+1;

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
    
    %Calculate SNR
    pre_snr = sum(signal_power) / noise_power_w;
    pre_snr_db = 10 * log10(pre_snr);
    fprintf(file_output, 'Pre-FFT SNR: %g dB\n', pre_snr_db);
    %Calculate expected SNR after Doppler Processing
    expected_snr = pre_snr * M;
    expected_snr_db = 10 * log10(expected_snr);
    fprintf(file_output, 'Predicted Post-FFT SNR: %g dB\n', expected_snr_db);
    
    %Perform FFT across every row
    for n = 1:N
        rdm(n,:,f) = fftshift(abs(fft(data(n, :, f))));
    end
    
    %Calculate noise and signal power
    %Estimate maximum noise value
    noise_limit = max(abs(noise)) * M * 1.2;
    post_signal_pwr = 0;
    post_noise_pwr = 0;
    num_samples = N * M;
    for n=1:N
       for m=1:M
          if(abs(rdm(n, m, f)) <= noise_limit) %Current bin is noise
              post_noise_pwr = post_noise_pwr + abs(rdm(n, m, f))^2;
          else %Current bin is signal
            post_signal_pwr = post_signal_pwr + abs(rdm(n, m, f))^2;
          end
       end
    end
    post_noise_pwr = post_noise_pwr / num_samples;
    post_snr = post_signal_pwr / post_noise_pwr;
    post_snr_db = 10 * log10(post_snr);
    fprintf(file_output, 'Measured Post-FFT SNR: %g dB\n', post_snr_db);
    
    %Convert Range Doppler Map to dB
    rdm(:, :, f) = 10 .* log10(rdm(:, :, f));
    
    %Generate Surface Plot
    plot_title = sprintf('RDM - %.2e MHz', (freq(f) / (10^6)));
    plot_file = sprintf('RDM-%d.png', f);
    axis_range = 0:range_resolution:(max_unambiguous_range-range_resolution);
    axis_velocity = (0:vel_resolution(f):(max_velocity(f)-vel_resolution(f)))-(max_velocity(f)/2);
    surf(axis_velocity, axis_range, rdm(:, :, f));
    xlabel('Velocity (m/s)');
    ylabel('Range (m)');
    zlabel('Signal (dB)');
    title(plot_title);
    saveas(gcf,plot_file);
end