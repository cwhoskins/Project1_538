function [noise_power] = NoisePower(bandwidth, rx_noise)

k = 1.3807 * 10^-23;
temp = 290;
noise_power = k .* temp .* bandwidth .* rx_noise;
end

