function [rx_power] = RadarRange(tx_power, tx_gain, rx_gain, wavelength, loss, rcs, range)
    numerator = (tx_power .* tx_gain .* rx_gain .* (wavelength.^2) * rcs);
    denominator = ((4*pi)^3 .* loss .* (range .^ 4));
    rx_power = numerator ./ denominator;
end

