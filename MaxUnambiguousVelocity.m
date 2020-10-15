function [max_velocity] = MaxUnambiguousVelocity(wavelength,prf)
max_velocity = (wavelength ./2) .* prf;
end

