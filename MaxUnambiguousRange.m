function [max_range] = MaxUnambiguousRange(prf)
c = 3 * 10^8;
max_range = c ./ (2 .* prf);
end

