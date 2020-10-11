function [db] = RatioToDb(ratio)
db = 10 .* log10(ratio);
end

