function [ratio] = DbToRatio(db)
ratio = 10 .^ (db ./ 10);
end

