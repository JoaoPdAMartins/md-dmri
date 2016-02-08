function t = ts_1x15_to_6x6(N)
% function t = ts_1x15_to_6x6(N)


% % t = [...
% 1%     N(:,1) .* N(:,1) .* N(:,1) .* N(:,1) * sqrt(1) ...
% 2%     N(:,1) .* N(:,1) .* N(:,1) .* N(:,2) * sqrt(4) ...
% 3%     N(:,1) .* N(:,1) .* N(:,1) .* N(:,3) * sqrt(4) ...
% 4%     N(:,1) .* N(:,1) .* N(:,2) .* N(:,2) * sqrt(6) ...
% 5%     N(:,1) .* N(:,1) .* N(:,2) .* N(:,3) * sqrt(12) ...
% 6%     N(:,1) .* N(:,1) .* N(:,3) .* N(:,3) * sqrt(6) ...
% 7%     N(:,1) .* N(:,2) .* N(:,2) .* N(:,2) * sqrt(4) ...
% 8%     N(:,1) .* N(:,2) .* N(:,2) .* N(:,3) * sqrt(12) ...
% 9%     N(:,1) .* N(:,2) .* N(:,3) .* N(:,3) * sqrt(12) ...
% 10%     N(:,1) .* N(:,3) .* N(:,3) .* N(:,3) * sqrt(4) ...
% 11%     N(:,2) .* N(:,2) .* N(:,2) .* N(:,2) * sqrt(1) ...
% 12%     N(:,2) .* N(:,2) .* N(:,2) .* N(:,3) * sqrt(4) ...
% 13%     N(:,2) .* N(:,2) .* N(:,3) .* N(:,3) * sqrt(6) ...
% 14%     N(:,2) .* N(:,3) .* N(:,3) .* N(:,3) * sqrt(4) ...
% 15%     N(:,3) .* N(:,3) .* N(:,3) .* N(:,3) * sqrt(1) ];

if (size(N,1) ~= 1), error('stop'); end

N = N ./ sqrt([...
    1 4 4 6 12 6 4 12 12 4 1 4 6 4 1]);

t = N([...
    1  4  6   5  3 2;
    4 11 13  12  8 7;
    6 13 15  14 10 9;
    ...
    5 12 14  13  9 8;
    3  8 10   9  6 5;
    2  7  9   8  5 4]) .* ...
    sqrt([
    1 1 1 2 2 2;
    1 1 1 2 2 2;
    1 1 1 2 2 2;
    2 2 2 4 4 4;
    2 2 2 4 4 4;
    2 2 2 4 4 4]);
    