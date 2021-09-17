function [Re, filtered_data] = noiseCovariance(data, long4mean)
% function Re = noiseCovariance(data, long4mean)
%
% Calculate de covariance of the noise of 'data'.
%
% 'long4mean is a measure of desired smoothing. That is, is the number of
% samples taken to calculate the average value around each point.

long4mean = ceil((long4mean-1)/2);
num_data = size(data, 1);
filtered_data = zeros(size(data));

for i = 1:long4mean
	filtered_data(i, :) = mean(data(1:i+long4mean, :));
end
for i = long4mean+1:num_data-long4mean
	filtered_data(i, :) = mean(data(i-long4mean:i+long4mean, :));
end
for i = num_data+1-long4mean:num_data
    filtered_data(i, :) = mean(data(i-long4mean:end, :));
end

Re = cov(data-filtered_data);

% % ##### To chech the function #####
% x = 1:num_data;
% for i = 1:size(data, 2)
% 	figure, subplot(211)
% 	plot(x, data(:, i), 'b-', x, filtered_data(:, i), 'r-')
% 	legend(strcat('data_{', num2str(i), '}'), 'filtered signal')
% 	subplot(212)
% 	plot(x, data(:, i)-filtered_data(:, i))
% 	title(strcat('Noise for data_{', num2str(i), '}'))
% end