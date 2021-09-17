function output = vectordelay(input, delays)
% output = VECTORDELAY(input, delay)
%
% Obtain a delayed signal across 'delay' values.
% 
% output = [ input(k-delays(1)), input(k-delays(2)), ..., input(k-delays(end)) ]
% 
% input  -> Input signal, in column form.
% output -> Delayed signal for Anfis.
% delays -> Delay values.

if all(size(delays) > 1) % Ensure 'delays' is a vector
    error('''delays'' must be a vector.');
end
if all(size(input) > 1) || all(size(input) == 1) % Ensure 'input' is a vector
    error('''input'' must be a vector.');
end
if size(input, 2) >1 % Ensure 'input' is a column vector
    input = input';
end

size_data = length(input);
output = zeros(size_data, length(delays));
for column = 1:length(delays)
    d = delays(column);
%     output(:, column) = [ones(d, 1)*input(1); input(1:size_data-d)];
	output(:, column) = [ones(d, 1)*0; input(1:size_data-d)];
end
