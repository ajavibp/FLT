function newData = cutData(data, firstPoint, lastPoint, refillPoints, columns, newFile, format)
% 
% function newData = cutData(data, firstPoint, lastPoint, refillPoints, columns)
% 
% function newData = cutData(data, firstPoint, lastPoint, refillPoints, columns, newFile)
% 
% function newData = cutData(data, firstPoint, lastPoint, refillPoints, columns, newFile, format)
% 
% newData -> Corrected data vector.
% 
% data -> Original data (a sample for row).
% 
% firstPoint -> Data before 'firstPoint' will be removed.
% 
% lastPoint -> Data after 'lastPoint' will be removed.
% 
% refillPoints -> Pairs of indices to fill with linear progression.
% 
% columns -> Columns to refill with linear progression.
% 
% newFile -> Nome of file to save new data in ASCII mode.
% 
% format -> Format to save data (ASCII by default). See save help.
% 
% Example:
% [x,data]=meshgrid(1:5, 1:10);
% data(5:6, 3) = 0;
% data(6:7, 5) = 0;
% newData = cutData(data, 2, 8, [4, 7; 5, 8], [3; 5])

dr = diff(refillPoints') + 1;

for c = columns
	for r = 1:size(refillPoints, 1)
		data(refillPoints(r,1):refillPoints(r, 2), c) = linspace(data(refillPoints(r, 1), c), data(refillPoints(r, 2), c), dr(r))';
	end
end

newData = data(firstPoint:lastPoint, :);

if exist('newFile', 'var')
	if ~exist('format', 'var')
		format = '-ascii';
	end
	save(newFile, 'newData', format)
end
