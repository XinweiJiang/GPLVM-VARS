function h = ncnmNoisePointPlot(noise, X, y, ...
                              fontName, fontSize, ...
                              markerSize, lineWidth);

% NCNMNOISEPOINTPLOT Plot the data-points for null category noise model.
%
%	Description:
%	h = ncnmNoisePointPlot(noise, X, y, ...
%                              fontName, fontSize, ...
%                              markerSize, lineWidth);
%% 	ncnmNoisePointPlot.m CVS version 1.1
% 	ncnmNoisePointPlot.m SVN version 29
% 	last update 2007-11-03T14:29:10.000000Z


pointsNeg = plot(X(find(y(:, 1)==-1), 1), ...
		 X(find(y(:, 1)==-1), 2), ...
		 'gx', 'erasemode', 'xor', ...
		 'markersize', markerSize+2, ...
		 'linewidth', lineWidth);
hold on
pointsPos = plot(X(find(y(:, 1)==1), 1), ...
		 X(find(y(:, 1)==1), 2), 'ro', ...
		 'erasemode', 'xor', ...
		 'markersize', markerSize, ...
		 'linewidth', lineWidth);

pointUn = plot(X(:, 1), ...
	       X(:, 2), 'm.', ...
	       'erasemode', 'xor', 'markersize', markerSize/2);

h = [pointsNeg; pointsPos; pointUn];