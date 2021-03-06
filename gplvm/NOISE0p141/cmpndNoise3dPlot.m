function gX = cmpndNoise3dPlot(noise, X)

% CMPNDNOISE3DPLOT Draws a 3D or contour plot for the CMPND noise model.
%
%	Description:
%
%	H = CMPNDNOISE3DPLOT(NOISE, PLOTTYPE, X, Y, MU, VARSIGMA, ...) draws
%	a 3D or contour plot for the compound noise model.
%	 Returns:
%	  H - the gradients of the diagonal with respect to each element of
%	   X. The returned matrix has the same dimensions as X.
%	 Arguments:
%	  NOISE - the noise structure for which the plot is required.
%	  PLOTTYPE - string containing the name of the plotting function
%	   (for example mesh, contour).
%	  X - the input X data in the form of a 'mesh' matrix.
%	  Y - the input Y data in the form of a 'mesh' matrix.
%	  MU - the input mean in the form of a 'mesh' matrix.
%	  VARSIGMA - the input variance in the form of a 'mesh' matrix.
%	  ... - optional additional arguments for the given plot type.
%	
%
%	See also
%	CMPNDNOISEPARAMINIT, NOISE3DPLOT, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence
% 	cmpndNoise3dPlot.m CVS version 1.1
% 	cmpndNoise3dPlot.m SVN version 29
% 	last update 2007-11-03T14:29:08.000000Z

