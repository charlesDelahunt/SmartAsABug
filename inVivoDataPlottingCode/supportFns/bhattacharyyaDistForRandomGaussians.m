% compare a 'true' gaussian to a randomly generated one with N draws, using Bhattacharyya distance
% that is, both distributions are gaussians.
% draw a histogram of the Bhattacharyya distances for the various generated samplings
% inputs:
% numDraws = number of samples to draw from the distribution (default = 100)
% numGen = number of samplings to generate (default = 200)
% sigmaParam = std of gaussian (default = 1)
% muParam = mean of gaussian (default = 0)
% showPlot = boolean, do you show the histogram of Bhat distances

% Disclaimer: This code has not been refactored or otherwise tidied up (!)

% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

function [] = bhattacharyyaDistForRandomGaussians(numDraws, numGen, sigmaParam, muParam, showPlot)

if nargin < 5, showPlot = true; end
if nargin < 4, muParam = 0; end
if nargin < 3, sigmaParam = 1; end
if nargin < 2, numGen = 200; end
if nargin < 1, numDraws = 100; end

% generate samplings and get their params:
allVals = muParam + sigmaParam*randn([numDraws, numGen]); % each col is a sampling
allSigmas = std(allVals); % row vector
allMus = mean(allVals);  % row vector

for i = 1:length(allSigmas),
    vals = allVals(:, i);
    
    % set up the bins for Bhat dist
    lowLim = min(min(vals),-4*max(sigmaParam) ); % assumes max(sigma) is a scalar 
    highLim = max(max(vals),4*max(sigmaParam) );
    
    bins = linspace(lowLim, highLim, 201);
    halfStep = ( bins(2) - bins(1) ) / 2;
    binCenters = bins(1:length(bins) - 1) + halfStep;
    
    % make a gaussian that matches this sampling's sigma and mu:
    gaussian = muParam + exp( -0.5* (binCenters.^2)/sigmaParam^2 ); % make the equiv gaussian
    gaussian = gaussian / sum(gaussian);
    gaussian = [binCenters; gaussian];
    % vals does not need to be made into a dist.
    
    % calc the bhattacharyya distance:
    showBhatPlot = 0;
    Bdistance(i) = bhattacharyya(vals, gaussian,showBhatPlot);
end

if showPlot,
    figure,
    % disp(num2str(Bdistance))
    hist(Bdistance,30)
    title(['hist of Bhat distance. numDraws = ' num2str(numDraws),...
        ', mu = ' num2str(muParam) ', sigma = ',num2str(sigmaParam) ] )
end

% MIT license:
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
% associated documentation files (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
% copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to 
% the following conditions: 
% The above copyright notice and this permission notice shall be included in all copies or substantial 
% portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
% PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
% AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
% WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.