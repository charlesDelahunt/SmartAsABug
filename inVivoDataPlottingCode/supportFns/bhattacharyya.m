% calculate Bhattacharyya distance between two distributions:
% inputs: 
% d1 and d2: 2 x N matrices or vectors. Row 1 = values, row 2 = counts at that values
%            if row 2 has negative values, nonsense will output.
% if d1 (or d2) is a vector, a row 2 = [1 1 1 ... 1 etc] will be added.
% d1 and d2 don't need to be normalized or binned.
% they will be binned according to the collective max and min of the values.
% also, showPlot: boolean, plot the 2 distributions and bhat dist, default = 0

% Disclaimer: This code has not been refactored or otherwise tidied up (!)

% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

function [bhat] = bhattacharyya(d1, d2, showPlot)

if nargin < 3, showPlot = 0; end

numBins = 200;

% put inputs in proper form:
vectorFlag = zeros(1,2);

% transpose if needed: (assumes > 2 samples in the sampling)
if size(d1,1) > 2, 
    d1 = d1'; 
end
% add counts as row 2 if needed:
if size(d1,1) == 1, 
    vectorFlag(1) = true; 
end

if size(d2,1) > 2, 
    d2 = d2'; 
end
if size(d2,1) == 1, 
    vectorFlag(2) = 1; 
end

% prepare the bins:
highLim = max( max(d1(1,:)), max(d2(1,:)) ); 
lowLim = min( min(d1(1,:)), min(d2(1,:)) );
binEdges = linspace(lowLim, highLim, numBins + 1);
halfStep = ( binEdges(2) - binEdges(1) ) / 2;
binCenters = binEdges(1:length(binEdges) - 1) + halfStep;

% bin and normalize:
if vectorFlag(1),  % can use hist( )
    hist1 = hist(d1, binCenters);
else   % need a loop:
    hist1 = zeros(size(binCenters));
    for j = 1:length(binEdges) - 1,
        inds = find(d(1,:) > binEdges(j) & d1(1,:) <= binEdges(j + 1) );
        if ~isempty(inds),
            hist1(j) = sum(d1(2, inds) );
        end
    end
end        
hist1 = hist1 / sum(hist1); % normalize

if vectorFlag(2),  % can use hist( )
    hist2 = hist(d2, binCenters);
else   % need a loop:
    hist2 = zeros(size(binCenters));
    for j = 1:length(binEdges) - 1,
        inds = find(d2(1,:) > binEdges(j) & d2(1,:) <= binEdges(j + 1) );
        if ~isempty(inds),
            hist2(j) = sum(d2(2, inds) );
        end
    end
end 
hist2 = hist2 / sum(hist2); % normalize
% hist1 and hist2 are now distributions

% calculate the bhat distance:
bhat = -log( sum (sqrt(hist1.*hist2) ) );

if showPlot,
    figure, 
    plot(binCenters, hist1,'b+'), 
    hold on, 
    plot(binCenters, hist2,'r+')
    title(['Bhat distance = ' num2str(bhat) ', d1 = blue, d2 = red'])
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
    
    
    