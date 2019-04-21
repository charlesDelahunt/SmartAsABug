% calculate probability-weighted Mahalanobis distance firing rates.
% Given FRs, stds of spont resp, and estimated spont resp:
% Find mahal distance of each stim response. Then use exp(mahal )  = exp(m) similar to but != exp(0.5*m^2) = 1/gaussian prob.
% exp(m) approx= 1/gaussian up to m = 2, but it has a much less steep exponential increase as m -> 3: exp(3) = 19, exp(0.5*3^2) = 90
% Also apply an upper max to avoid extremely high FR distances (19 for exp(mahal)
% apply to data from one neuron with N stims

% Disclaimer: This code has not been refactored or otherwise tidied up (!)

% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

function [mahalProbFRs] = calcMahalProbFRs( stimResp, stimTimes, estSpontResp,...
    spontRespStd, timePointsForStd, spread, upperMahalRespLim, scaleFlag, controlNormFlag, controlInds)

% inputs:
% 1. stimResp = raw stim response FRs: vector 1 x N
% 2. stimTimes = 1 x N vector
% 3. estSpontResp = estimated spont response at timePointsForStd (below): vector 1 x M
% 4. spontRespStd = std of spont responses using a moving window: vector 1 x M
% 5. timePointsForStd = centers of moving windows used to find stds: vector 1 x M
% 6. spread = half of window width over which to average spont resp std's for
%             a given stim. scalar. This is in indices, not seconds, so it
%             should probably be changed at some point :)
% 7. upperMahalRespLim = upper limit for prob-weighted distances . scalar. optional argin
% 8. scaleFlag = indicate whether to scale mahal dists to better represent probability. optional argin.
% 9. controlNormFlag = subtract the median mahal dist of the control (mineral oil) from each stim resp. optional argin.
% 10. controlInds  = indices of the controls in stimResp and stimTimes
%
% outputs:
% 1. mahalProbFRs = 1 x N vector = probability-weighted Mahal distance values for FRs
% 
% method: 
% 1. a = (rawFRs - estSpontResp) / (mean(spontRespStd) % Gives Mahal dist
% 2. b = exp(a) - 1 % b(0) = 0; b(1) = 1.7; b(2) = 6.4; b(3) = 19.1
% 3. output = min(b, upperMahalRespLim).

if nargin < 10, % since this step requires 2 argins
    controlNormFlag = 0;
end
if nargin < 8,
    scaleFlag = 0; % default: return straight mahal distances
end
if nargin < 7,
    upperMahalRespLim = 19.1; % so 3 std's is max 'distance'
end
if nargin < 6,
    spread = 3; % 11;  % not necessarily a good catch-all value.
end

for k = 1:length(stimResp),
    closest = find( abs( timePointsForStd - stimTimes(k) ) == min( abs( timePointsForStd - stimTimes(k) ) ) ); % ind of closest window
    range = [ max( closest - spread, 1): min( closest + spread, length(timePointsForStd) ) ]; % use the stds from these windows
    mahalStimResponses(k) = ( stimResp(k) - median(estSpontResp(range) ) ) / median( spontRespStd(range) ) ;
end

% if wished, normalize by subtracting the median mahal of the control:
if controlNormFlag,
    normFactor = median(mahalStimResponses(controlInds));
    mahalStimResponses = mahalStimResponses - normFactor;
end
% if wished, scale the std distances by exp to weight the dist by probabilities
if scaleFlag,
    mahalProbFRs = sign(mahalStimResponses).*( exp(abs(mahalStimResponses) ) - 1 );
else
    mahalProbFRs = mahalStimResponses;
end

% also, apply a max limit, eg to make anything > 3 std about the same:
mahalProbFRs = min(mahalProbFRs, upperMahalRespLim*ones(size(mahalProbFRs)) );
mahalProbFRs = max(mahalProbFRs, -upperMahalRespLim*ones(size(mahalProbFRs)) );

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
