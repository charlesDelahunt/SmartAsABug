% get spontaneous responses for a single neuron, given
% inputs:
%   spikes: vector of spike timestamps
%   timePoints: vector of points to check. values = time in sec
%   stimTimes: vector of times where stimuli were given. these need to be
%       avoided by the spont response calculator. They are assumed to be 200 mSec stims.
%   window: n-length window (n = # mSec)
%   The next 2 argin are optional: 
%   starts: when the sensors are turned on and off, the time between
%         starts and stops is no good
%   stops: when the sensors are turned off
%   excludeNumWindows: how many window lengths after a stim should not be used for calculating spont responses
% outputs:
%   spontResp = the spike counts in the window. If the timepoint was invalid,
%       spontResp = -1. there is a value in spontResp for each timePoint.
%       So size(spontResp) = size(timePoints) 

% Disclaimer: This code has not been refactored or otherwise tidied up (!)

% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

function [spontResponse] = getSpontaneousResponses(spikes, timePoints, stimTimes, window, postStimExclusion)

winSize = length(window)/1000; % length of window in sec
winArea = 1000/sum(window); % = how many weighted window lengths fit into 1 second, ie windows/sec. 
                            % We multiply the windowed spike count by this to get the FR in spikes/sec rather than spikes/window.
N = spikes; % rename for ease

if isempty(N),    % case: this neuron has no spikes at all over all timePoints, so spontResponse is easy.
    spontResponse = zeros(size(timePoints));
else
    for k = 1:length(timePoints),
        point = timePoints(k);
        % case: there are some spikes. Process all the timepoints.
        % for each timePoint, assume it is at the very start of a window. Get a spike count. 
        % shift the spike time stamps for this point:
        shiftedN = N - point;
        thisN = shiftedN(shiftedN > 0 & shiftedN < winSize); % the timestamps captured by the window
        if isempty(thisN),
            spontResponse(k) = 0;
        else
            % remove anything that rounds to 0:
            thisN = thisN(thisN*1000 > 0.5); % *1000 turns timestamp (sec) into index(1:#mSec)
            % each spike gets weighted by the hamming window value at its timestamp:
            theseSpikes = zeros(size(window));
            theseSpikes(round(thisN*1000)) = 1;  % 1's and 0's, 1's where-ever there was a spike.
            wtThisN = window.*theseSpikes;  % weight the spikes, so those on ramps count less
            % wtThisN = window(round(thisN*1000) ); 
            spontResponse(k) = winArea*sum (wtThisN );
        end
    end % for k
end % if isempty(N)

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