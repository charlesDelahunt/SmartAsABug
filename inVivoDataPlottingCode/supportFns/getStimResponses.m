% get windowed FR responses to stimuli for a single neuron
% inputs:
%   spikes: vector of spike timestamps
%   stimTimes: vector of times when stimuli were given. This is NOT in time order
%   window: n-length window (n = # mSec), with central values = 1 and ramps at each end
%   delay: scalar = # mSec to wait before applying window
% outputs:
%   stimResponses = the spike counts in the window. 
 
% Disclaimer: This code has not been refactored or otherwise tidied up (!)

% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

function [stimResponses] = getStimResponses(spikes, stimTimes, window, delay)

winSize = length(window) / 1000; % length of window in sec
winArea = sum(window)/1000; % this is the weighted length of the window. We divide by this to get the FR in 
			    % spikes/sec rather than spikes/window.

% for each timePoint, see if it is far from stimTimes and also in a valid
% region in terms of starts and stops. If both these hold, get a spike count:
 for k = 1:length(stimTimes),
        point = stimTimes(k);
         % shift the spike time stamps for this point:
        shiftedSpikes = spikes - (point + delay/1000);  % convert 'delay' to sec
        theseSpikes = shiftedSpikes(shiftedSpikes > 0 & shiftedSpikes < winSize); % the timestamps captured by the window
                
        % remove any shifted timestamps that round to 0 mSecs, ie are super-close to the start of the window.
        % (I don't remember why this matters):
        theseSpikes = theseSpikes(theseSpikes*1000 > 0.5); % might be empty, but it still works
        % each spike gets weighted by the hamming window value at its timestamp 
        wtTheseSpikes = window(round(theseSpikes*1000) ); % *1000 turns timestamp (sec) into index(1:#mSec)
        stimResponses(k) = sum(wtTheseSpikes) / winArea ;        
 end  % for k 

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