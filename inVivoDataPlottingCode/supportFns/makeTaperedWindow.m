% make a tapered window with hamming-style tapers.
% inputs:
% windowCenterSection = width in mSec of central (=1) section. default = 300
% windowSide = width of tapered sections which go on each end. default = 150
% outputs:
% vector of values between 0 and 1.

% Disclaimer: This code has not been refactored or otherwise tidied up (!)

% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

function [hamWin] = makeTaperedWindow(windowCenterSection, windowSide)

if nargin == 0,
    windowCenterSection = 300;
    windowSide = 150;
end

% comment: for some reason, 'hamming( )' prints a figure (the wvtool) when
% called inside a function :(
hamWinSides = hamming(2*windowSide);
hamWinSides = hamWinSides';
hamWin = [ hamWinSides(1:length(hamWinSides)/2), ...
    ones(1, windowCenterSection), hamWinSides(length(hamWinSides)/2+1:end) ];

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