% get median trajectories, ie assume each set of closely-spaced puffs is an
% independent trial, and combine their FRs to get a platonic response to a
% set of puffs

% Disclaimer: This code has not been refactored or otherwise tidied up (!)

% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

function [ trajFRsOdor, trajFRsControl, odorStd] = calcTrajectoryResponseMedian( FRmatrix, odorStimInds )

% inputs:
% 1. FRmatrix = N x M matrix of response FRs for one neuron: each col is response to a set of closely-spaced puffs,
%               each row is the response to the puffs in a particular slot
%               (eg 1st of 5, 2nd of 5, etc). N puffs per set, M sets of puffs.
%               Note that in the input one set of puffs is a col, while in
%               in the output the median trajectory is given as a row vector
% 2. odorStimInds = vector that gives indices of responses corresponding to odors (versus controls) in 
%                   the col vector FRmatrix(:).
%
% outputs:
% 1. trajFRsOdor = 1 x N row vector, the median response to a set of odor puffs
% 2. trajFRsControl = 1 x N row vector, the median response to a set of control puffs.
% 3. odorStd = 1 x N vector of stds for odors only, to give an indication of how noisy the estimate is

% First make a matrix that labels odor vs control:
stimTypeMatrix = zeros(size(FRmatrix(:) ) );
stimTypeMatrix(odorStimInds) = 1;
stimTypeMatrix = reshape(stimTypeMatrix, size(FRmatrix) );
% stimTypeMatrix is the same size as FRmatrix, and has entries = 1 (odor) or 0 (control).

% Not needed: medianStimType(j,:) = median( reshape (stimTypeMatrix, size(stim) ) ); % 1 = odor, 0 = control

% Now calc medians per position in the string of puffs, ie combine all the puff trajectories into one trajectory.
% Note that odor responses are combined separately from control responses:
% (Usually the controls are at the end, ie last col(s), but this cannot be assumed, nor can
%      # of controls. So the method needs to be general.)

% Odor:
sizeOdor = [size(stimTypeMatrix,1), sum(stimTypeMatrix(1,:) ) ]; % how big the odor resp matrix should be
odorResp = FRmatrix(:);
odorResp = odorResp(stimTypeMatrix(:) == 1);
odorResp = reshape(odorResp, sizeOdor);
trajFRsOdor = median( odorResp, 2 ); % col vector
trajFRsOdor = trajFRsOdor'; % row vector

% Control:
sizeControl = [ sizeOdor(1), size(stimTypeMatrix,2) - sizeOdor(2) ]; 
controlResp = FRmatrix(:);
controlResp = controlResp(stimTypeMatrix(:) == 0);
controlResp = reshape(controlResp, sizeControl);
trajFRsControl = median( controlResp, 2 ); 
trajFRsControl = trajFRsControl';  % row vector

% std of odor trajectories:
odorStd = std(odorResp, 1, 2);
odorStd = odorStd'; 

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



