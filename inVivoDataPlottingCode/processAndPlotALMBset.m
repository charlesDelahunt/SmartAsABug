% Code to plot timecourses of a single experiment with odor + octopamine via sugar at proboscis, from Jeff Riffell's lab at U. Washington. Seattle.
% Experimental protocols: 
% “We applied a 5-component mixture as the stimulus, and the duration of each odor pulse was 500 ms, with a 5 sec inter-pulse interval (5 total pulses for each training session).”
% Concurrent octopamine (sugar at proboscis).

% Disclaimer: This code has not been refactored or otherwise tidied up (!)

% Running this script loads 'AL_MB_octo.mat', and plots spontaneous
% firing rates (blue line), spontaneous stds (dotted lines),  responses to odor (red *), and
% responses to mineral oil control (green *).

% The first 10 neurons are AL (5, 9, 10 are duds). The last 12 (11 - 22) are MB (though not sure what MB neurons they are).

% Dependencies: Matlab, Statistics and machine learning toolbox, Signal processing toolbox

% Copyright (c) 2018 Charles B. Delahunt.  delahunt@uw.edu
% MIT License

%-------------------------------------------------------------------------------------------

% Processing steps: 
% get windowed spike rates for spontaneous bursting sections.
% use points at 1 sec intervals between starts and finishes of recording,
% as long as the points are not just after a stim (just before a stim is
% fine)
% compare to control spike rates
% all the various data matrices are indexed by neuron number. Eg, if
% plotting neurons 11 to 22, the first 10 rows (or cols) of a data matrix
% will be zero. neuron 16 data is found in row (or col) 16.


clear all


% specify which neurons to plot
neuronsToPlot = [1:4 6:8];  % Note: duds = 5 9 10   % [1:4 6:8 11:22]; 
% plotting specs, needs to be edited so that plotR*plotC >= max neuron index:
plotR = 4; 
plotC = 6;


plotRawSpont = 0;   % plots spont FRs in each window as blue +'s
plotStdLines = 1; 
calcAndPlotBhatFlag = 0;
plotSpontRunningMedians = 0; % instead of fitting a cubic, just plot the mean of each section of Spont Resps.
closePreviousFigs =  1;
normalizeToInitialFR =  0;  % this sets mean(pre-first-stim windowed FRs) = 1
calcFRMahalDist = 1;
    scaleFlag = 0; % 0 means return straight mahal dists. 1 means scale them to reflect probabilities
    controlNormFlag = 0;  % 1 -> subtract the mahal dist of the control from each odor response.
    spread = 11; % std will be mean of nearest 2*spread + 1 values
plotMedianTrajectories = 0; % plots the platonic trajectory of a set of puffs (ie take median for 1st puff over all sets of puffs, etc)

collectStats = 1; % Note: Requires that plotStdLines = 1 and calcFRMahalDist = 1.
    startOcto = -1;  % ie entire run = octo
    stopOcto = 1e9;
    statsFilename = 'spontAndOdorResponseStats_ALMB';
%%----------------------------------------

% make a hamming window of width windowCenter + 2*windowSide:
windowCenter = 300; % the part == 1
windowSide = 100; % the tapered sides
% make a tapered window of correct size:
window = makeTaperedWindow(windowCenter,windowSide);

delayAL = 300;  % in mS. delay due to 1. odor going down tube; and 2. processing in antennae and RNs
delayMB = delayAL + 30;

%%-------------------------------------------- 
% fit quadratic or cubic:
fitDegree = 3; % can be 2 (fit quadratic) or 3 (fit cubic)

% size of window to calculate std's of spontaneous responses = 2*dist:
dist = 100;  % reducing to 50 adds wiggle but does not appreciably change the main behavior

if closePreviousFigs
    close all
end

% plot each neuron on its own figure to compare delay effects:
% for thisNeuron = [1:4 6:8 11:22];  % ie exclude duds 5 9 10

load AL_MB_octo.mat
% gives the following vars:
% 'MB_AL_ensembles' = 1 x 22 cell, each cell a vertical vector 250 - 6500 elements long;
%  matrix called 'stim' (5 x 4), each col the timestamps of a set of 5 odors, the last col = timestamps of 5 control puffs


%% process prep (there is only one, so no loop over 'ind')

odorStimInds = 1:15;  % ie first 3 cols of 'stim' are odor, last col is control
controlStimInds = 16:20;
ind = 1; % refers to which prep. There is only one prep in this dataset
spikes = MB_AL_ensembles; % ma = MB and AL neurons
% st = stim;

numNeurons = length(neuronsToPlot);
timePoints = [1:850];  % points to sample a window for spont response.

%% process each neuron in this prep:
for j = neuronsToPlot
    % is this neuron in AL or MB:
    AL = j <= 10; % true if from AL
    if AL 
        delay = delayAL;  % to allow travel time for odor
    else
        delay = delayMB;  % allow for 30 mSec travel time AL -> MB
    end
    
    N = spikes{j};     % col vector = the spike times for the j'th neuron
    endTime = 850;
    N = N(N < endTime);  
    stimTimes = stim(:);
    
    % for each point in pointsList, if there is at least 600 mSec before the next stim or finish,
    % get a windowed spike rate. 
 
    % go through this pointsList, getting spont response if it's a good point
    spontResponses(j,:) = getSpontaneousResponses(N,timePoints,stimTimes,window);
    
    % go through stimTimes, getting stim response for each:    
    stimResponses(j,:) = getStimResponses(N, stimTimes, window, delay); % 'delay' is in mSec 
    
    if normalizeToInitialFR
        firstStim = min(stimTimes);
        preStimTimePointInds = find(timePoints < firstStim - length(window)/1000);
        preStimResp = spontResponses(j,preStimTimePointInds);
        normFactor(j) = mean(preStimResp);
        normFactor(j) = ( round(normFactor(j)*100) )/100;
        spontResponses(j,:) = spontResponses(j,:) / normFactor(j);
        stimResponses(j,:) = stimResponses(j,:) / normFactor(j);
    end
    
end % for j = neuronsToPlot 
   
 %% plot the spont responses and stim responses:
 if plotRawSpont
     figure,
     for j = neuronsToPlot
         % index p is used for subplot location:
        p = find(neuronsToPlot == j);
        % special case: if plotting all neurons, move the last 2 MB's to pos 9,10.
        if numNeurons > 16, if j > 20, p = j - 12;else p = j; end, end
        subplot(plotR, plotC, p)

        % get rid of bad start/finish times:
        good = find(spontResponses (j,:) >= 0 );         
        spontResp = spontResponses( j,good );
        spontTimes = timePoints(good);
        % order spontResp by timestamp:
        [spT, orderInds] = sort(spontTimes,'ascend');
        spontResp = spontResp(orderInds);
        spontTimes = spontTimes(orderInds);
        plot( spontTimes, spontResp, '+' );  
        title(num2str(j),'fontweight','bold')
        if normalizeToInitialFR
            title(['#' num2str(j) ', initial FR = ' num2str(normFactor(j))],'fontsize',14,'fontweight','bold'),
        end
        if p == 3 
          title(['MB-AL octo, ' 'delay ' num2str(delay) ' mSec, window ' num2str(length(window)),...
              ' mSec. Blue = spont value, red = response to odor, green = response to control']),
        end % title shows delay and window width
        hold on, grid on,  
        % now plot the stim responses:
        plot(stimTimes, stimResponses(j,:),  'r*'),
	% the next line assumes that the last col only, with 5 rows, is the control (mineral oil):
        plot(stimTimes(end-4:end), stimResponses(j, end-4:end), 'g*')
     end
 end
 
 %% fit a quadratic or cubic (controlled by 'fitDegree') to the spont responses:
 for j = neuronsToPlot
     
    % remove bad start/finish times:
    good = find(spontResponses (j,:) >= 0 ); % ie not = -1        
    spontResp = spontResponses( j,good );
    spontTimes = timePoints(good);
    % order spontResp by timestamp:
    [spT, orderInds] = sort(spontTimes,'ascend');
    spontResp = spontResp(orderInds);
    spontTimes = spontTimes(orderInds);
    
    B = polyfit(spontTimes', spontResp', fitDegree);
    % get the fitted curve:
    x2 = spontTimes;
    if fitDegree == 2
       y2 = B(1)*x2.^2 + B(2)*x2 + B(3);
       stimY2(j,:) = B(1)*stimTimes.^2 + B(2)*stimTimes + B(3);
    end
    if fitDegree == 3
       y2 = B(1)*x2.^3 + B(2)*x2.^2 + B(3)*x2 + B(4); 
       stimY2(j,:) = B(1)*stimTimes.^3 + B(2)*stimTimes.^2 + B(3)*stimTimes + B(4);
    end
    
    thisMeanSubtrSpResp = spontResp - y2; % the actual spont responses, minus the interpolated spont response
   % meanSubtrSpResp{j} = spontResp - y2; % in case needed later
    
    % calc std for different time windows:
    timePointsForStd = 1:10:endTime;

    for t = 1:length(timePointsForStd)
        time = timePointsForStd(t);
        
        % remove the two highest spikes (assumed to be undesirable outliers)
        top = sort(thisMeanSubtrSpResp,'descend');
        rejectLine = top(2);
        
        vals = thisMeanSubtrSpResp(abs(spontTimes - time) < dist & thisMeanSubtrSpResp < rejectLine ) ;
        centeredVals{t,j} = vals;
        std2windows(t,j) = std(vals);
        % note there is some assymetry at the start and end (ie the spontTimes values are
        % mostly from one side of the timepoint)
              
        % find the interpolated spont resp value close to the timePoint:
        if fitDegree == 2
            estSpontResp(t,j) = B(1)*time.^2 + B(2)*time + B(3);
        end
        if fitDegree == 3
            estSpontResp(t,j) = B(1)*time.^3 + B(2)*time.^2 + B(3)*time + B(4);
        end
    end % for t = 1:length(timePointsForStd)

    % calc the est response +/- stds:
    pos1Std(:,j) = estSpontResp(:,j) + std2windows(:,j);
    neg1Std(:,j) = estSpontResp(:,j) - std2windows(:,j);
    pos2Std(:,j) = estSpontResp(:,j) + 2*std2windows(:,j);
    neg2Std(:,j) = estSpontResp(:,j) - 2*std2windows(:,j);

%     % collect the std-normalized spont responses for plotting:
%     % there are lots and lots of duplicates here, but see how it goes
%     % ignore boundary times:
%     stdNormedSpontRespTemp = []; % initialize storage for the std-normed spont responses
%     % this will be used later to check whether the spread of spont resp around the
%     % fitted quadratic is gaussian.        
%     if time > dist && time < endTime - dist,
%         temp = meanSubtrSpResp(abs(spontTimes-time) < dist);
%         temp = temp / std2windows(t,j);
%         stdNormedSpontRespTemp = [stdNormedSpontRespTemp, temp];
%     end     
%     stdNormedSpontResp{j} = stdNormedSpontRespTemp;
    
 end % for j = 1:numNeurons, fitting cubic and std lines loop
 
%% plot this neuron's std2windows: 
if plotStdLines     
    figure,    
    for j = neuronsToPlot
        % index p is used for subplot location:
        p = find(neuronsToPlot == j);
        
        % special case: if plotting all neurons, move the last 2 MB's to pos 9,10.
        if numNeurons > 16, if j > 20, p = j - 12;else p = j; end, end
        
        subplot(plotR, plotC, p),
        xlim([0, endTime - 100])
        % to set yLim, need some gyrations to account for Nan results:
        validStimResp = stimResponses(j,:);
        validStimResp = validStimResp (~isnan(validStimResp) ) ;
        if ~isempty(validStimResp)
            upperLim = max(validStimResp);
        else 
            upperLim = 1;
        end
        ylim([-0.2, upperLim + 1 ]), % finally ylim defined
        title(num2str(j), 'fontweight','bold')
        if normalizeToInitialFR
            title(['#' num2str(j) ', initial FR = ' num2str(normFactor(j))],'fontsize',14,'fontweight','bold'),
        end
        hold on,

      % plot the estimated spontaneous response +/- 2 stds, which vary with time:
      % because timePoints is now stimTimes, which is not ordered (control
      % times are given last), reorder for plotting:
        [~,ord] = sort(timePointsForStd, 'ascend');
        plot (timePointsForStd(ord), estSpontResp(ord,j))
        plot(timePointsForStd(ord), pos1Std(ord, j), ':');
        plot(timePointsForStd(ord), pos2Std(ord, j), ':');
        plot(timePointsForStd(ord), neg1Std(ord, j), ':');
        plot(timePointsForStd(ord), neg2Std(ord, j), ':');

        % now plot the stim responses:
        plot(stimTimes, stimResponses(j,:),  'r+'),
        plot(stimTimes(end-4:end), stimResponses(j, end-4:end), 'g+')     
    end % for j = neuronsToPlot   % put the figure title in the lower right
             subplot(plotR,plotC,plotR*plotC - 1),
    title(['MB-AL octo, ' 'delay ' num2str(delay) ' mSec, window ' num2str(length(window)),...
        ' mSec. Blue = mean spont (also 1,2 std dotted lines); red = odor response; green = control response.']),
end

  %% normalize FRs by subtracting the estSpontResp and dividing by the std:

    if calcFRMahalDist
        
        figure,
        
        clear medianStimTimes medianMahalResponses medianStimType medianTrajectoryMahalRespOdor medianTrajectoryMahalRespControl
        mahalStimResponses = zeros(size (stimResponses) ); % initialize
        
        medianStimTimes(j,:) = median(stim);  % used to plot spont responses
        
        upperMahalRespLim = 19.08;
        for j = neuronsToPlot
            
           mahalStimResponses(j,:) =...
               calcMahalProbFRs(stimResponses(j,:), stimTimes, (estSpontResp(:,j))', (std2windows(:,j))', timePointsForStd, spread, upperMahalRespLim);
           
            % calc medians in 2 directions:
            % the most relevant is the second, ie finding a platonic form
            % of response trajectory to a set of fast puffs.
            
            % 1a. make a vector that labels odor vs control:
            stimTypeMatrix = zeros(size(stim(:) ) );  % temp orarily a col vector, not a matrix
            stimTypeMatrix(odorStimInds) = 1;
            stimTypeMatrix = reshape(stimTypeMatrix, size(stim) ); % now it's a matrix
            % now stimTypeMatrix has entries = 1 (odor) or 0 (control).
            medianStimType(j,:) = median( reshape (stimTypeMatrix, size(stim) ) ); % 1 = odor, 0 = control
            
            % 1b. calc medians per set of puffs. Ie, combine a set of 5 or so closely-spaced puffs
            % into one response:
            medianStimTimes(j,:) = median(reshape(stimTimes,size(stim) ) );
            medianMahalResponses(j,:) = median( reshape( (mahalStimResponses(j,:))', size(stim) ) );
            
            % 2. calculate the platonic ideal of a response trajectory given a set of puffs:
            FRmatrix = reshape( (mahalStimResponses(j,:))', size(stim) );
            [ trajFRsOdor, trajFRsControl, odorStd ] = calcTrajectoryResponseMedian( FRmatrix, odorStimInds );
            medianTrajectoryMahalRespOdor (j,: ) = trajFRsOdor;
            medianTrajectoryMahalRespControl(j,:) = trajFRsControl; 

            % index p is used for subplot location:
            p = find(neuronsToPlot == j);

            subplot(plotR, plotC, p),

            title(num2str(j), 'fontweight','bold')
            if normalizeToInitialFR
                title(['#' num2str(j) ',initial FR = ' num2str(normFactor(j)) ] ) %,...
                    % 'fontsize',14,'fontweight','bold'),
            end
            hold on,
            if plotMedianTrajectories % case: we want to see the platonic response to a set of fast puffs
                plot(medianTrajectoryMahalRespOdor(j,:),'r*')
                plot(medianTrajectoryMahalRespControl(j,:),'g*')
                plot(medianTrajectoryMahalRespOdor(j,:) - odorStd,'b+')
                plot(medianTrajectoryMahalRespOdor(j,:) + odorStd,'b+')
                
                lowY = floor( min( [ medianTrajectoryMahalRespControl(:); medianTrajectoryMahalRespOdor(:) ] ) - max(odorStd) );
                highY = ceil( max( [ medianTrajectoryMahalRespControl(:); medianTrajectoryMahalRespOdor(:) ] ) + max(odorStd) );
                ylim( [ lowY, highY ] )
                xlim ( [0.5, length(odorStd) + 0.5])
                
            else   % case: we want to see all the responses, visually grouped into sets of fast puffs:
                xlim([0, max(timePoints) - 100])
                % use same y scale for all neurons in a prep:
                lowY = floor(min (mahalStimResponses(j,:)) );
                highY = ceil(max (mahalStimResponses(j,:)) ); 
                ylim([lowY, highY] )
                % plot stim responses:
                plot(stimTimes(odorStimInds), mahalStimResponses(j,odorStimInds),  'r*'),
                % plot control responses:
                plot(stimTimes(controlStimInds), mahalStimResponses(j, controlStimInds), 'g*') 
                % plot median responses:
                plot(medianStimTimes(j,:), medianMahalResponses(j,:),'b*')
            end
            
%             plot( [ timePointsForStd(closestInds(j,2) - spread), timePointsForStd (closestInds(j,2) + spread ) ], [1 1],'k+')
%             plot( timePointsForStd(closestInds(j,2) - spread), 1,'k+'),
%             plot( timePointsForStd (closestInds(j,2) + spread ), 1,'k+')
            grid on

        end % for j = neuronsToPlot

            % put the figure title in the lower right
             subplot(plotR,plotC,plotR*plotC - 1),
            
         title(['MB-AL octo, prep ' num2str(ind) , '. Mahal distance. blue = median spont, red = odor resp, green = control resp'])
         % title(['MAHAL. prep ' num2str(ind) '. delay ' num2str(delay) ' mSec, window ' num2str(length(window)), ' mSec'])
    end   

 %% calculate Bhattacharyya distance from a gaussian made with the
    % estimated std for this window (assume mean = 0 since we subtracted the estResp).
    % calc the hypothetical gaussian:
if calcAndPlotBhatFlag
   
    for j = neuronsToPlot
        for t = 1:size(std2windows,1) % should be the same as length(timePointsForStd)
            s = std2windows(t,j);
            theseVals = centeredVals{t,j};
            bins = linspace(-4*s, 4*s, 201);
            halfStep = ( bins(2) - bins(1) ) / 2;
            binCenters = bins(1:length(bins) - 1) + halfStep;
            gaussian = exp( -0.5* (binCenters.^2)/s^2 ); % make the equiv gaussian
            % set the gaussian = 0 for any negative bins, ie bins corresponding to < 0 spontFR:
            gaussian(binCenters < min(vals)) = 0;
            formattedGaussian = [binCenters; gaussian]; % 'bhatt...' needs this gaussian as a 2 x N matrix
            if isempty(theseVals)
                Bdist(t,j) = -1;
            else
                Bdist(t,j) = bhattacharyya(theseVals, formattedGaussian);
            end
        end % for t
    end

    % for comparison: for N = 180 (= approx number of spontResps in each window), 
    % a randomly generated gaussian is 0.15 +/- 0.012 Bdist from its theoretical
    % version (given the estimated mu, sigma)
    % for N = 100, a randomly generated gaussian is 0.28 +/- 0.3 Bdist from its theoretical
    % version (given the estimated mu, sigma)
    % these numbers found using bhattacharyyaDistForRandomGaussians.m
    
    % plot the bhat distances
     figure,
     for j = neuronsToPlot
         % index p is for subplot location:
        p = find(neuronsToPlot == j);
        % special case: if plotting all neurons, move the last 2 MB's to pos 9,10.
        if numNeurons > 16, if j > 20, p = j - 12;else p = j; end, end
        
        subplot(plotR, plotC, p),
        this = Bdist(:,j);
        plot(timePointsForStd(this >= 0), this(this >= 0) ),
        title(num2str(j))
        grid on   
     end
     mtit(['Bhattacharyya distance: delay ' num2str(delay) ' mSec, window ' num2str(length(window)), ' mSec']),     
end

%% calculate and plot the running means of the spontaneous responses.

if plotSpontRunningMedians
    % since the std calc was for spontResp - cubicEstimate, re-do the windowing operation:
    for j = neuronsToPlot
        for t = 1:length(timePointsForStd)
             % remove bad start/finish times:
            good = find(spontResponses (j,:) >= 0 ); % ie not = -1        
            spontResp = spontResponses( j,good );
            spontTimes = timePoints(good);
            % order spontResp by timestamp:
            [spT, orderInds] = sort(spontTimes,'ascend');
            spontResp = spontResp(orderInds);
            spontTimes = spontTimes(orderInds);
            
            time = timePointsForStd(t);
            vals = spontResp(abs(spontTimes - time) < dist ) ;            
            
            if isempty(vals)
                runningMedian(t,j) = -1;
            else
                runningMedian(t,j) = median(vals);
            end
            if isempty(vals)
                runningMean(t,j) = -1;
            else
                runningMean(t,j) = mean(vals);
            end
        end % for t
    end
    
    % plot the running medians:
     figure,
     for j = neuronsToPlot
         % index p is for subplot location:
        p = find(neuronsToPlot == j);
        % special case: if plotting all neurons, move the last 2 MB's to pos 9,10.
        if numNeurons > 16, if j > 20, p = j - 12;else p = j; end, end
        
        subplot(plotR, plotC, p),
        plot(timePointsForStd, runningMedian(:,j),'linewidth',2 ),
        hold on,
        plot(timePointsForStd, runningMean(:,j),'k','linewidth',2 ),
        % plot the stim times:
        plot(stimTimes, stimResponses(j,:),  'r+'),
        plot(stimTimes(end-4:end), stimResponses(j, end-4:end), 'g+')
        ylim([0,2])
        title([ num2str(j) ', initial FR = ' num2str(normFactor(j)) ] )
        grid on   
     end
     mtit(['MB-AL octo,' '  running median (blue), mean (black): delay ' num2str(delay) ' mSec, window ' num2str(length(window)), ' mSec']),     
end
    %% collect stats
    
    if collectStats
%   goal: for each neuron, calculate the following:
%         CAUTION: if octo input is mixed, we need to add code to segregate
%         on the basis of octo vs no octo. The tricky part is medS, since
%         this requires the spont response to be partitioned.
%         1. median Spont
%         2. std Spont
%         3. as a check: coeff of variation of spont ie std(median spont values)/median Spont
%         (let R = median of odor responses for each cluster, so R is a vector with an entry for each cluster)
%         4. median R 
%         as checks:
%         5. max R
%         6. min R
%         7. coeff of var of R
%         for preps with both octo and no octo:
%         8. medSpont(octo) - medSpont(noOcto) / std Spont(noOcto)
        % make list of odor stim indices that is <= # cols of stim matrix
        %       (since each col is one type (odor or control) and since
        %       medianMahalResponses condenses each col into one value:
        % at the same time, separate these into octo and noOcto lists of inds:
        % the code below assumes that octo is added in the middle of a run, once, window defined by startOcto and stopOcto.
        noOctoOdorStimInds = find(stimTimes < startOcto | stimTimes >= stopOcto);
        noOctoOdorRespInds = unique(ceil(noOctoOdorStimInds/size(stim,1)));
        octoOdorStimInds = find(stimTimes >= startOcto & stimTimes < stopOcto);
        octoOdorRespInds = unique(ceil(octoOdorStimInds/size(stim,1)));
        
        % separate spont inds into octo and noOcto:
        noOctoSpontInds = find(timePointsForStd < startOcto | timePointsForStd >= stopOcto);
        octoSpontInds = find(timePointsForStd >= startOcto & timePointsForStd < stopOcto);
        noOctoStimInds = find(stimTimes < startOcto | stimTimes >= stopOcto);
        octoStimInds = find(stimTimes < startOcto | stimTimes >= stopOcto);
        
        for j = neuronsToPlot
            % noOcto stats:
            if ~isempty(noOctoSpontInds)
                medS(ind,j) = median(estSpontResp(noOctoSpontInds,j));
                sigmaS(ind,j) = median(std2windows(noOctoSpontInds,j));
                coeffVarS(ind,j) = abs( std(estSpontResp(noOctoSpontInds,j))/medS(ind,j) );
                R(ind,j) = median(medianMahalResponses(j,noOctoOdorRespInds));
                maxR(ind,j) = max(medianMahalResponses(j,noOctoOdorRespInds));
                minR(ind,j) = min(medianMahalResponses(j,noOctoOdorRespInds));
                coeffVarR(ind,j) = abs( std(medianMahalResponses(j,noOctoOdorRespInds)) / R(ind,j) );
            end
            % octo stats:
            if ~isempty(octoSpontInds)
                octoMedS(ind,j) = median(estSpontResp(octoSpontInds,j));
                octoSigmaS(ind,j) = median(std2windows(octoSpontInds,j));
                octoCoeffVarS(ind,j) = abs( std(estSpontResp(octoSpontInds,j))/octoMedS(ind,j) );
                octoR(ind,j) = median(medianMahalResponses(j,octoOdorRespInds));
                octoMaxR(ind,j) = max(medianMahalResponses(j,octoOdorRespInds));
                octoMinR(ind,j) = min(medianMahalResponses(j,octoOdorRespInds));
                octoCoeffVarR(ind,j) = abs( std(medianMahalResponses(j,octoOdorRespInds)) / octoR(ind,j) );
            end
            
            bothOctoAndNoOctoExist = length(noOctoSpontInds) > 1 && length(octoSpontInds) > 1; % ie there are both noOcto and octo regions.
            if bothOctoAndNoOctoExist % ie there are both octo and no octo responses
                deltaSDueToOcto(ind,j) = ( octoMedSo(ind,j) - medS(ind,j) ) / sigmaS(ind,j);
            end
        end
    end     
 
%% save stats if collected    
if collectStats
    % now save all this in a struct:
    if exist('medS'), s.medS = medS; else s.medS = []; end
    if exist('sigmaS'), s.sigmaS = sigmaS; else s.sigmaS = []; end
    if exist('coeffVarS'), s.coeffVarS = coeffVarS; else s.coeffVarS = []; end
    if exist('R'), s.R = R; else s.R = []; end
    if exist('maxR'), s.maxR = maxR; else s.maxR = []; end
    if exist('minR'), s.minR = minR; else s.minR = []; end
    if exist('coeffVarR'), s.coeffVarR = coeffVarR; else s.coeffVarR = []; end

    if exist('octoMedS'), s.octoMedS = octoMedS; else s.octoMedS = []; end
    if exist('octoSigmaS'), s.octoSigmaS = octoSigmaS; else s.octoSigmaS = []; end
    if exist('octoCoeffVarS'), s.octoCoeffVarS = octoCoeffVarS; else s.octoCoeffVarS = []; end
    if exist('octoR') , s.octoR = octoR; else s.octoR = []; end
    if exist('octoMaxR') , s.octoMaxR = octoMaxR; else s.octoMaxR = []; end
    if exist('octoMinR') , s.octoMinR = octoMinR; else s.octoMinR = []; end
    if exist('octoCoeffVarR'), s.octoCoeffVarR = octoCoeffVarR; else s.octoCoeffVarR = []; end

    if exist('deltaSDueToOcto'), s.deltaSDueToOcto = deltaSDueToOcto; else s.deltaSDueToOcto = []; end
    
    spontAndOdorResponseStats = s;
    save ( statsFilename, 'spontAndOdorResponseStats' )
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
