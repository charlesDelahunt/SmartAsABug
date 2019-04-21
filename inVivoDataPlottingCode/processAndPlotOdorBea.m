% Code to plot timecourses of experiments with various odors only (no octopamine or sugar reward), from Jeff Riffell's lab at U. Washington. Seattle.

% Disclaimer: This code has not been refactored or otherwise tidied up (!)
 
% Running this script loads 'BEA_data.mat', and plots spontaneous
% firing rates (blue line), spontaneous stds (dotted lines),  responses to odor (red, blue, black, or magenta +'s), and
% responses to mineral oil control (green +'s). 
% Odor-to-color map: P9 = red, Pn = blue, Mn = magenta, misc = black, control = green

% data file 'BEA_data.mat' = 1 x 12 struct array with fields: 
%     prepName
%     numNeurons: scalar
%     stimTimestamps: <5 x 23 double> Last col = control puffs (not chronological)
%     spikes: 1 x numNeurons struct. each field = 'timestamps', 1 x N double = spike times for that neuron, values between 0 and 6500, N ~ 2000.
%     starts: 1 x num double
%     finishes: 1 x num double
%     odorColumns:  0 = control (min oil), Odors: 1 = misc, 2 = P9, 3 = Pn (behavior) 4 = Mn (no behavior)     %% These are the different  odor types used
%     octoColumns: all 0 = no octo
 
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

% specify which preparations to plot
preparationsToPlot = [2:12];     % there are 12 total. #1 is mostly dead neurons

plotRawSpont =  0;
plotStdLines =  1;
calcAndPlotBhatFlag =  0; % 1;
plotSpontRunningMedians = 0; % 1; % instead of fitting a cubic, just plot the median of each section of Spont Resps.
closePreviousFigs = 0;
normalizeToInitialFR =   0; % 1;% this sets mean(pre-first-stim windowed FRs) = 1
calcFRMahalDist = 0;
plotMedianTrajectories = 0; % ie combine all groups of closely-spaced puffs into one trajectory, 
%                           by combining all the 1st puffs in the groups,
%                           all the 2nd puffs in the groups, etc
%                           This flag prevents the total mahal FR responses from being plotted

% fit quadratic or cubic:
fitDegree = 0; % can be 0 (spline fit to median), 2 (fit quadratic) or 3 (fit cubic)
% on the PN dataset, cubics are poorly conditioned but the curves still look good.

%-----------------------------------------------------------------

% make a hamming window of width windowCenter + 2*windowSide:
windowCenter = 300; % the part == 1
windowSide = 100; % the tapered sides
% make a tapered window of correct size:
window = makeTaperedWindow(windowCenter,windowSide);

delayAL = 300;  
delayMB = delayAL + 30; % in mS. delay due to 1. odor going down tube; and 2. processing in antennae and RNs

%-------------------------------------------- 
    
% the size of window to calculate std's of spontaneous responses will be 2*dist:
dist =  100; % 50 is too noisy for this dataset due to regions without recordings
stdTimePointSpacing = 50; % secs between estimates of std and median calcs

if closePreviousFigs,
    close all
end

% plot each neuron on its own figure to compare delay effects:
% for thisNeuron = [1:4 6:8 11:22];  % ie exclude duds 5 9 10

load BEA_data.mat
% data = 1 x 12 struct array with fields: 
%     prepName
%     numNeurons: scalar
%     stimTimestamps: <5 x 23 double> Last col = control puffs (not chronological)
%     spikes: 1 x numNeurons struct. each field = 'timestamps', 1 x N double =
%           spike times for that neuron, values between 0 and 6500, N ~ 2000.
%     starts: 1 x num double
%     finishes: 1 x num double
%     odorColumns:  0 = control (min oil), Odors: 1 = misc, 2 = P9, 3 = Pn
%     (behavior) 4 = Mn (no behavior)
%     octoColumns: all 0 = no octo

%% start looping through preps

for ind = preparationsToPlot,
    % clear old variables
    clear spontResponses stimResponses normFactor centeredVals std2windows estSpontResp 
    clear rawSpontResponses rawStimResponses
    clear pos1Std pos2Std neg1Std neg2Std Bdist 
    clear runningMean runningMedian pos1StdMed pos2StdMed neg1StdMed neg2StdMed
    clear stimTimes stimY2 
    clear odorStimInds controlStimInds spontStimInds pheromoneStimInds airStimInds
    clear preOctoInds octoInds postOctoInds
    
    % assign the variables based on this preparation:
    prepName = data(ind).prepName;
    stim = data(ind).stimTimestamps;
    sp = data(ind).spikes;
    
    neuronsToPlot = 1:data(ind).numNeurons; % default, process all the neurons.
    
    numNeurons = length(neuronsToPlot);
    % we want the neuron spike times in cols, each col one neuron:
    for i = 1:length(sp),
        spikes{i} = ( sp(i).timestamps )';
        lastStamp(i) = max(spikes{i});
    end
    timePoints = [ 1: floor( max(lastStamp) ) ];  % points at which to sample a window for spont response.
    
    starts = data(ind).starts; % 'starts' and 'stops' are assumed to have the same length
    stops = data(ind).finishes;
    
    % now remove timePoints that are in dead zones:
    temp = [];
    for i = 1:length(starts),
        temp = [temp, timePoints(timePoints > starts(i) & timePoints < stops(i) ) ];
    end
    timePoints = temp;    
  
    %%
    
    % plotting specs, may need to be edited:
    if data(ind).numNeurons == 1 % pheromone
        plotR = 1;
        plotC = 2;
    elseif data(ind).numNeurons <=11,
        plotR = 3; 
        plotC = 4;
    else                 % assume numNeurons < 16 always
        plotR = 4;
        plotC = 4;
    end

    %% process each neuron:
    
    for j = neuronsToPlot,
        % is this neuron in AL or MB:
         AL = true; % all this dataset are PNs in AL
        if AL, 
            delay = delayAL;  % to allow travel time for odor
        else
            delay = delayMB;  % not used here since all PNs. allows for 30 mSec travel time AL -> MB
        end

        N = spikes{j};     % col vector = the spike times for the j'th neuron

%         endTime = max(timePoints); % these 2 lines are not useful anymore
%         N = N(N < endTime);

        stimTimes = stim(:);
        % track which stims are odor, which are control:
        numPuffsPerTrial = size(stim,1);
        odorColumns = data(ind).odorColumns; % 1 means that col of stim = odor; 0 means that col is control
        odorFlag = [];
        for a = 1:length(odorColumns),
            odorFlag = [odorFlag, odorColumns(a)*ones(1,numPuffsPerTrial)]; % adds some 1s or zeros
        end
        octoColumns = data(ind).octoColumns; % 1 means that col of stim = odor; 0 means that col is control
        octoFlag = [];
        for a = 1:length(octoColumns),
            octoFlag = [octoFlag octoColumns(a)*ones(1,numPuffsPerTrial)]; % adds some 1s or zeros
        end
        MnStimInds = find(odorFlag == 4);
        PnStimInds = find(odorFlag == 3); 
        miscStimInds = find(odorFlag == 1);
        controlStimInds = find(odorFlag == 0);
        P9StimInds = find(odorFlag == 2);
        
        preOctoInds = find(octoFlag == -1); % irrelevant
        octoInds = find(octoFlag == 0);
        postOctoInds = find(octoFlag == 1);        % irrelevant

        % for each point in pointsList, if there is at least 600 mSec before the next stim or finish,
        % get a windowed spike rate.

        % go through this pointsList, getting spont response if it's a good point
        % note that one prep has spont responses pre-octo. These are
        % spontaneous responses, so they should not be part of the
        % 'stimTimes' sent into this function. If an issue, easy to remove them.
        rawSpontResponses(j,:) = getSpontaneousResponses(N, timePoints, stimTimes, window); % , starts, stops);

        % go through stimTimes, getting stim response for each:    
        rawStimResponses(j,:) = getStimResponses(N, stimTimes, window, delay); % 'delay' is in mSec 
        
        
        %% Do calcs for initial FR. 
        % Instead of using spike counts from
        % before the FIRST stim, use spike counts before the first ODOR
        % stim. This is because the first control stim is often too close
        % to the start and gives distorted initial FRs:
        
        firstStim = min(stimTimes(miscStimInds));
       
        preStimTimePointInds = find(timePoints < firstStim - length(window)/1000);
        preStimResp = rawSpontResponses(j,preStimTimePointInds);
        preStimResp = preStimResp(preStimResp >= 0); % need to avoid the 'bad' timepoints, which have stimResp = -10k
        normFactor(j) = mean(preStimResp);
        normFactor(j) = ( round(normFactor(j)*100) )/100; 
        normFactor = max(normFactor, 0.1); % in a few cases normFactor(j) = 0, resulting in NaN data for subplot
        
        % apply normalization if indicated:
        if normalizeToInitialFR,
            spontResponses(j,:) = rawSpontResponses(j,:) / normFactor(j);
            stimResponses(j,:) = rawStimResponses(j,:) / normFactor(j);
        else 
            spontResponses(j,:) = rawSpontResponses(j,:);
            stimResponses(j,:) = rawStimResponses(j,:);
        end
        %-----------------------
        
        %% medians:
        
        % calculate and plot the running medians of the spontaneous responses.
        % first list the time points that will be used to calc std:
        timePointsForStd = 1:stdTimePointSpacing:max(timePoints);
        % now remove timePoints that are in dead zones:
        temp = [];
        for i = 1:length(starts),
            temp = [temp, timePointsForStd(timePointsForStd > starts(i) & timePointsForStd < stops(i) ) ];
        end
        timePointsForStd = temp; 
        
        % windowing:
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

            if isempty(vals),
                runningMedian(t,j) = -1;
            else
                runningMedian(t,j) = median(vals);
            end
            if isempty(vals),
                runningMean(t,j) = -1;
            else
                runningMean(t,j) = mean(vals);
            end
        end % for t
        
    end % for j = neuronsToPlot 

     %% plot the spont responses and stim responses:
     
     if plotRawSpont, 
         figure,
         for j = neuronsToPlot,
             % index p is used for subplot location:
            p = find(neuronsToPlot == j);
            
            subplot(plotR, plotC, p)

            % get rid of bad start/finish times:
            good = find(rawSpontResponses (j,:) >= 0 );         
            rawSpontResp = rawSpontResponses( j,good );
            spontTimes = timePoints(good);
            % order spontResp by timestamp:
            [spT, orderInds] = sort(spontTimes,'ascend');
            rawSpontResp = rawSpontResp(orderInds);
            spontTimes = spontTimes(orderInds);
            plot( spontTimes, rawSpontResp, '+' );  
            title(num2str(j),'fontsize',14,'fontweight','bold')
            if normalizeToInitialFR,
                title(['#' num2str(j) ', initial FR = ' num2str(normFactor(j))] ) % ,...
                    % 'fontsize',14,'fontweight','bold'),
            end
           
            hold on, grid on,  
            % now plot the stim responses:
            % do them one at a time to get all the color/shape variations:
            for k = 1:length(stimTimes),
                shape = '+'; 
                col = 'k'; % s
                if ismember(k,controlStimInds), col = 'g'; end
                if ismember(k,PnStimInds), col = 'b'; end
                if ismember(k,P9StimInds), col = 'r'; end
                if ismember(k,MnStimInds), col = 'm'; end
                plot(stimTimes(k), rawStimResponses(j,k),  [col shape] ),
                xlim([0, max(timePoints) ])
            end
         end
         % put the figure title in the lower right
         subplot(plotR,plotC,plotR*plotC),
            title(['BEA. prep ' num2str(ind) ]) % ', delay ' num2str(delay) ' mSec, window ' num2str(length(window)), ' mSec']),
            
     end % if plotRawSpont

     %% fit a quadratic or cubic (controlled by 'fitDegree') to the spont responses:
     
     for j = neuronsToPlot,

        % remove bad start/finish times:
        good = find(spontResponses (j,:) >= 0 ); % ie not = -1        
        spontResp = spontResponses( j,good );
        spontTimes = timePoints(good);
        % order spontResp by timestamp:
        [spT, orderInds] = sort(spontTimes,'ascend');
        spontResp = spontResp(orderInds);
        spontTimes = spontTimes(orderInds);
        
        % the cubic throws ill-conditioned warnings for this dataset. Turn
        % off this warning:
        id = 'MATLAB:polyfit:RepeatedPointsOrRescale';
        warning('off',id);
        
        if fitDegree == 0,
            y2 = spline(timePointsForStd, runningMean(:,j), spontTimes); % note this uses mean, not median. 
                                 % There's not a big difference in stds either way
            stimY2(j,:) = spline(timePointsForStd, runningMean(:,j),stimTimes);
        else   % case: polynomial fit
            B = polyfit(spontTimes', spontResp', fitDegree);
        end
        
        warning('on',id);
        
        % get the fitted curve:
        x2 = spontTimes;
        if fitDegree == 2,
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
        timePointsForStd = 1:stdTimePointSpacing:max(timePoints);
        % now remove timePoints that are in dead zones:
        temp = [];
        for i = 1:length(starts),
            temp = [temp, timePointsForStd(timePointsForStd > starts(i) & timePointsForStd < stops(i) ) ];
        end
        timePointsForStd = temp;
        

        for t = 1:length(timePointsForStd),
            time = timePointsForStd(t);

%             % remove the two highest spikes (assumed to be undesirable outliers)
%             top = sort(thisMeanSubtrSpResp,'descend');
%             if length(top) > 10, 
%                 rejectLine = top(2);
%             else rejectLine = max(top) + 1;
%             end
            
            vals = thisMeanSubtrSpResp(abs(spontTimes - time) < dist); % & thisMeanSubtrSpResp < rejectLine ) ;            
            centeredVals{t,j} = vals;
            std2windows(t,j) = std(vals);
            % note there is some assymetry at the start and end (ie the spontTimes values are
            % mostly from one side of the timepoint)

           % find the interpolated spont resp value close to the timePoint:
            switch fitDegree
                case 0,
                    estSpontResp(t,j) = spline(timePointsForStd, runningMean(:,j), time);
                case 2,
                    estSpontResp(t,j) = B(1)*time.^2 + B(2)*time + B(3);
                case 3,
                    estSpontResp(t,j) = B(1)*time.^3 + B(2)*time.^2 + B(3)*time + B(4);
%                 case 4,
%                     estSpontResp(t,j) = B(1)*time.^4 + B(2)*time.^3 + B(3)*time.^2 + B(4)*time + B(5); 
%                 case 5,
%                     estSpontResp(t,j) = B(1)*time.^5 + B(2)*time.^4 + B(3)*time.^3 + B(4)*time.^2 + B(5)*time + B(6);
            end % switch
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
        
    if plotStdLines,     
        figure,    
        for j = neuronsToPlot
            % index p is used for subplot location:
            p = find(neuronsToPlot == j);
            
            subplot(plotR, plotC, p),
            xlim([0, max(timePoints) ])
            title(num2str(j), 'fontsize',14, 'fontweight','bold')
            if normalizeToInitialFR,
                title(['#' num2str(j) ', initial FR = ' num2str(normFactor(j))] ) %,...
                    % 'fontsize',14,'fontweight','bold'),
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
            plot(timePointsForStd(ord), neg2Std(ord, j),  ':');

            % now plot the stim responses:
            % do them one at a time to get all the color/shape variations:
            for k = 1:length(stimTimes),
                shape = '+'; 
                col = 'k'; % s
                if ismember(k,controlStimInds), col = 'g'; end
                if ismember(k,PnStimInds), col = 'b'; end
                if ismember(k,P9StimInds), col = 'r'; end
                if ismember(k,MnStimInds), col = 'm'; end
                plot(stimTimes(k), stimResponses(j,k),  [col shape] ),
            end  
            yVals = [pos2Std(:,j); stimResponses(j,:)' ];
            ylim([0, max(yVals) ] )
        end % for j = neuronsToPlot
        
        % put the figure title in the lower right
         subplot(plotR,plotC,plotR*plotC),
            title(['BEA. prep ' num2str(ind) ]) % '. delay ' num2str(delay) ' mSec, window ' num2str(length(window)), ' mSec'])
            xlabel('P9 = r, Pn = b, Mn = m, misc = k, ctl = g')
    end

    %% calculate Bhattacharyya distance from a gaussian made with the
        % estimated std for this window (assume mean = 0 since we subtracted the estResp).
        % calc the hypothetical gaussian:

    if calcAndPlotBhatFlag,
         
        for j = neuronsToPlot,
            for t = 1:size(std2windows,1), % should be the same as length(timePointsForStd)
                s = std2windows(t,j);
                theseVals = centeredVals{t,j};
                bins = linspace(-4*s, 4*s, 201);
                halfStep = ( bins(2) - bins(1) ) / 2;
                binCenters = bins(1:length(bins) - 1) + halfStep;
                gaussian = exp( -0.5* (binCenters.^2)/s^2 ); % make the equiv gaussian
                % set the gaussian = 0 for any negative bins, ie bins corresponding to < 0 spontFR:
                gaussian(binCenters < min(vals)) = 0;
                formattedGaussian = [binCenters; gaussian]; % 'bhatt...' needs this gaussian as a 2 x N matrix
                if isempty(theseVals),
                    Bdist(t,j) = -1;
                else
                    Bdist(t,j) = bhattacharyya(theseVals, formattedGaussian, rand(1) < 0.002) ; % show a few spont FR window hists
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
         for j = neuronsToPlot,
             % index p is for subplot location:
            p = find(neuronsToPlot == j);
            
            subplot(plotR, plotC, p),
            this = Bdist(:,j);
            plot(timePointsForStd(this >= 0), this(this >= 0), 'b+' ),
            title(num2str(j))
            grid on   
         end
         % put the figure title in the lower right
         subplot(plotR,plotC,plotR*plotC),
         title(['prep ' num2str(ind) '. Bhattacharyya distance ']) %: delay ' num2str(delay),...
             %' mSec, window ' num2str(length(window)), ' mSec']),     
    end

    %% calculate and plot the running means of the spontaneous responses.

    if plotSpontRunningMedians,
        
        % since the std calc was for spontResp - cubicEstimate, re-do the windowing operation:
        for j = neuronsToPlot,
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

                if isempty(vals),
                    runningMedian(t,j) = -1;
                else
                    runningMedian(t,j) = median(vals);
                end
                if isempty(vals),
                    runningMean(t,j) = -1;
                else
                    runningMean(t,j) = mean(vals);
                end
            end % for t
             % calc the median +/- stds:
            pos1StdMed(:,j) = runningMedian(:,j) + std2windows(:,j);
            neg1StdMed(:,j) = runningMedian(:,j) - std2windows(:,j);
            pos2StdMed(:,j) = runningMedian(:,j) + 2*std2windows(:,j);
            neg2StdMed(:,j) = runningMedian(:,j) - 2*std2windows(:,j);
        end

       
        
        % plot the running medians:
         figure,
         for j = neuronsToPlot,
             % index p is for subplot location:
            p = find(neuronsToPlot == j);

            subplot(plotR, plotC, p),
            plot(timePointsForStd, runningMedian(:,j),'b', 'linewidth',2 ),
            hold on,
          % plot(timePointsForStd, runningMean(:,j),'k','linewidth',2 ),
           % [~,ord] = sort(timePointsForStd, 'ascend');
            plot(timePointsForStd, pos1StdMed(:, j), ':');
            plot(timePointsForStd, pos2StdMed(:, j), ':');
            plot(timePointsForStd, neg1StdMed(:, j), ':');
            plot(timePointsForStd, neg2StdMed(:, j), ':');
            
           % now plot the stim responses:
            % do them one at a time to get all the color/shape variations:
            for k = 1:length(stimTimes),
                shape = '+'; 
                col = 'k'; % s
                if ismember(k,controlStimInds), col = 'g'; end
                if ismember(k,PnStimInds), col = 'b'; end
                if ismember(k,P9StimInds), col = 'r'; end
                if ismember(k,MnStimInds), col = 'm'; end
                plot(stimTimes(k), stimResponses(j,k),  [col shape] ),
            end
            
%             if normalizeToInitialFR,
%                 ylim([0,20])
%             end
            title([ num2str(j) ', initial FR = ' num2str(normFactor(j)) ] )
            grid on 
            yVals = [pos2StdMed(:,j); stimResponses(j,:)' ];             
            ylim([0, max(yVals) ] )
         end
         % put the figure title in the lower right
         subplot(plotR,plotC,plotR*plotC),
            title(['BEA. prep ' num2str(ind) ]) % title(['prep ' num2str(ind) '. running median (blue), mean (black) ' ] ) % : delay ' num2str(delay),...
             %' mSec, window ' num2str(length(window)), ' mSec']),     
    end

    %% normalize FRs by subtracting the estSpontResp and dividing by the std:

    if calcFRMahalDist,
        
        figure,
        
        clear medianStimTimes medianMahalResponses medianStimType medianTrajectoryMahalRespOdor medianTrajectoryMahalRespControl
        clear meanStimTimes meanMahalResponses
        clear stimTypeMatrix medianOctoType
        mahalStimResponses = zeros(size (stimResponses) ); % initialize
        
        spread = 3; % std will be mean of nearest 2*spread + 1 values
        % to max out at mahal = 3: 19 if using dist = exp(mahal) - 1; 70 to 100 if using 1/normal dist prob
        scaleFlag = 1; % 0 means return straight mahal dists. 1 means scale them to reflect probabilities
        upperMahalRespLim = 19.08; % will only matter if scaleFlag = 1.
        controlNormFlag = 1;
        for j = neuronsToPlot
            
           mahalStimResponses(j,:) =...
               calcMahalProbFRs(stimResponses(j,:), stimTimes, (estSpontResp(:,j))', (std2windows(:,j))',...
               timePointsForStd, spread, upperMahalRespLim,scaleFlag, controlNormFlag, controlStimInds);
           
            % calc medians in 2 directions:
            % the most relevant is the second, ie finding a platonic form
            % of response trajectory to a set of fast puffs.
            
            % 1a. make a vector that labels odor vs control:
            stimTypeMatrix = zeros(size(stim(:) ) );
            stimTypeMatrix(miscStimInds) = 1;
            stimTypeMatrix(MnStimInds) = -2;
            stimTypeMatrix(PnStimInds) = -1;
            stimTypeMatrix(controlStimInds) = 0;
            stimTypeMatrix(P9StimInds) = 2;
            stimTypeMatrix = reshape(stimTypeMatrix, size(stim) );
            % now stimTypeMatrix has entries = 1 (odor), 0 (control), -2
            % (spont),  -1 (air), 2 (pheromone)
            
            octoMatrix = zeros( size ( stim(:) ) );
            octoMatrix(preOctoInds) = -1;
            octoMatrix(postOctoInds) = 1;
            medianOctoType(j,:) = median( reshape (octoMatrix, size(stim) ) ); % 1 = odor, 0 = control
            
            % 1b. calc medians per set of puffs. ie combine a set of 5 or so closely-spaced puffs into one response:
            medianStimTimes(j,:) = median(reshape(stimTimes,size(stim) ) );
            meanStimTimes(j,:) = mean(reshape(stimTimes, size(stim) ) );
            medianMahalResponses(j,:) = median( reshape( (mahalStimResponses(j,:))', size(stim) ) );
            meanMahalResponses(j,:) = mean( reshape( (mahalStimResponses(j,:))', size(stim) ) );
            
            % 2. calculate the platonic ideal of a response trajectory given a set of puffs:
            FRmatrix = reshape( (mahalStimResponses(j,:))', size(stim) );
            [ trajFRsOdor, trajFRsControl, odorStd ] = calcTrajectoryResponseMedian( FRmatrix, miscStimInds );
            medianTrajectoryMahalRespOdor (j,: ) = trajFRsOdor;
            medianTrajectoryMahalRespControl(j,:) = trajFRsControl; 

            % index p is used for subplot location:
            p = find(neuronsToPlot == j);

            subplot(plotR, plotC, p),

            title(num2str(j), 'fontsize',14, 'fontweight','bold')
            if 1, % normalizeToInitialFR,
                title(['#' num2str(j) ',initial FR = ' num2str(normFactor(j)) ] ) %,...
                    % 'fontsize',14,'fontweight','bold'),
            end
            hold on,
            if plotMedianTrajectories, % case: we want to see the platonic time-series response to a set of fast puffs
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
                % now plot the stim responses:
                % do them one at a time to get all the color/shape variations:
                for k = 1:length(stimTimes),
                    shape = '+'; 
                    col = 'k'; % misc
                    if ismember(k,controlStimInds), col = 'g'; end
                    if ismember(k,PnStimInds), col = 'b'; end
                    if ismember(k,P9StimInds), col = 'r'; end
                    if ismember(k,MnStimInds), col = 'm'; end
                    plot(stimTimes(k), mahalStimResponses(j,k),  [col shape] ),                    
                    % plot median responses:
                    [~, equivInd] = ind2sub(size(stim),k);
                    plot(medianStimTimes(j,equivInd), medianMahalResponses(j,equivInd), [col 'o'] ),
                end
               % % plot mean stim responses:
               % plot(meanStimTimes(j,:), meanMahalResponses(j,:),'k*')
            end
            
%             plot( [ timePointsForStd(closestInds(j,2) - spread), timePointsForStd (closestInds(j,2) + spread ) ], [1 1],'k+')
%             plot( timePointsForStd(closestInds(j,2) - spread), 1,'k+'),
%             plot( timePointsForStd (closestInds(j,2) + spread ), 1,'k+')
            grid on

        end % for j = neuronsToPlot

            % put the figure title in the lower right
             subplot(plotR,plotC,plotR*plotC),
            title(['BEA. prep ' num2str(ind) ]) %'. delay ' num2str(delay) ' mSec, window ' num2str(length(window)), ' mSec'])
            xlabel('P9 = r, Pn = b, Mn = m, misc = k, ctl = g')
    end   
    
%     % save results for future use:
%     medianMahalResults(ind).medianStimTimes = medianStimTimes;
%     medianMahalResults(ind).medianMahalResponses = medianMahalResponses;
%     medianMahalResults(ind).medianStimType = medianOctoType;
%     medianMahalResults(ind).initialFR = normFactor;
%     medianMahalResults(ind).medianTrajMahalRespOdor = medianTrajectoryMahalRespOdor;
%     medianMahalResults(ind).medianTrajMahalRespControl = medianTrajectoryMahalRespControl;
            
            
end % for ind

% save medianMahalResultsForOdorNoReward medianMahalResults
 
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
                  
            
            
            
            
            
            
            
            
    
    