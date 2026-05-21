clear all
close all
clc

%% Initialization of some values

fsEMG = 4096 ;
ptilethresh = 98 ;

%% Load Data

dataFolder = '/Users/amyneeson/Desktop/PHD/Acute Study 1/Data/Post CKC/7' ;
content = dir(dataFolder) ;
preFiles = content(contains({content.name},'Pre')) ;
postFiles = content(contains({content.name},'Post')) ;

%% Go Through Files

allClusters = {} ;
allCouples = {} ;

preNames = {} ;
postNames = {} ;

preFilesSave = {} ;
postFilesSave = {} ;

numFiles = 1 ;

% Code first computes all possible clusters then shows them to user to
% track accross trials.

for scanContentPre = 1:size(preFiles,1)
    progBar = waitbar(0,'Tracking') ;

    %% Pre Fatigue

    preData = load([dataFolder '/' preFiles(scanContentPre).name]) ;
    preTxx = preData.Txx ; % Txx is a 3D matrix that contains MUAP shapes in the format [Samples X Channels X NumberOfUnits]
    NbMUsPre = size(preTxx,3) ;

    for scanContentPost = 1:size(postFiles,1)

        %% Post Fatigue

        postData = load([dataFolder '/' postFiles(scanContentPost).name]) ;
        postTxx = postData.Txx ; % Txx is a 3D matrix that contains MUAP waveforms in the format [Samples X Channels X NumberOfUnits]
        NbMUsPost = size(postTxx,3) ;

        %% Tracking

        startEndPre = findStartEnd_Jer(preTxx,fsEMG) ; % Code identifies start and end of MUAP waveform
        startEndPost = findStartEnd_Jer(postTxx,fsEMG) ;
        whichChansToUse = true(1,size(preTxx,2)) ; % If some channels need to be discard, specify here
        [ClustCompareNew,MUCouples,~,~,allHotSpots] = ClusterCompareTrialsHighDim(cat(3,preTxx,postTxx),cat(1,startEndPre,startEndPost),1./fsEMG,[],ptilethresh,size(preTxx,3),0,whichChansToUse) ;
        MUCouples(:,2) = MUCouples(:,2) - NbMUsPre ;

        allClusters{numFiles,1} = ClustCompareNew ;
        allCouples{numFiles,1} = MUCouples ;

        preNames{numFiles,1} = preFiles(scanContentPre).name(1:end-3) ;
        postNames{numFiles,1} = postFiles(scanContentPost).name(1:end-3) ;

        preFilesSave{numFiles,1} = preTxx ;
        postFilesSave{numFiles,1} = postTxx ;

        numFiles = numFiles + 1 ;

        waitbar(scanContentPost/size(postFiles,1),progBar,'Tracking') ;

    end
end

for scanClusters = 1:size(allClusters,1)
    %% Plotting Results

    ClustCompareNew = allClusters{scanClusters,1} ;
    MUCouples = allCouples{scanClusters,1} ;

    % preTxx = preFilesSave{scanClusters,1} ;
    % postTxx = postFilesSave{scanClusters,1} ;

    Inds1 = num2str(MUCouples(:,1)) ;
    Inds1Num = MUCouples(:,1) ;

    Inds2 = num2str(MUCouples(:,2)) ;
    Inds2Num = MUCouples(:,2) ;
    commas = repmat(',',size(Inds1,1),1) ;
    AllString = [Inds1,commas,Inds2] ;

    featuresList = [8 11 12] ;

    [PCACoeffs,PCAScoresCompareNew,~,~,~,~] = pca((ClustCompareNew(:,featuresList)./std(ClustCompareNew(:,featuresList),[],1)),'Centered','on','NumComponents',size(ClustCompareNew(:,featuresList),2)) ;

    invCoeffs = inv(PCACoeffs) ;
    zeroPosInPCA = (-mean((ClustCompareNew(:,featuresList)./std(ClustCompareNew(:,featuresList),[],1))))*invCoeffs' ;
    zeroPosInPCA = zeroPosInPCA(1:3) - mean(PCAScoresCompareNew(:,1:3),1) ;

    newRefPoints = PCAScoresCompareNew(:,1:3) - prctile(PCAScoresCompareNew(:,1:3),75) ;
    thetaAngle = acosd(dot(newRefPoints,repmat(zeroPosInPCA,[size(newRefPoints,1) 1]),2)./(vecnorm(newRefPoints')'.*vecnorm(repmat(zeroPosInPCA,[size(newRefPoints,1) 1])')')) ;
    coneTowardsZero = dot(newRefPoints,repmat(zeroPosInPCA,[size(newRefPoints,1) 1]),2) > 0 & thetaAngle <= 45 ;

    distFromZero = sqrt(sum((PCAScoresCompareNew(:,1:3)-zeroPosInPCA).^2,2)) ;
    distFromZeroZScore = (distFromZero - mean(distFromZero))./std(distFromZero) ;

    outliersFirstPass = coneTowardsZero & distFromZeroZScore < prctile(distFromZeroZScore,20) ;

    removeCone = coneTowardsZero & distFromZeroZScore < prctile(distFromZeroZScore,20) ;
    sortedOutliers = sortrows(cat(2,Inds1Num(removeCone),Inds2Num(removeCone),distFromZero(removeCone),zeros(sum(removeCone),1)),3) ;

    bestCouples = sortedOutliers(1,1:2) ;
    bestCouplesPos = find(distFromZero == sortedOutliers(1,3)) ;
    removeCouples = [] ;
    for scanCouples = 2:size(sortedOutliers,1)
        if ~((sum(bestCouples(:,1) == sortedOutliers(scanCouples,1)) == 1) || (sum(bestCouples(:,2) == sortedOutliers(scanCouples,2)) == 1))
            bestCouples = cat(1,bestCouples,sortedOutliers(scanCouples,1:2)) ;
            bestCouplesPos  = cat(1,bestCouplesPos,find(distFromZero == sortedOutliers(scanCouples,3))) ;
            sortedOutliers(scanCouples,4) = 1 ;
        else
            removeCouples = cat(1,removeCouples,find(distFromZero == sortedOutliers(scanCouples,3))) ;
        end
    end

    removeConeDist = coneTowardsZero ;
    sortedOutliersDist = sortrows(cat(2,Inds1Num(removeConeDist),Inds2Num(removeConeDist),distFromZero(removeConeDist),zeros(sum(removeConeDist),1)),3) ;

    bestCouplesDist = sortedOutliersDist(1,1:2) ;
    secondBestMatch =  [] ;
    secondBestMatchPos = [] ;
    removeCouplesDist = [] ;
    for scanCouples = 2:size(sortedOutliersDist,1)
        if ~((sum(bestCouplesDist(:,1) == sortedOutliersDist(scanCouples,1)) == 1) || (sum(bestCouplesDist(:,2) == sortedOutliersDist(scanCouples,2)) == 1))
            bestCouplesDist = cat(1,bestCouplesDist,sortedOutliersDist(scanCouples,1:2)) ;
        elseif isempty(secondBestMatch)
            secondBestMatch = cat(1,secondBestMatch,sortedOutliersDist(scanCouples,1:2)) ;
            secondBestMatchPos = cat(1,secondBestMatchPos,find(distFromZero == sortedOutliersDist(scanCouples,3))) ;
        elseif (sum(secondBestMatch(:,1) == sortedOutliersDist(scanCouples,1)) < 1) && (sum(secondBestMatch(:,2) == sortedOutliersDist(scanCouples,2)) < 1)
            secondBestMatch = cat(1,secondBestMatch,sortedOutliersDist(scanCouples,1:2)) ;
            secondBestMatchPos = cat(1,secondBestMatchPos,find(distFromZero == sortedOutliersDist(scanCouples,3))) ;
        else
            removeCouplesDist = cat(1,removeCouplesDist,find(distFromZero == sortedOutliersDist(scanCouples,3))) ;
        end
    end

    if ~isempty(secondBestMatch)
        distancesSecondMatches = distFromZero(secondBestMatchPos) ;
        confVal = prctile(distancesSecondMatches,10) ;
        tooFar = distFromZero(bestCouplesPos) > confVal ;
        removeCouples = cat(1,removeCouples,bestCouplesPos(tooFar)) ;
        bestCouplesPos(tooFar) = [] ;
    end

    distFromZeroStep2_Removed = distFromZero ;
    distFromZeroStep2_Removed(removeCouples) = [] ;
    coneTowardsZero_Removed = coneTowardsZero ;
    coneTowardsZero_Removed(removeCouples) = [] ;
    distFromZeroStep2ZScore = ((distFromZeroStep2_Removed - mean(distFromZeroStep2_Removed))./std(distFromZeroStep2_Removed)) ;

    %% REMOVE THE ZERO VALUES

    outliersSecondPass = cat(2,outliersFirstPass,Inds1Num,Inds2Num) ;
    outliersSecondPass(removeCouples,:) = [] ;

    newThresh = 100*((size(preTxx,3)+size(postTxx,3))/(2*size(distFromZeroStep2ZScore,1))) ;
    outliersSecondPass(:,1) = outliersSecondPass(:,1) & distFromZeroStep2ZScore < prctile(distFromZeroStep2ZScore(distFromZeroStep2ZScore ~=0),newThresh) ;
    PCAScoresCut = PCAScoresCompareNew ;
    PCAScoresCut(removeCouples,:) = [] ;

    clustCenter = mean(PCAScoresCompareNew(:,1:3),1) ;
    clusterFig = figure('WindowState','maximized') ;
    hold on
    plot3(PCAScoresCompareNew(:,1),PCAScoresCompareNew(:,2),PCAScoresCompareNew(:,3),'ro')
    plot3(PCAScoresCompareNew(bestCouplesPos,1),PCAScoresCompareNew(bestCouplesPos,2),PCAScoresCompareNew(bestCouplesPos,3),'bx')
    plot3(PCAScoresCompareNew(outliersFirstPass,1),PCAScoresCompareNew(outliersFirstPass,2),PCAScoresCompareNew(outliersFirstPass,3),'g+')
    plot3(PCAScoresCompareNew(removeCouples,1),PCAScoresCompareNew(removeCouples,2),PCAScoresCompareNew(removeCouples,3),'rd')
    plot3(PCAScoresCut(outliersSecondPass(:,1) == 1,1),PCAScoresCut(outliersSecondPass(:,1) == 1,2),PCAScoresCut(outliersSecondPass(:,1) == 1,3),'ks','MarkerSize',15)
    text(PCAScoresCompareNew(:,1),PCAScoresCompareNew(:,2),PCAScoresCompareNew(:,3),AllString) ;
    % close all ;

    % User can press "a" for automatic thresholding at 95% confidence
    % interval or press "m" for manual selection of pairs. If pressing "m"
    % pairs need to be specified in the format [1 2;3 4;5 6...].

    keyPressed = 0 ;
    while ~keyPressed
        waitforbuttonpress() ;
        whoPressed = clusterFig.CurrentCharacter ;
        if strcmp(whoPressed,'m')
            keyPressed = 1 ;
            autoSelect = 0 ;
        elseif strcmp(whoPressed,'a')
            keyPressed = 1 ;
            autoSelect = 1 ;
        end
    end

    if autoSelect == 0
        prompt = 'Which MUs do you want to track ?' ;
        trackList = inputdlg(prompt,'MUs To Track') ;
        trackThese = str2num(trackList{1,1}) ;
    else
        trackThese = outliersSecondPass(outliersSecondPass(:,1) == 1,2:3) ;
    end
    close ;

    if ~isempty(trackThese)
        %% Check Couples Found

        MUSet1 = trackThese(:,1) ;
        MUSet2 = trackThese(:,2) ;

        couplesAgreed = [] ;

        for scanMUs = 1:length(MUSet1)
     
            gridSTA1 = permute(reshape(preTxx(:,:,MUSet1(scanMUs)).',5,13,[]),[2 1 3]) ;
            gridSTA2 = permute(reshape(postTxx(:,:,MUSet2(scanMUs)).',5,13,[]),[2 1 3]) ;

            plotGridComparison(gridSTA1,gridSTA2,'Grid',1,1) ;

            pause(1)
            agreedOrNot = questdlg('Is this a couple','Couple','Yes','No','Yes') ;
            if length(agreedOrNot) == 3
                couplesAgreed = [couplesAgreed; scanMUs] ;
            end
            close ;
        end

        MUSet1 = MUSet1(couplesAgreed) ;
        MUSet2 = MUSet2(couplesAgreed) ;

        fileNamePre = preNames{scanClusters} ;
        fileNamePost = postNames{scanClusters} ;
        save([dataFolder '/dataTracked_' fileNamePre '_' fileNamePost '.mat'],"MUSet1","MUSet2","NbMUsPre","NbMUsPost") ;

    end
end

