function [ClustCompare,MU_Couple_Orig,ScaleFactor,AllRatios,allHotSpots]=ClusterCompareTrialsHighDim(Txx_Smooth,StartEndAllCNoDown,dt,StabilityScore,ptile,MUsInSet1,simData,whichChansToUse)

Txx_All = Txx_Smooth ;
Txx_Smooth = Txx_Smooth(:,whichChansToUse,:) ;

GoodMus = 1:size(Txx_Smooth,3) ;

% SmoothingIx = 0.01 ;
% might have to change to match the delsys
% Matrices for the clusters

NbMus = size(Txx_Smooth,3) ;
MeasureDistMin = zeros(NbMus,NbMus,4) ;
MeasureDistMinNew= zeros(NbMus,NbMus) ;
MeasureDistMinZtrans= zeros(NbMus,NbMus) ;
MeasureVeloMin = zeros(NbMus,NbMus) ;
MeasureAccelMin = zeros(NbMus,NbMus) ;
MeasureCurvatureMin = zeros(NbMus,NbMus,4) ;
MeasureVectorMin = zeros(NbMus,NbMus,4) ;
KeepMU_CVchange = zeros(NbMus,NbMus) ;
MeasureCMCMax = zeros(NbMus,NbMus,4) ;

AngleVeloMin = zeros(NbMus,NbMus) ;
AngleDistMin = zeros(NbMus,NbMus) ;
AngleAccelMin = zeros(NbMus,NbMus) ;
AngleCurvMin = zeros(NbMus,NbMus,4) ;
AngleVectorMin = zeros(NbMus,NbMus,4) ;

MaxAmps = zeros(NbMus,NbMus,4) ;
DistSeuc= zeros(NbMus,NbMus);
DistCos= zeros(NbMus,NbMus);
ProSimilarity= zeros(NbMus,NbMus);

allHotSpots = cell(NbMus,NbMus) ;

for scanMUs = 1:NbMus
    if scanMUs < MUsInSet1 + 1
        secondSet{scanMUs} =  [MUsInSet1+1:NbMus] ;
    else
        secondSet{scanMUs} = [1:MUsInSet1] ;
    end
end

%% Raw MUAPs Clustering

f = @(x,a,b)sum((vecnorm((a./repmat(x(1),size(a,1),size(a,2)))-b,2,2).*(vecnorm(a,2,2).*vecnorm(b,2,2))./max(vecnorm(a,2,2).*vecnorm(b,2,2)))./mean([max(vecnorm((a./repmat(x(1),size(a,1),size(a,2))),2,2)),max(vecnorm(b,2,2))])) ;
x0=0.05;
x0end=20;

for Mu1Ind = 1:NbMus
    for Mu2Ind = secondSet{Mu1Ind}

        Mu1Avg = Txx_Smooth(:,:,Mu1Ind) ;
        Mu2Avg = Txx_Smooth(:,:,Mu2Ind) ;

        if simData && size(whichChansToUse,2) > 4

            rmsMu1 = rms(Mu1Avg,1) ;
            rmsMu2 = rms(Mu2Avg,1) ;

            listChans = find(whichChansToUse) ;
            listChans(~(rmsMu1 > max(rmsMu1)/10 | rmsMu2 > max(rmsMu2)/10 | rmsMu1 > max(rmsMu1)/10 | rmsMu2 > max(rmsMu2)/10)) = [] ;
            Mu1Avg = Mu1Avg(:,rmsMu1 > max(rmsMu1)/10 | rmsMu2 > max(rmsMu2)/10,:) ;
            Mu2Avg = Mu2Avg(:,rmsMu1 > max(rmsMu1)/10 | rmsMu2 > max(rmsMu2)/10,:) ;

        else 

            listChans = 1:size(Txx_Smooth,2) ;

        end
        
        %% Scaling Motor units

        AllDistRotateStart = zeros(size(Mu1Avg,1),1) ;
        Mu1AvgT=[];Mu2AvgT=[];
        Window = zeros(size(Mu1Avg,1),1);
        if abs(StartEndAllCNoDown(Mu2Ind,2)-StartEndAllCNoDown(Mu2Ind,1))<75||abs(StartEndAllCNoDown(Mu1Ind,2)-StartEndAllCNoDown(Mu1Ind,1))<75
            Window(round(((size(Mu1Avg,1)/2)-0.005/dt):((size(Mu1Avg,1)/2)+0.005/dt))) = tukeywin(length(round(((size(Mu1Avg,1)/2)-0.005/dt):((size(Mu1Avg,1)/2)+0.005/dt))),0.5) ;
        else
            Window(round(((size(Mu1Avg,1)/2)-0.01/dt):((size(Mu1Avg,1)/2)+0.01/dt))) = tukeywin(length(round(((size(Mu1Avg,1)/2)-0.01/dt):((size(Mu1Avg,1)/2)+0.01/dt))),0.5) ;
        end
        Mu1AvgT=Mu1Avg.*repmat(Window,1,size(Mu1Avg,2));
        Mu2AvgT=Mu2Avg.*repmat(Window,1,size(Mu2Avg,2));
        Window = [];
        for RotateSmooth = 1:size(Mu1Avg,1)
            Mu2RotateSmooth = repmat(Mu2AvgT,2,1) ;
            Mu2RotateSmooth = Mu2RotateSmooth(RotateSmooth:RotateSmooth+size(Mu2AvgT,1)-1,:) ;
            AllDistRotateStart(RotateSmooth) = nanmean(sqrt(sum(((Mu1AvgT-Mu2RotateSmooth)).^2,2))) ;
        end

        % Restrict the amount it can Rotate by
        [~,RotateSmooth] = min(AllDistRotateStart([1:0.005/dt round(((size(Mu1Avg,1)-1)-0.005/dt):((size(Mu1Avg,1))))])) ;
        winC=[];
        winC=floor(RotateSmooth/2);

        if RotateSmooth > 0.005/dt
            RotateSmooth = RotateSmooth + (round((size(Mu1Avg,1)-1)-0.005/dt)-1)-0.005/dt ;
            winC=floor((RotateSmooth-size(Mu1Avg,1))/2);
        end

        Mu2RotateSmooth = repmat(Mu2Avg,2,1) ;
        Mu2RotateSmooth = Mu2RotateSmooth(round(RotateSmooth):round(RotateSmooth+size(Mu2Avg,1)-1),:) ;

        % Only take the centre of the window
        % take a narrow window if durations are narrow
        if abs(StartEndAllCNoDown(Mu2Ind,2)-StartEndAllCNoDown(Mu2Ind,1))<75||abs(StartEndAllCNoDown(Mu1Ind,2)-StartEndAllCNoDown(Mu1Ind,1))<75
            fun = @(x)f(x,Mu1Avg(round(((size(Mu1Avg,1)/2-winC)-0.0025/dt):((size(Mu1Avg,1)/2-winC)+0.0025/dt)),:),Mu2RotateSmooth(round(((size(Mu1Avg,1)/2-winC)-0.0025/dt):((size(Mu1Avg,1)/2-winC)+0.0025/dt)),:)) ;
        else
            fun = @(x)f(x,Mu1Avg(round(((size(Mu1Avg,1)/2-winC)-0.01/dt):((size(Mu1Avg,1)/2-winC)+0.01/dt)),:),Mu2RotateSmooth(round(((size(Mu1Avg,1)/2-winC)-0.01/dt):((size(Mu1Avg,1)/2-winC)+0.01/dt)),:)) ;
        end
        options = optimset('Display','off') ;
        [scaleNorm,~,exitflag] = fminbnd(fun,x0,x0end,options) ;

        Mu1Avg_NoScale = Mu1Avg ;
        Mu2Avg_NoScale = Mu2Avg ;

        if exitflag == 0 || (scaleNorm < 0.2 || scaleNorm > 5) || (max(vecnorm(Mu1Avg'./scaleNorm))/max(vecnorm(Mu1Avg')) < 0.2 || max(vecnorm(Mu1Avg'./scaleNorm))/max(vecnorm(Mu1Avg')) > 10)

            NormeMu1 = sqrt(sum(Mu1Avg.^2,2)) ;
            NormeMu1 = max(NormeMu1) ;
            Mu1Avg = Mu1Avg./NormeMu1 ;

            NormeMu2 = sqrt(sum(Mu2Avg.^2,2)) ;
            NormeMu2 = max(NormeMu2) ;
            Mu2Avg = Mu2Avg./NormeMu2 ;

        else

            Mu1Avg = Mu1Avg./scaleNorm ;
            NormeMu1 = sqrt(sum(Mu1Avg.^2,2)) ;
            NormeMu1 = max(NormeMu1) ;

            NormeMu2 = sqrt(sum(Mu2Avg.^2,2)) ;
            NormeMu2 = max(NormeMu2) ;
            Mu1Avg = Mu1Avg./max([NormeMu1 NormeMu2]) ;
            Mu2Avg = Mu2Avg./max([NormeMu1 NormeMu2]) ;

        end

        %% Computing Weights

        rmsMu1 = rms(Mu1Avg(StartEndAllCNoDown(Mu1Ind,1):StartEndAllCNoDown(Mu1Ind,2),:)) ;
        rmsMu2 = rms(Mu2Avg(StartEndAllCNoDown(Mu2Ind,1):StartEndAllCNoDown(Mu2Ind,2),:)) ;

        hotSpot1 = rmsMu1 > median(rmsMu1) ;
        hotSpot2 = rmsMu2 > median(rmsMu2) ;

        for scanSpots = 1:2
            if scanSpots == 1
                currentHotSpot = hotSpot1 ;
            else
                currentHotSpot = hotSpot2 ;
            end

            currentHotSpot_GridForm = zeros(14,9) ;
            chanCount = 1 ;
            hotSpotCount = 1 ;

            for scanRows = 1:size(currentHotSpot_GridForm,1)
                for scanCols = 1:size(currentHotSpot_GridForm,2)
                    if any(listChans == chanCount)
                        currentHotSpot_GridForm(scanRows,scanCols) =  currentHotSpot(hotSpotCount) ;
                        hotSpotCount = hotSpotCount + 1;
                    end
                    chanCount = chanCount + 1 ;
                end
            end

            currentHotSpot_GridForm = bwlabel(currentHotSpot_GridForm,4) ;
            if max(currentHotSpot_GridForm,[],"all") > 1
                objSize = zeros(max(currentHotSpot_GridForm,[],"all"),1) ;
                for scanObjs = 1:max(currentHotSpot_GridForm,[],"all")
                    objSize(scanObjs) = sum(currentHotSpot_GridForm == scanObjs,"all") ;
                end
                [~,biggestObj] = max(objSize) ;
                currentHotSpot_GridForm(currentHotSpot_GridForm ~= biggestObj) = 0 ;
                currentHotSpot_GridForm = currentHotSpot_GridForm./max(currentHotSpot_GridForm,[],"all") ;
            end

            currentHotSpot_GridForm = reshape(permute(currentHotSpot_GridForm,[2 1 3]),126,[]).';
            if scanSpots == 1
                hotSpot1 = logical(currentHotSpot_GridForm(listChans)) ;
            else
                hotSpot2 = logical(currentHotSpot_GridForm(listChans)) ;
            end
        end

        hotSpotFinal = hotSpot1 | hotSpot2 ;
        hotSpotWeights = hotSpot1 & hotSpot2 ;

        allHotSpots{Mu1Ind,Mu2Ind} = hotSpotFinal ;

        Mu1Avg = Mu1Avg(:,hotSpotFinal) ;
        Mu2Avg = Mu2Avg(:,hotSpotFinal) ;

        Mu1Avg_NoScale = Mu1Avg_NoScale(:,hotSpotFinal) ;
        Mu2Avg_NoScale = Mu2Avg_NoScale(:,hotSpotFinal) ;


        numChans(Mu1Ind,Mu2Ind) = 100*sum(hotSpotWeights)/sum(hotSpotFinal) ;
       
        %% CV Changes and Rotate Stuff

        StartEnd = [StartEndAllCNoDown(Mu1Ind,:);
                    StartEndAllCNoDown(Mu2Ind,:)] ;

        [CVRatio,whoToResample] = scanCVs(Mu1Avg,Mu2Avg) ;

        CVRatio = round(CVRatio,4) ;

        AllRatios(Mu1Ind,Mu2Ind) = CVRatio ;

        % Change to 1 to not allow a CV change, typically 1.5
        ratioLimit = 2 ;

        if CVRatio <= ratioLimit

            if whoToResample == 1
                flagCV=2;
                Center = round(mean(StartEnd(1,:))*CVRatio) ;

                MuToResample = Mu1Avg ;
                MuToResample = resample(MuToResample,round(CVRatio*10000),10000) ;
                MuToResample = MuToResample(Center - (round(mean(StartEnd(1,:)))-1):Center + round(size(Mu1Avg,1)-round(mean(StartEnd(1,:)))),:) ;

                MuToResample_NoScale = Mu1Avg_NoScale ;
                MuToResample_NoScale = resample(MuToResample_NoScale,round(CVRatio*10000),10000) ;
                MuToResample_NoScale = MuToResample_NoScale(Center - (round(mean(StartEnd(1,:)))-1):Center + round(size(Mu1Avg_NoScale,1)-round(mean(StartEnd(1,:)))),:) ;

                for Rotate = 1:size(Mu2Avg,1)

                    MuRotate = repmat(Mu1Avg,2,1) ;
                    MuRotate = MuRotate(Rotate:Rotate+size(Mu2Avg,1)-1,:) ;
                    DistBefore(Rotate) = nanmean(sqrt(sum(((Mu2Avg-MuRotate)).^2,2))) ;

                    MuRotate = repmat(MuToResample,2,1) ;
                    MuRotate = MuRotate(Rotate:Rotate+size(Mu2Avg,1)-1,:) ;
                    DistAfter(Rotate) = nanmean(sqrt(sum(((Mu2Avg-MuRotate)).^2,2))) ;

                end
            else
                flagCV = 1 ; % MU to keep for scaling
                Center = round(mean(StartEnd(2,:))*CVRatio) ;

                MuToResample = Mu2Avg ;
                MuToResample = resample(MuToResample,round(CVRatio*10000),10000) ;
                MuToResample = MuToResample(Center - (round(mean(StartEnd(2,:)))-1):Center + round(size(Mu2Avg,1)-round(mean(StartEnd(2,:)))),:) ;

                MuToResample_NoScale = Mu2Avg ;
                MuToResample_NoScale = resample(MuToResample_NoScale,round(CVRatio*10000),10000) ;
                MuToResample_NoScale = MuToResample_NoScale(Center - (round(mean(StartEnd(2,:)))-1):Center + round(size(Mu2Avg,1)-round(mean(StartEnd(2,:)))),:) ;

                for Rotate = 1:size(Mu1Avg,1)

                    MuRotate = repmat(Mu2Avg,2,1) ;
                    MuRotate = MuRotate(Rotate:Rotate+size(Mu1Avg,1)-1,:) ;
                    DistBefore(Rotate) = nanmean(sqrt(sum(((Mu1Avg-MuRotate)).^2,2))) ;

                    MuRotate = repmat(MuToResample,2,1) ;
                    MuRotate = MuRotate(Rotate:Rotate+size(Mu1Avg,1)-1,:) ;
                    DistAfter(Rotate) = nanmean(sqrt(sum(((Mu1Avg-MuRotate)).^2,2))) ;

                end
            end

            if whoToResample == 1 && min(DistBefore) > min(DistAfter)

                Mu1Avg = MuToResample ;
                Mu1Avg_NoScale = MuToResample_NoScale ;
                NewWinLength = CVRatio*diff(StartEnd(1,:)) ;
                if ~(NewWinLength > size(Mu2Avg,1))
                    StartEnd(1,:) = [(round(mean(StartEnd(1,:))))-ceil(NewWinLength/2) (round(mean(StartEnd(1,:))))+floor(NewWinLength/2)] ;

                    if StartEnd(1,2) > size(Mu1Avg,1)
                        TooMuch = StartEnd(1,2) - size(Mu1Avg,1) ;
                        StartEnd(1,:) = [max(1,StartEnd(1,1)-TooMuch) size(Mu1Avg,1)] ;
                    elseif StartEnd(1,1) < 1
                        TooMuch = 1 - StartEnd(1,1) ;
                        StartEnd(1,:) = [1 min(size(Mu1Avg,1),StartEnd(1,2)+TooMuch)] ;
                    end
                end

            elseif min(DistBefore) > min(DistAfter)

                Mu2Avg = MuToResample ;
                Mu2Avg_NoScale = MuToResample_NoScale ;
                NewWinLength = CVRatio*diff(StartEnd(2,:)) ;
                if ~(NewWinLength > size(Mu2Avg,1))
                    StartEnd(2,:) = [(round(mean(StartEnd(2,:))))-floor(NewWinLength/2) (round(mean(StartEnd(2,:))))+ceil(NewWinLength/2)] ;

                    if StartEnd(2,2) > size(Mu1Avg,1)
                        TooMuch = StartEnd(2,2) - size(Mu1Avg,1) ;
                        StartEnd(2,:) = [max(1,StartEnd(2,1)-TooMuch) size(Mu1Avg,1)] ;
                    elseif StartEnd(2,1) < 1
                        TooMuch = 1 - StartEnd(2,1) ;
                        StartEnd(2,:) = [1 min(size(Mu1Avg,1),StartEnd(2,2)+TooMuch)] ;
                    end
                end

            end
        else
            flagCV=0;
        end

        WinSizes = StartEnd(:,2) - StartEnd(:,1) + 1 ;

        [~,Smallest] = min(WinSizes) ;

        if WinSizes(1) ~= WinSizes(2)

            AllDistRotateOrig = zeros(size(Mu1Avg,1),1) ;
            SizeDiff = abs(WinSizes(1) - WinSizes(2)) ;

            if Smallest==1
                Mu2nd=Mu1Avg;
                MuOrig=Mu2Avg;
            else
                Mu2nd=Mu2Avg;
                MuOrig=Mu1Avg;
            end

            for Rotate = 1:size(MuOrig,1)
                MuRotate = repmat(Mu2nd,2,1) ;
                MuRotate = MuRotate(Rotate:Rotate+size(Mu2nd,1)-1,:) ;
                AllDistRotateOrig(Rotate) = nanmean(sqrt(sum(((MuOrig-MuRotate)).^2,2))) ;
            end

            CentersDiff = mean(StartEnd(Smallest,:)) - mean(StartEnd(mod(Smallest,2)+1,:)) ;
            if CentersDiff > 0
                SearchArea = max(1,round(CentersDiff-SizeDiff/2)):min(size(Mu1Avg,1),round(CentersDiff+SizeDiff/2)) ;
            else
                SearchArea = max(1,round(size(Mu1Avg,1)+CentersDiff-SizeDiff/2)):min(size(Mu1Avg,1),round(size(Mu1Avg,1)+CentersDiff+SizeDiff/2)) ;
            end

            [~,RotateOrig] = min(AllDistRotateOrig(SearchArea)) ;

            if abs(RotateOrig - CentersDiff) > SizeDiff/2
                Direction = sign(RotateOrig - CentersDiff) ;
                if ~(SizeDiff*(Direction+1)/2+StartEnd(Smallest,2) > size(Mu1Avg,1))
                    StartEnd(Smallest,:) = [SizeDiff*(Direction-1)/2+StartEnd(Smallest,1) SizeDiff*(Direction+1)/2+StartEnd(Smallest,2)] ;
                else
                    StartEnd(Smallest,:) = [SizeDiff*(Direction-1)/2+StartEnd(Smallest,1)-(SizeDiff*(Direction+1)/2+StartEnd(Smallest,2)-size(Mu1Avg,1)) size(Mu1Avg,1)] ;
                end
            else
                Direction = sign(RotateOrig - CentersDiff) ;
                if ~(abs(RotateOrig - CentersDiff)*(Direction+1)+StartEnd(Smallest,2) > size(Mu1Avg,1))
                    StartEnd(Smallest,:) = [abs(RotateOrig - CentersDiff)*(Direction-1)+StartEnd(Smallest,1) abs(RotateOrig - CentersDiff)*(Direction+1)+StartEnd(Smallest,2)] ;
                else
                    StartEnd(Smallest,:) = [abs(RotateOrig - CentersDiff)*(Direction-1)+StartEnd(Smallest,1)-(abs(RotateOrig - CentersDiff)*(Direction+1)+StartEnd(Smallest,2)-size(Mu1Avg,1)) size(Mu1Avg,1)] ;
                end
                if ~(floor((SizeDiff-2*abs(RotateOrig - CentersDiff))/2)+StartEnd(Smallest,2) > size(Mu1Avg,1))
                    StartEnd(Smallest,:) = [-ceil((SizeDiff-2*abs(RotateOrig - CentersDiff))/2)+StartEnd(Smallest,1) floor((SizeDiff-2*abs(RotateOrig - CentersDiff))/2)+StartEnd(Smallest,2)] ;
                else
                    StartEnd(Smallest,:) = [-ceil((SizeDiff-2*abs(RotateOrig - CentersDiff))/2)+StartEnd(Smallest,1)-(floor((SizeDiff-2*abs(RotateOrig - CentersDiff))/2)+StartEnd(Smallest,2)-size(Mu1Avg,1)) size(Mu1Avg,1)] ;
                end
            end
        end

        WinSizes = StartEnd(:,2) - StartEnd(:,1) + 1 ;

        maxRotation = 5 ;

        possibleRotates1 = max(1,StartEnd(1,1)-maxRotation):min(floor(size(Mu1Avg,1)-WinSizes(1)/2),min(StartEnd(1,1)+maxRotation,size(Mu1Avg,1)-WinSizes(1))) ;
        if ~isempty(possibleRotates1)
            DistDoubleRotate = zeros(length(possibleRotates1),1) ;

            Mu1RotateInd = 1 ;
            for RotateMu1 = possibleRotates1
                Mu1Rotated = Mu1Avg(RotateMu1:RotateMu1+WinSizes(1)-1,:) ;
                Mu2Rotated = Mu2Avg(StartEnd(2,1):StartEnd(2,2),:) ;
                DistDoubleRotate(Mu1RotateInd) = mean(sqrt(sum(((Mu1Rotated-Mu2Rotated)).^2,2)),"omitmissing") ;
                Mu1RotateInd = Mu1RotateInd + 1 ;
            end

            DistDoubleRotate = DistDoubleRotate - min(DistDoubleRotate,[],'all') ;
            minRotate1 = find(DistDoubleRotate == 0) ;
            StartEnd(1,:) = [possibleRotates1(minRotate1), possibleRotates1(minRotate1)+WinSizes(1)-1] ;
        end

        Window = tukeywin(WinSizes(1),0.25) ;
        indtukey = find(Window<0.01) ;

        Window = repmat(Window,1,size(Mu1Avg,2)) ;

        Mu1Avg = Mu1Avg(StartEnd(1,1):StartEnd(1,2),:) ;
        Mu1Avg = Mu1Avg.*Window ;
        Mu1Avg(indtukey,:) = 1e-4 ;

        Mu2Avg = Mu2Avg(StartEnd(2,1):StartEnd(2,2),:) ;
        Mu2Avg = Mu2Avg.*Window ;
        Mu2Avg(indtukey,:) = 1e-4 ;

        Mu1Avg_NoScale = Mu1Avg_NoScale(StartEnd(1,1):StartEnd(1,2),:) ;
        Mu1Avg_NoScale = Mu1Avg_NoScale.*Window ;
        Mu1Avg_NoScale(indtukey,:) = 1e-4 ;

        Mu2Avg_NoScale = Mu2Avg_NoScale(StartEnd(2,1):StartEnd(2,2),:) ;
        Mu2Avg_NoScale = Mu2Avg_NoScale.*Window ;
        Mu2Avg_NoScale(indtukey,:) = 1e-4 ;

        %% Compute Features (Described at the end of the code)

        AllDistRotate = zeros(size(Mu1Avg,1),1) ;

        for Rotate = 1:size(Mu1Avg,1)

            Mu2Rotate = repmat(Mu2Avg,2,1) ;
            Mu2Rotate = Mu2Rotate(Rotate:Rotate+size(Mu2Avg,1)-1,:) ;
            AllDistRotate(Rotate) = nanmean(sqrt(sum(((Mu1Avg-Mu2Rotate)).^2,2))) ;

        end

        AllDistRotateSmooth = zeros(size(Mu1Avg,1),1) ;

        for RotateSmooth = 1:size(Mu1Avg,1)

            Mu2RotateSmooth = repmat(Mu2Avg,2,1) ;
            Mu2RotateSmooth = Mu2RotateSmooth(RotateSmooth:RotateSmooth+size(Mu2Avg,1)-1,:) ;
            AllDistRotateSmooth(RotateSmooth) = nanmean(sqrt(sum(((Mu1Avg-Mu2RotateSmooth)).^2,2))) ;

        end

        [~,Rotate] = min(AllDistRotate) ;
        [~,RotateSmooth] = min(AllDistRotateSmooth) ;

        % Mu2RotateSmooth = repmat(Mu2Avg,2,1) ;
        % Mu2RotateSmooth = Mu2RotateSmooth(RotateSmooth:RotateSmooth+size(Mu2Avg,1)-1,:) ;
        % Mu2Rotate = repmat(Mu2Avg,2,1) ;
        % Mu2Rotate = Mu2Rotate(Rotate:Rotate+size(Mu2Avg,1)-1,:) ;

        Mu2RotateSmooth_NoScale = repmat(Mu2Avg_NoScale,2,1) ;
        Mu2RotateSmooth_NoScale = Mu2RotateSmooth_NoScale(RotateSmooth:RotateSmooth+size(Mu2Avg_NoScale,1)-1,:) ;
        Mu2Rotate_NoScale = repmat(Mu2Avg_NoScale,2,1) ;
        Mu2Rotate_NoScale = Mu2Rotate_NoScale(Rotate:Rotate+size(Mu2Avg_NoScale,1)-1,:) ;

        % calculate in 4D rather than 3D
        [Velocity2Rotate,~,~,~,~]  = lineParameters([Mu2RotateSmooth(size(Mu2RotateSmooth,1),:);Mu2RotateSmooth;Mu2RotateSmooth(1,:)],dt) ;
        [Velocity1,~,~,~,~]  = lineParameters([Mu1Avg(size(Mu1Avg,1),:);Mu1Avg;Mu1Avg(1,:)],dt) ;

        [~,Acceleration2Rotate,~,~,~]  = lineParameters([Mu2RotateSmooth(size(Mu2RotateSmooth,1),:);Mu2RotateSmooth;Mu2RotateSmooth(1,:)],dt) ;
        [~,Acceleration1,~,~,~]  = lineParameters([Mu1Avg(size(Mu1Avg,1),:);Mu1Avg;Mu1Avg(1,:)],dt) ;

        Velocity1 = Velocity1(:,2:end-1) ;
        Velocity2Rotate = Velocity2Rotate(:,2:end-1) ;

        Acceleration1 = Acceleration1(:,2:end-1) ;
        Acceleration2Rotate = Acceleration2Rotate(:,2:end-1) ;

        VeloTheta = acosd(dot(Velocity1,Velocity2Rotate,1)./(sqrt(sum(Velocity1.^2,1)).*sqrt(sum(Velocity2Rotate.^2,1)))) ;
        AccelTheta = acosd(dot(Acceleration1,Acceleration2Rotate,1)./(sqrt(sum(Acceleration1.^2,1)).*sqrt(sum(Acceleration2Rotate.^2,1)))) ;
        DistTheta = acosd(dot(Mu1Avg',Mu2RotateSmooth',1)./(sqrt(sum(Mu1Avg'.^2,1)).*sqrt(sum(Mu2RotateSmooth'.^2,1)))) ;

        % distDiffFinal=(sqrt(sum(((Mu1Avg-Mu2RotateSmooth)).^2,2)))';
        distDiffFinal=(sqrt(sum(((Mu1Avg_NoScale-Mu2RotateSmooth_NoScale)).^2,2)))';

        % [d,Mu2RotateSmoothZ,~] = procrustes(Mu1Avg,Mu2RotateSmooth,'reflection',false);
        [d,Mu2RotateSmoothZ_NoScale,~] = procrustes(Mu1Avg_NoScale,Mu2RotateSmooth_NoScale,'reflection',false);

        % distZDiffFinal = (sqrt(sum(((Mu1Avg-Mu2RotateSmoothZ)).^2,2)))';
        distZDiffFinal = (sqrt(sum(((Mu1Avg-Mu2RotateSmoothZ_NoScale)).^2,2)))';
        veloDiffFinal = sqrt(sum((Velocity1-Velocity2Rotate).^2,1));
        AccelDiffFinal = abs(sqrt(sum(Acceleration1.^2,1)) - sqrt(sum(Acceleration2Rotate.^2,1)));

        MeasureDistMinNew(Mu1Ind,Mu2Ind)=mean(distDiffFinal,'omitnan') ; %#ok<UDIM>
        ProSimilarity(Mu1Ind,Mu2Ind)=d;
        % DistSeuc(Mu1Ind,Mu2Ind)=mean(diag(pdist2(Mu1Avg,Mu2RotateSmoothZ,'seuclidean')),'omitnan');
        % DistCos(Mu1Ind,Mu2Ind)=mean(diag(pdist2(Mu1Avg(vecnorm(Mu1Avg')>=2.5*1e-4&vecnorm(Mu2RotateSmoothZ')>=2.5*1e-4,:),Mu2RotateSmoothZ(vecnorm(Mu1Avg')>=2.5*1e-4&vecnorm(Mu2RotateSmoothZ')>=2.5*1e-4,:),'cosine')),'omitnan');
        DistSeuc(Mu1Ind,Mu2Ind)=mean(diag(pdist2(Mu1Avg_NoScale,Mu2RotateSmoothZ_NoScale,'seuclidean')),'omitnan');
        DistCos(Mu1Ind,Mu2Ind)=mean(diag(pdist2(Mu1Avg_NoScale(vecnorm(Mu1Avg_NoScale')>=2.5*1e-4&vecnorm(Mu2RotateSmoothZ_NoScale')>=2.5*1e-4,:),Mu2RotateSmoothZ_NoScale(vecnorm(Mu1Avg_NoScale')>=2.5*1e-4&vecnorm(Mu2RotateSmoothZ_NoScale')>=2.5*1e-4,:),'cosine')),'omitnan');

        MeasureAccelMin(Mu1Ind,Mu2Ind)=mean(AccelDiffFinal,'omitnan') ;
        MeasureVeloMin(Mu1Ind,Mu2Ind)=mean(veloDiffFinal,'omitnan') ;

        AngleVeloMin(Mu1Ind,Mu2Ind) = mean(VeloTheta,'omitnan');
        AngleAccelMin(Mu1Ind,Mu2Ind) = mean(AccelTheta,'omitnan');
        AngleDistMin(Mu1Ind,Mu2Ind) = real(mean(DistTheta,'omitnan')) ;
        MeasureDistMinZtrans(Mu1Ind,Mu2Ind) = mean(distZDiffFinal,'omitnan'); %#ok<UDIM> 
        KeepMU_CVchange(Mu1Ind,Mu2Ind)=flagCV;
       
        %% Feature Calculations (Old features not used in paper)

        NewMu1 = [] ;
        NewMu2 = [] ;

        NewSmooth1 = [] ;
        NewSmooth2 = [] ;

        for scanChannels = 1:size(Mu1Avg,2)
            GrabThese = 1:size(Mu1Avg,2) ;
            GrabThese(scanChannels) = [] ;
            NewMu1 = cat(3,NewMu1,Mu1Avg(:,GrabThese)) ;
            NewMu2 = cat(3,NewMu2,Mu2RotateSmooth(:,GrabThese)) ;
        end

        Mu1AvgOld = Mu1Avg ;
        Mu2AvgOld = Mu2RotateSmooth ;
        
        for possibleChannels = 1:size(NewMu1,3)

            MaxAmps(Mu1Ind,Mu2Ind,possibleChannels) = 1./mean([max(Mu1AvgOld(:,possibleChannels))-min(Mu1AvgOld(:,possibleChannels)) max(Mu2AvgOld(:,possibleChannels))-min(Mu2AvgOld(:,possibleChannels))]) ;
            Mu1Avg = NewMu1(:,:,possibleChannels) ;
            Mu2Avg = NewMu2(:,:,possibleChannels) ;

            Mu2Rotate = Mu2Avg;
            MeasureDistMin(Mu1Ind,Mu2Ind,possibleChannels) = nanmean(sqrt(sum(((Mu1Avg-Mu2Rotate)).^2,2))) ;
            % MeasureDistMin(Mu1Ind,Mu2Ind,possibleChannels) = nanmean(sqrt(sum(((Mu1Avg_NoScale-Mu2Rotate_NoScale)).^2,2))) ;

            MeasureCurvatureMin(Mu1Ind,Mu2Ind,possibleChannels)  = zeros;
            MeasureVectorMin(Mu1Ind,Mu2Ind,possibleChannels)= zeros;
            MeasureCMCMax(Mu1Ind,Mu2Ind,possibleChannels) =zeros;

            AngleCurvMin(Mu1Ind,Mu2Ind,possibleChannels) =zeros;
            AngleVectorMin(Mu1Ind,Mu2Ind,possibleChannels) = zeros;
    
        end
    end
end

assignin('base','AllRatios',AllRatios) ;

%% Averaging over the different possible Channels (Test)

for scanX = 1:size(MaxAmps,2)
    for scanY = 1:size(MaxAmps,1)
        if scanY ~= scanX && sum(MaxAmps(scanY,scanX,:)) ~= 0
            MaxAmps(scanY,scanX,:) = MaxAmps(scanY,scanX,:)/sum(MaxAmps(scanY,scanX,:)) ;
        end
    end
end

MeasureDistMin = MeasureDistMin.*MaxAmps ;
MeasureCurvatureMin = MeasureCurvatureMin.*MaxAmps ;
MeasureVectorMin = MeasureVectorMin.*MaxAmps ;
MeasureCMCMax = MeasureCMCMax.*MaxAmps ;
AngleCurvMin = AngleCurvMin.*MaxAmps ;
AngleVectorMin = AngleVectorMin.*MaxAmps ;

MeasureDistMin =  mean(MeasureDistMin,3) ;
MeasureCurvatureMin = mean(MeasureCurvatureMin,3) ;
MeasureVectorMin = mean(MeasureVectorMin,3) ;
MeasureCMCMax = mean(MeasureCMCMax,3) ;
AngleCurvMin = mean(AngleCurvMin,3) ;
AngleVectorMin = mean(AngleVectorMin,3) ;

%% Test Scaling

cornerNaNs = zeros(size(numChans)) ;
cornerNaNs(1:MUsInSet1,1:MUsInSet1) = NaN ;
cornerNaNs(MUsInSet1+1:size(numChans,1),MUsInSet1+1:size(numChans,2)) = NaN ;

if simData
    for scanMat = 1:MUsInSet1
        numChans(scanMat,scanMat+MUsInSet1) = min(numChans(numChans ~= 0),[],'all') ;
        numChans(scanMat+MUsInSet1,scanMat) = min(numChans(numChans ~= 0),[],'all') ;
    end
end

numChans(numChans < mean(numChans(~isnan(cornerNaNs)),"all","omitmissing")) = mean(numChans(~isnan(cornerNaNs)),"all","omitmissing") ;
numChans = numChans./max(numChans,[],"all") ;

assignin('base','numChans',numChans) ;

%% Add NaNs to remove Zeros

MeasureDistMin = MeasureDistMin + cornerNaNs ;
MeasureVeloMin = MeasureVeloMin + cornerNaNs ;
MeasureAccelMin = MeasureAccelMin + cornerNaNs ;
MeasureCurvatureMin = MeasureCurvatureMin + cornerNaNs ;
MeasureVectorMin = MeasureVectorMin + cornerNaNs ;
MeasureCMCMax = MeasureCMCMax + cornerNaNs ;
AngleVeloMin = AngleVeloMin + cornerNaNs ;
AngleAccelMin = AngleAccelMin + cornerNaNs ;
AngleCurvMin = AngleCurvMin + cornerNaNs ;
AngleVectorMin = AngleVectorMin + cornerNaNs ;
MeasureDistMinNew = MeasureDistMinNew + cornerNaNs ;
AngleDistMin = AngleDistMin + cornerNaNs ;
MeasureDistMinZtrans = MeasureDistMinZtrans + cornerNaNs ;
ProSimilarity = ProSimilarity + cornerNaNs ;
DistSeuc = DistSeuc + cornerNaNs ;
DistCos = DistCos + cornerNaNs ;

DistUp = triu(MeasureDistMin./numChans) ;
DistDown = tril(MeasureDistMin./numChans) ;
MeasureDistMin = (DistUp+DistDown')/2 ;

DistUp2 = triu(MeasureDistMinNew./numChans) ;
DistDown2 = tril(MeasureDistMinNew./numChans) ;
MeasureDistMinNew = (DistUp2+DistDown2')/2 ;

DistZtransUp = triu(MeasureDistMinZtrans./numChans) ;
DistDistZtransDown = tril(MeasureDistMinZtrans./numChans) ;
MeasureDistMinZtrans = (DistZtransUp+DistDistZtransDown')/2 ;

AccelUp = triu(MeasureAccelMin./numChans) ;
AccelDown = tril(MeasureAccelMin./numChans) ;
MeasureAccelMin = (AccelUp+AccelDown')/2 ;

VeloUp = triu(MeasureVeloMin./numChans) ;
VeloDown = tril(MeasureVeloMin./numChans) ;
MeasureVeloMin = (VeloUp+VeloDown')/2 ;

CurvUp = triu(MeasureCurvatureMin./numChans) ;
CurvDown = tril(MeasureCurvatureMin./numChans) ;
MeasureCurvatureMin = (CurvUp+CurvDown')/2 ;

VectorUp = triu(MeasureVectorMin./numChans) ;
VectorDown = tril(MeasureVectorMin./numChans) ;
MeasureVectorMin = (VectorUp+VectorDown')/2 ;

CMCUp = triu(MeasureCMCMax./numChans) ;
CMCDown = tril(MeasureCMCMax./numChans) ;
MeasureCMCMax = (CMCUp+CMCDown')/2 ;

AngleVeloUp = triu(AngleVeloMin./numChans) ;
AngleVeloDown = tril(AngleVeloMin./numChans) ;
AngleVeloMin = (AngleVeloUp+AngleVeloDown')/2 ;

AngleAccelUp = triu(AngleAccelMin./numChans) ;
AngleAccelDown = tril(AngleAccelMin./numChans) ;
AngleAccelMin = (AngleAccelUp+AngleAccelDown')/2 ;

AngleDistUp = triu(AngleDistMin./numChans) ;
AngleDistDown = tril(AngleDistMin./numChans) ;
AngleDistMin = (AngleDistUp+AngleDistDown')/2 ;

AngleCurvUp = triu(AngleCurvMin./numChans) ;
AngleCurvDown = tril(AngleCurvMin./numChans) ;
AngleCurvMin = (AngleCurvUp+AngleCurvDown')/2 ;

AngleVectorUp = triu(AngleVectorMin./numChans) ;
AngleVectorDown = tril(AngleVectorMin./numChans) ;
AngleVectorMin = (AngleVectorUp+AngleVectorDown')/2 ;

ProSimilarityUp = triu(ProSimilarity./numChans) ;
ProSimilarityDown = tril(ProSimilarity./numChans) ;
ProSimilarity = (ProSimilarityUp+ProSimilarityDown')/2 ;

DistSeucUp = triu(DistSeuc./numChans) ;
DistSeucDown = tril(DistSeuc./numChans) ;
DistSeuc = (DistSeucUp+DistSeucDown')/2 ;

DistCosUp = triu(DistCos./numChans) ;
DistCosDown = tril(DistCos./numChans) ;
DistCos = (DistCosUp+DistCosDown')/2 ;

ValuesGrabber = tril(NaN(size(MeasureVeloMin))) ;
MeasureDistMin = MeasureDistMin + ValuesGrabber ;
MeasureDistMinZtrans = MeasureDistMinZtrans + ValuesGrabber ;
MeasureDistMinNew = MeasureDistMinNew + ValuesGrabber ;
MeasureVeloMin = MeasureVeloMin + ValuesGrabber ;
MeasureAccelMin = MeasureAccelMin + ValuesGrabber ;
MeasureCurvatureMin = MeasureCurvatureMin + ValuesGrabber ;
MeasureVectorMin = MeasureVectorMin + ValuesGrabber ;
MeasureCMCMax = MeasureCMCMax + ValuesGrabber ;
AngleVeloMin = AngleVeloMin + ValuesGrabber ;
AngleDistMin = AngleDistMin + ValuesGrabber ;
AngleCurvMin = AngleCurvMin + ValuesGrabber  ;
AngleVectorMin = AngleVectorMin + ValuesGrabber ;
AngleAccelMin = AngleAccelMin + ValuesGrabber ;
KeepMU_CVchange = KeepMU_CVchange + ValuesGrabber ;
DistCos = DistCos + ValuesGrabber ;
DistSeuc = DistSeuc + ValuesGrabber ;
ProSimilarity = ProSimilarity + ValuesGrabber ;

%% Code if using code on the same sets twice 

if sum(Txx_Smooth(:,:,1) == Txx_Smooth(:,:,MUsInSet1+1),"all") == size(Txx_Smooth,1)*size(Txx_Smooth,2) 
    for scanMat = 1:MUsInSet1
        MeasureDistMin(scanMat,scanMat+MUsInSet1) = 0 ;
        MeasureVeloMin(scanMat,scanMat+MUsInSet1) = 0 ;
        MeasureAccelMin(scanMat,scanMat+MUsInSet1) = 0 ;
        MeasureCurvatureMin(scanMat,scanMat+MUsInSet1) = 0 ;
        MeasureVectorMin(scanMat,scanMat+MUsInSet1) = 0 ;
        MeasureCMCMax(scanMat,scanMat+MUsInSet1) = 0 ;
        AngleVeloMin(scanMat,scanMat+MUsInSet1) = 0 ;
        AngleAccelMin(scanMat,scanMat+MUsInSet1) = 0 ;
        AngleCurvMin(scanMat,scanMat+MUsInSet1) = 0 ;
        AngleVectorMin(scanMat,scanMat+MUsInSet1) = 0 ;
        MeasureDistMinNew(scanMat,scanMat+MUsInSet1) = 0 ;
        AngleDistMin(scanMat,scanMat+MUsInSet1) = 0 ;
        MeasureDistMinZtrans(scanMat,scanMat+MUsInSet1) = 0 ;
        ProSimilarity(scanMat,scanMat+MUsInSet1) = 0 ;
        DistSeuc(scanMat,scanMat+MUsInSet1) = 0 ;
        DistCos(scanMat,scanMat+MUsInSet1) = 0 ;
    end
end

%% Removing unstable MUs and scaling cluster

v1Compare = MeasureDistMin(~isnan(MeasureDistMin)) ;
v2Compare = MeasureVeloMin(~isnan(MeasureVeloMin)) ;
v3Compare = MeasureAccelMin(~isnan(MeasureAccelMin)) ;
v4Compare = MeasureCurvatureMin(~isnan(MeasureCurvatureMin)) ;
v5Compare = MeasureVectorMin(~isnan(MeasureVectorMin)) ;
v6Compare = MeasureCMCMax(~isnan(MeasureCMCMax)) ;
v7Compare = AngleVeloMin(~isnan(AngleVeloMin)) ;
v8Compare = AngleAccelMin(~isnan(AngleAccelMin)) ;
v9Compare = AngleCurvMin(~isnan(AngleCurvMin)) ;
v10Compare = AngleVectorMin(~isnan(AngleVectorMin)) ;
v11Compare = MeasureDistMinNew(~isnan(MeasureDistMinNew)) ;
v12Compare = AngleDistMin(~isnan(AngleDistMin)) ;
v13Compare = MeasureDistMinZtrans(~isnan(MeasureDistMinZtrans)) ;
v14Compare = ProSimilarity(~isnan(ProSimilarity));
v15Compare = DistSeuc(~isnan(DistSeuc));
v16Compare = DistCos(~isnan(DistCos));

KeepMU_CVchange = KeepMU_CVchange(~isnan(KeepMU_CVchange)) ;
ClustCompare = [v1Compare v2Compare v3Compare v4Compare v5Compare v6Compare v7Compare v8Compare v9Compare v10Compare v11Compare v12Compare v13Compare v14Compare v15Compare v16Compare] ;
MU_Couple=zeros(size(v2Compare,1),2);
FindIndMatrix = v2Compare ;
SavedValues = [] ;

% ClustCompare is the output with features used in the paper:
    % 1 - Euclidean Distance - v11Compare
    % 2 - Standardized Euclidean Distance - v15Compare
    % 3 - Procrustes - v14Compare
    % 4 - Standardized Procruses - v13Compare
    % 5 - Velocity Magnitude - v2Compare
    % 6 - Acceleration Magnitude - v3Compare
    % 7 - Position Angle - v12Compare
    % 8 - Velocity Angle - v7Compare
    % 9 - Acceleration Angle - v8Compare

for indexOfPairs=1:size(v2Compare,1)
    if ~length(find(MeasureVeloMin == FindIndMatrix(indexOfPairs))) > 1
        [MU_Couple(indexOfPairs,1),MU_Couple(indexOfPairs,2)]=find(MeasureVeloMin == FindIndMatrix(indexOfPairs));
    else
        if isempty(find(SavedValues == FindIndMatrix(indexOfPairs)))
            [whoToTake1,whoToTake2] = find(MeasureVeloMin == FindIndMatrix(indexOfPairs)) ;
            MU_Couple(indexOfPairs,1) = whoToTake1(1) ;
            MU_Couple(indexOfPairs,2) = whoToTake2(1) ;
            SavedValues = [SavedValues;FindIndMatrix(indexOfPairs)] ;
        else
            [whoToTake1,whoToTake2] = find(MeasureVeloMin == FindIndMatrix(indexOfPairs)) ;
            MU_Couple(indexOfPairs,1) = whoToTake1(length(find(SavedValues == FindIndMatrix(indexOfPairs)))+1) ;
            MU_Couple(indexOfPairs,2) = whoToTake2(length(find(SavedValues == FindIndMatrix(indexOfPairs)))+1) ;
        end
    end
end

MU_Couple_Orig=GoodMus(MU_Couple);

if ~isempty(StabilityScore)

    ScaleFactor=ones(size(MU_Couple_Orig,1),size(ClustCompare,2));

    for ind=1:size(MU_Couple_Orig,1)
        range1 = [] ;
        range2 = [] ;

        if ~KeepMU_CVchange(ind,1)
            ScaleFactor(ind,1:size(StabilityScore{MU_Couple_Orig(ind,1),11},2))=prctile([real(StabilityScore{MU_Couple_Orig(ind,1),11}(1:end,:));real(StabilityScore{MU_Couple_Orig(ind,2),11}(1:end,:))],ptile)- mean([real(StabilityScore{MU_Couple_Orig(ind,1),12}(1:end,:));real(StabilityScore{MU_Couple_Orig(ind,2),12}(1:end,:))]);
        elseif KeepMU_CVchange(ind,1)==1
            ScaleFactor(ind,1:size(StabilityScore{MU_Couple_Orig(ind,1),11},2))=prctile([real(StabilityScore{MU_Couple_Orig(ind,1),11}(1:end,:))],ptile)- mean([real(StabilityScore{MU_Couple_Orig(ind,1),12}(1:end,:))]);
        else
            ScaleFactor(ind,1:size(StabilityScore{MU_Couple_Orig(ind,1),11},2))=prctile([real(StabilityScore{MU_Couple_Orig(ind,2),11}(1:end,:))],ptile)- mean([real(StabilityScore{MU_Couple_Orig(ind,2),12}(1:end,:))]);
        end
    end

else

    ScaleFactor = [] ;

end

end