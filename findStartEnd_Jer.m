function [startEnd] = findStartEnd_Jer(Txx,sampFreq)

for scanMUs = 1:size(Txx,3)

    if size(Txx,2) > 4
        rmsAllChans = rms(Txx(:,:,scanMUs),1) ;
        chansToUse = Txx(:,rmsAllChans > median(rmsAllChans),scanMUs) ;
    else
        chansToUse = Txx(:,:,scanMUs) ;
    end
    normMU = sqrt(sum(chansToUse.^2,2)) ;
    if ~(max(normMU) == 0)
        startBigPeak = find(normMU > 0.25*max(normMU),1,'first') ;
        endBigPeak = find(normMU > 0.25*max(normMU),1,'last') ;
        normMU(startBigPeak:endBigPeak) = max(normMU(startBigPeak:endBigPeak)) ;

        %% Method 1

        [bLow,aLow] = butter(4,250/(sampFreq/2),'low') ;
        normMUFilt = filtfilt(bLow,aLow,normMU) ;
        derivNorm = diff(normMUFilt) ;
        startMethod1 = find(derivNorm(1:startBigPeak)<0,1,'last') ;
        endMethod1 = endBigPeak + find(derivNorm(endBigPeak:end)<0,1,'first') ;

        %% Method 2

        medStart = median(normMU(1:startBigPeak)) ;
        startMethod2 = find(normMU(1:startBigPeak)<=medStart+0.01*max(normMU),1,'last') ;
        medEnd = median(normMU(endBigPeak:end)) ;
        endMethod2 = endBigPeak + find(normMU(endBigPeak:end)<=medEnd,1,'first') ;

        startEnd(scanMUs,1) = round(mean([startMethod1 startMethod2])) ;
        startEnd(scanMUs,2) = round(mean([endMethod1 endMethod2])) ;
    else
        startEnd(scanMUs,1) = 1 ;
        startEnd(scanMUs,2) = length(normMU) ;
    end
end

startEnd(startEnd > size(Txx,1)) = size(Txx,1) ;
startEnd(startEnd < 1) = 1 ;

end