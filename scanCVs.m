function [cvRatio,whoToResample] = scanCVs(Mu1,Mu2)

Mu1 = movmean(Mu1,3) ;
Mu2 = movmean(Mu2,3) ;

Mu1 = resample(Mu1,1000,1) ;
Mu2 = resample(Mu2,1000,1) ;

NormeMu1 = sqrt(sum(Mu1.^2,2)) ;
[~,maxPos1] = max(NormeMu1) ;

NormeMu2 = sqrt(sum(Mu2.^2,2)) ;
[~,maxPos2] = max(NormeMu2) ;

otherPointDist1 = pdist2(Mu1(maxPos1,:),Mu1) ;
[~,otherPointPos1] =  max(otherPointDist1) ;

otherPointDist2 = pdist2(Mu2(maxPos2,:),Mu2) ;
[~,otherPointPos2] =  max(otherPointDist2) ;

mu1Set1 = round(linspace(min(otherPointPos1,maxPos1),max(otherPointPos1,maxPos1),11))' ;
distMu1Set1ToMu2 = pdist2(Mu1(mu1Set1,:),Mu2) ;
[~,mu1Set2] = min(distMu1Set1ToMu2,[],2) ;

mu2Set1 = round(linspace(min(otherPointPos2,maxPos2),max(otherPointPos2,maxPos2),11))' ;
distMu2Set1ToMu1 = pdist2(Mu2(mu2Set1,:),Mu1) ;
[~,mu2Set2] = min(distMu2Set1ToMu1,[],2) ;

segLenghtsMu1Set1 = abs(mu1Set1(2:end) - mu1Set1(1:end-1)) ;
segLenghtsMu1Set2 = abs(mu1Set2(2:end) - mu1Set2(1:end-1)) ;
segLenghtsMu1Set2 = segLenghtsMu1Set2 + 1*(segLenghtsMu1Set2 == 0) ;

segLenghtsMu2Set1 = abs(mu2Set1(2:end) - mu2Set1(1:end-1)) ;
segLenghtsMu2Set2 = abs(mu2Set2(2:end) - mu2Set2(1:end-1)) ;
segLenghtsMu2Set2 = segLenghtsMu2Set2 + 1*(segLenghtsMu2Set2 == 0) ;

ratioMu1 = median(segLenghtsMu1Set1./segLenghtsMu1Set2) ;
ratioMu2 = median(segLenghtsMu2Set1./segLenghtsMu2Set2) ;

if ratioMu1 < 1
    ratioMu1 = 1./ratioMu1 ;
    whoToResample = 1 ;
end
if ratioMu2 < 1
    ratioMu2 = 1./ratioMu2 ;
    whoToResample = 2 ;
end

if  ~exist('whoToResample','var')
    whoToResample = 0 ;
    cvRatio = 10 ;
else
    cvRatio = (ratioMu1+ratioMu2)./2 ;
end

end

