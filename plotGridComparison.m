function plotGridComparison(emgGrid1,emgGrid2,norm,hm,maxSize)
%plotGrid Function to plot the EMG and motor units from a grid array
%
% INPUTS:
%           emgGrid         Data from a grid electrode
%           norm            Normalisation method - 'Grid', each channel is
%                           normalised to the entire grid, 'Channel', each
%                           channel is normalised to itself.
%           fs              Sampling Rate
%           hm              Heatmap
%
% Ben O'Callaghan 2020

switch nargin
    case 2
        norm = 'Grid';
        hm = 0;
        maxSize = 0 ;
    case 3
        hm = 0;
        maxSize = 0 ;
    case 4
        maxSize = 0 ;
end

if length(size(emgGrid1)) == 2
    emgGrid1 = permute(reshape(emgGrid1.',9,14,[]),[2 1 3]); % FDI Grid in the correct order if specified as 2D input (BioSemi Flex-print).
end

if length(size(emgGrid2)) == 2
    emgGrid2 = permute(reshape(emgGrid2.',9,14,[]),[2 1 3]); % FDI Grid in the correct order if specified as 2D input (BioSemi Flex-print).
end

if hm
    heatmap = true;
    hMap1 = rms(emgGrid1,3);
    hMap2 = rms(emgGrid2,3);
else
    heatmap = false;
end

nRows = size(emgGrid1,1);
nCols = size(emgGrid1,2);

x1 = linspace(-0.4,0.4,size(emgGrid1,3));
x2 = linspace(-0.4,0.4,size(emgGrid2,3));

if maxSize
    figure('units','normalized','outerposition',[0 0 1 1])
else
    figure
end
subplot(121)
if heatmap
    imagesc(hMap1)
end
hold on

for row = 1:nRows
    for col = 1:nCols
        if strcmp(norm,'Grid')
            plot(x1+col,0.5*squeeze(emgGrid1(row,col,:))./max(abs(emgGrid1(:)))+row,'k')
        elseif strcmp(norm,'Channel')
            plot(x1+col,0.5*squeeze(emgGrid1(row,col,:))./max(abs(emgGrid1(row,col,:)))+row,'k')
        else
            error('Normalisation should be of either ''Grid'' or ''Channel''');
        end
    end
end
xlabel('Column')
ylabel('Row')
xticks(1:nCols)
yticks(1:nRows)
set(gca,'YDir','reverse')
axis square

subplot(122)
if heatmap
    imagesc(hMap2)
end
hold on

for row = 1:nRows
    for col = 1:nCols
        if strcmp(norm,'Grid')
            plot(x2+col,0.5*squeeze(emgGrid2(row,col,:))./max(abs(emgGrid2(:)))+row,'k')
        elseif strcmp(norm,'Channel')
            plot(x2+col,0.5*squeeze(emgGrid2(row,col,:))./max(abs(emgGrid2(row,col,:)))+row,'k')
        else
            error('Normalisation should be of either ''Grid'' or ''Channel''');
        end
    end
end
xlabel('Column')
ylabel('Row')
xticks(1:nCols)
yticks(1:nRows)
set(gca,'YDir','reverse')
axis square
end

