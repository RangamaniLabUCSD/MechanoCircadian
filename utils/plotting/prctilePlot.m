function prctilePlot(xVals, yMat, varargin)
% plots shaded regions associated with interquartile range and range from
% 2.5 quantile to 97.5 quantile
% xVals is a 1 by numTests vector
% yMat is a matrix of data values, numTests by numSamples
    numTests = length(xVals);
    if size(yMat, 1) ~= numTests
        error('Rows of yMat should be the same as the length of numTests')
    end
    % keepCol = true(1, size(yMat,2));
    % for i = 1:size(yMat,2)
    %     if any(isnan(yMat(:,i)))
    %         keepCol(i) = false;
    %     end
    % end
    % yMat = yMat(:,keepCol);
    color = 'b';
    edges = 2.5;
    if length(varargin)==1
        color = varargin{1};
    elseif length(varargin)==2
        color = varargin{1};
        edges = varargin{2};
    elseif length(varargin) > 2
        error('Undefined input')
    end
    yPercentiles = zeros(numTests,5);
    for i = 1:size(yMat,1)
        yCur = yMat(i,:);
        yPercentiles(i,:) = prctile(yCur(~isnan(yCur)),[edges,25,50,75,100-edges]);
    end
    yPercentiles = prctile(yMat,[edges,25,50,75,100-edges],2);
    testsVecLoop = [xVals, xVals(end:-1:1), xVals(1)];
    upperRegion = [yPercentiles(:,5); yPercentiles(end:-1:1,3); yPercentiles(1,3)];
    upperQuartile = [yPercentiles(:,4); yPercentiles(end:-1:1,3); yPercentiles(1,4)];
    lowerQuartile = [yPercentiles(:,3); yPercentiles(end:-1:1,2); yPercentiles(1,3)];
    lowerRegion = [yPercentiles(:,3); yPercentiles(end:-1:1,1); yPercentiles(1,3)];
    plot(xVals, yPercentiles(:,3), 'Color', color, 'LineWidth', 1, 'LineStyle', '-', 'Marker', 'none')
    hold on
    fill(testsVecLoop, lowerRegion, color, 'FaceAlpha', .1, 'LineStyle', 'none','Marker','none')
    fill(testsVecLoop, lowerQuartile, color, 'FaceAlpha', .15, 'LineStyle', 'none','Marker','none')
    fill(testsVecLoop, upperRegion, color, 'FaceAlpha', .1, 'LineStyle', 'none','Marker','none')
    fill(testsVecLoop, upperQuartile, color, 'FaceAlpha', .15, 'LineStyle', 'none','Marker','none')
end
        