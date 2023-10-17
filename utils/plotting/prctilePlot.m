function prctilePlot(xVals, yMat, varargin)
% plots shaded regions associated with interquartile range and range from
% 2.5 quantile to 97.5 quantile
% xVals is a 1 by numTests vector
% yMat is a matrix of data values, numTests by numSamples
    numTests = length(xVals);
    if size(yMat, 1) ~= numTests
        error('Rows of yMat should be the same as the length of numTests')
    end
    color = 'b';
    if length(varargin)==1
        color = varargin{1};
    elseif length(varargin) > 1
        error('Undefined input')
    end
    yPercentiles = prctile(yMat,[2.5,25,50,75,97.5],2);
    testsVecLoop = [xVals, xVals(end:-1:1), xVals(1)];
    upperRegion = [yPercentiles(:,5); yPercentiles(end:-1:1,3); yPercentiles(1,3)];
    upperQuartile = [yPercentiles(:,4); yPercentiles(end:-1:1,3); yPercentiles(1,4)];
    lowerQuartile = [yPercentiles(:,3); yPercentiles(end:-1:1,2); yPercentiles(1,3)];
    lowerRegion = [yPercentiles(:,3); yPercentiles(end:-1:1,1); yPercentiles(1,3)];
    plot(xVals, yPercentiles(:,3), color, 'LineWidth', 1)
    hold on
    fill(testsVecLoop, lowerRegion, color, 'FaceAlpha', .1, 'LineStyle', 'none')
    fill(testsVecLoop, lowerQuartile, color, 'FaceAlpha', .15, 'LineStyle', 'none')
    fill(testsVecLoop, upperRegion, color, 'FaceAlpha', .1, 'LineStyle', 'none')
    fill(testsVecLoop, upperQuartile, color, 'FaceAlpha', .15, 'LineStyle', 'none')

end
        