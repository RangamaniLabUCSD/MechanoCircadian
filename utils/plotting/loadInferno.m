function infernoMap = loadInferno()
%LOADINFERNO load inferno colormap, assuming inferno.json is available in
%directory
% (this is tailored to a single json (inferno) currently,
% would need to be edited for loading in other color maps)
% [json_file,json_path] = uigetfile('*.json','Select json file with color map');
% json_file = fullfile(json_path,json_file);
json_file = 'inferno.json';
load_colors = readcell(json_file,'FileType','text');
RGBPoints = cell2mat(load_colors(16:1039));
RGBPoints = reshape(RGBPoints',[4, length(RGBPoints)/4]);
dataVals = RGBPoints(1,:);
dataInterp = dataVals(1):(dataVals(end)-dataVals(1))/1000:dataVals(end);
RGBVals = RGBPoints(2:4,:);
RInterp = interp1(dataVals,RGBVals(1,:),dataInterp);
GInterp = interp1(dataVals,RGBVals(2,:),dataInterp);
BInterp = interp1(dataVals,RGBVals(3,:),dataInterp);
RGBInterp = vertcat(RInterp,GInterp,BInterp);
infernoMap = RGBInterp(:,150:end)';
end

