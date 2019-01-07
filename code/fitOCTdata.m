function fitOCTdata(varargin )
% Fits the Curcio (2011) macular thickness measures to empirical OCT data
%
% Description:
%   Curcio and colleagues measured the histologic thickness of several
%   layers of the human retina in ex-vivo specimens. They report the
%   average thickness as a function of retinal eccentricity (mm) for each
%   of these layers. The report notes that these measurements may be
%   subject to some degree of tissue shrinkage.
%
%       
%
% Inputs:
%
% Optional key / value pairs:
%
% Outputs:
%
% Examples:
%


%% input parser
p = inputParser;

% Optional analysis params
p.addParameter('polarAngle',180,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal' 'superior' 'temporal' 'inferior'},@iscell);
p.addParameter('octDataFileName', ...
    fullfile(tbLocateProject('rgcPopulationModel','verbose',false),'data','rgcIplThicknessMap.mat'),...
    @ischar);

% parse
p.parse(varargin{:})


%% Obtain the TOME OCT thickness measurements
dataLoad = load(p.Results.octDataFileName);
rgcIplThicknessMap = dataLoad.rgcIplThicknessMap;
clear dataLoad

% Extract the horizontal meridian
rgciplOCTthickness = rgcIplThicknessMap(round(size(rgcIplThicknessMap,1)/2),:);

% Construct a data structure to hold the OCT thickness values and adjust
% the thickness measure to reflect just the RGC component
rgcOCTMicrons.supportDeg.temporal = ((1:round(length(rgciplOCTthickness)/2))./round(length(rgciplOCTthickness)/2)).*15;
rgcOCTMicrons.supportDeg.nasal = ((1:round(length(rgciplOCTthickness)/2))./round(length(rgciplOCTthickness)/2)).*15;
rgcOCTMicrons.thickMM.temporal = fliplr(rgciplOCTthickness(1:round(length(rgciplOCTthickness)/2))).*rgcIplThicknessRatio( rgcOCTMicrons.supportDeg.temporal )';
rgcOCTMicrons.thickMM.nasal = rgciplOCTthickness(round(length(rgciplOCTthickness)/2)+1:end).*rgcIplThicknessRatio( rgcOCTMicrons.supportDeg.nasal )';
% Obtain a spline fit to the thickness measurements
for mm = [1 3]
    tmpSupport = rgcOCTMicrons.supportDeg.(p.Results.cardinalMeridianNames{mm})';
    tmpVals = rgcOCTMicrons.thickMM.(p.Results.cardinalMeridianNames{mm})';
    nonNanIdx = ~isnan(tmpVals);
    rgcOCTMicrons.fitDeg.(p.Results.cardinalMeridianNames{mm}) = ...
        fit( tmpSupport(nonNanIdx), tmpVals(nonNanIdx), 'smoothingspline');
end


figure
supportDeg = 0:0.1:14;
subplot(1,2,1)
plot(supportDeg, rgcOCTMicrons.fitDeg.temporal(supportDeg))
xlim([0 15]);
ylim([0 75]);
supportDeg = 0:0.1:10;
subplot(1,2,2)
plot(supportDeg, rgcOCTMicrons.fitDeg.nasal(supportDeg))
xlim([0 15]);
ylim([0 75]);

end % fitOCTdata function

