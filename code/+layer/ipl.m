function ipl = ipl( cardinalMeridianAngles, cardinalMeridianNames )


%% RGC and IPL thickness
% Determined by 18 histologically sectioned maculas. Data from Curcio 2011,
% Figure 7B and 7C
%
%   Curcio, Christine A., et al. "Human chorioretinal layer thicknesses 
%   measured in macula-wide, high-resolution histologic sections." Investigative 
%   ophthalmology & visual science 52.7 (2011): 3943-3954.
%


% Tissue shrinkage. The Curcio RGC/IPL layer thickness measurements were
% subjected to a 29% reduction in size due to histological prep
tissueShrinkage = 0.29;


ipl.supportMM.temporal = [0, 0.1, 0.23, 0.38, 0.59, 0.79, 1.04, 1.31, 1.48, 1.76, 1.98, 2.19, 2.45, 2.75, 3];
ipl.thickMM.temporal = (1/(1-tissueShrinkage)).*[.00303, .00758, .01667, .02273, .02879, .03333, .03636, .03939, .03939, .04091, .04091, .03939, .03788, .03636, .03333];
ipl.supportMM.nasal = [3, 2.86, 2.69, 2.55, 2.34, 2.13, 1.92, 1.69, 1.56, 1.34, 1.15, 0.84, 0.55, 0.37, 0.13, 0.06, 0];
ipl.thickMM.nasal = (1/(1-tissueShrinkage)).*[.0303, .03182, .03182, .03182, .03333, .03485, .03636, .03788, .03939, .03939, .03788, .03485, .02879, .02121, .00909, .00455, .00303];

% Tissue shrinkage. The Curcio RGC/IPL layer thickness measurements were
% subjected to a 29% reduction in size due to histological prep
tissueShrinkage = 0.29;

% Obtain a spline fit to the thickness measurements
for mm = [1 3]
    tmpSupport = ipl.supportMM.(cardinalMeridianNames{mm})';
    tmpVals = ipl.thickMM.(cardinalMeridianNames{mm})';
    ipl.fitMM.(cardinalMeridianNames{mm}) = ...
        fit( tmpSupport, tmpVals, 'smoothingspline');
end

end

