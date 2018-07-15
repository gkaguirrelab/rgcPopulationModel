


function rgc = rgc( cardinalMeridianAngles, cardinalMeridianNames )


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

rgc.supportMM.temporal = [0, 0.09, 0.18, 0.26, 0.35, 0.4, 0.56, 0.69, 0.82, 0.91, 1.11, 1.26, 1.45, 1.7, 1.88, 2.06, 2.27, 2.48, 2.7, 2.94, 2.99];
rgc.thickMM.temporal = (1/(1-tissueShrinkage)).*[.00455, .00455, .01212, .02273, .03485, .04848, .05303, .05758, .06061, .06061, .05909, .05606, .05303, .05, .04545, .03636, .03485, .0303, .02727, .02121, .02121];
rgc.supportMM.nasal = [3, 2.88, 2.72, 2.54, 2.35, 2.15, 1.99, 1.86, 1.68, 1.53, 1.33, 1.14, 1.03, 0.87, 0.7, 0.6, 0.48, 0.31, 0.19, 0.11, 0.05, 0];
rgc.thickMM.nasal = (1/(1-tissueShrinkage)).*[.02121, .02121, .02273, .02576, .0303, .03485, .03788, .04394, .05, .05455, .07212, .06818, .07273, .07273, .07121, .0697, .05909, .05303, .03939, .0197, .00909, .00455];


% Obtain a spline fit to the thickness measurements
for mm = [1 3]
    tmpSupport = rgc.supportMM.(cardinalMeridianNames{mm})';
    tmpVals = rgc.thickMM.(cardinalMeridianNames{mm})';
    rgc.fitMM.(cardinalMeridianNames{mm}) = ...
        fit( tmpSupport, tmpVals, 'smoothingspline');
end

end

