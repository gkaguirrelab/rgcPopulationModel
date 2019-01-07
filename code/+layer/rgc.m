function rgc = rgc()
% Determined by 18 histologically sectioned maculas. Data from Curcio 2011,
% Figure 7B and 7C
%
%   Curcio, Christine A., et al. "Human chorioretinal layer thicknesses 
%   measured in macula-wide, high-resolution histologic sections." Investigative 
%   ophthalmology & visual science 52.7 (2011): 3943-3954.
%

%% RGC thickness

rgc.supportMM.temporal = [0	0.05	0.1	0.2	0.4	0.6	0.8	1	1.5	2	2.5	3];
rgc.thickMM.temporal = [1.741	6.278	11.139	24.376	53.589	58.073	59.321	56.696	49.066	36.2	29.051	19.608]./1000;
rgc.supportMM.nasal = [3.5	3	2.5	2	1.5	1	0.8	0.6	0.4	0.2	0.1	0.05	0];
rgc.thickMM.nasal = [13.763	18.447	24.214	36.521	54.329	72.797	71.259	68.945	52.759	18.282	8.541	2.624	1.741]./1000;

% Convert from mm to visual deg
rgc.supportDeg.temporal = convert_mmRetina_to_degRetina(rgc.supportMM.temporal);
rgc.supportDeg.nasal = convert_mmRetina_to_degRetina(rgc.supportMM.nasal);

% Obtain a spline fit to the thickness measurements
cardinalMeridianNames = {'temporal','nasal'};
for mm = 1:2
    tmpSupport = rgc.supportDeg.(cardinalMeridianNames{mm})';
    tmpVals = rgc.thickMM.(cardinalMeridianNames{mm})';
    rgc.fitDeg.(cardinalMeridianNames{mm}) = ...
        fit( tmpSupport, tmpVals, 'smoothingspline');
end

end

