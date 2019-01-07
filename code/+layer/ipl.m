function ipl = ipl()
% Determined by 18 histologically sectioned maculas. Data from Curcio 2011,
% Figure 7B and 7C
%
%   Curcio, Christine A., et al. "Human chorioretinal layer thicknesses 
%   measured in macula-wide, high-resolution histologic sections." Investigative 
%   ophthalmology & visual science 52.7 (2011): 3943-3954.
%

%% IPL thickness

ipl.supportMM.temporal = [0	0.05	0.1	0.2	0.4	0.6	0.8	1	1.5	2	2.5	3];
ipl.thickMM.temporal = [0.919	2.841	5.922	14.218	21.359	27.068	31.653	35.258	38.958	39.224	36.958	32.184]./1000;
ipl.supportMM.nasal = [3.5	3	2.5	2	1.5	1	0.8	0.6	0.4	0.2	0.1	0.05	0];
ipl.thickMM.nasal = [19.97	29.002	31.349	34.267	38.472	36.911	33.145	28.573	22.142	11.219	5.934	2.702	0.919]./1000;

% Convert mm to deg visual
ipl.supportDeg.temporal = convert_mmRetina_to_degRetina(ipl.supportMM.temporal);
ipl.supportDeg.nasal = convert_mmRetina_to_degRetina(ipl.supportMM.nasal);

% Obtain a spline fit to the thickness measurements
cardinalMeridianNames = {'temporal','nasal'};
for mm = 1:2
    tmpSupport = ipl.supportDeg.(cardinalMeridianNames{mm})';
    tmpVals = ipl.thickMM.(cardinalMeridianNames{mm})';
    ipl.fitDeg.(cardinalMeridianNames{mm}) = ...
        fit( tmpSupport, tmpVals, 'smoothingspline');
end

end

