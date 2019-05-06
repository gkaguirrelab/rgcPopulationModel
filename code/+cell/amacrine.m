function amacrine = amacrine( cardinalMeridianAngles, cardinalMeridianNames )

%% Displaced amacrine cells
% Curcio & Allen 1990 distinguished displaced amacrine cells from retinal
% ganglion cells through imaging and evalutating their morphology, and
% determined their soma size and densities at eccentricities across the
% human retina. Amacrine cell densities are averages across four meridians.
%
%   Curcio, Christine A., and Kimberly A. Allen. "Topography of ganglion
%   cells in human retina." Journal of comparative Neurology 300.1 (1990):
%   5-25.
%
% In Figure 3, Curcio & Allen note that displaced amacrine cells have soma
% diameters of 9-14 microns

for mm = 1:length(cardinalMeridianAngles)
    % Data from Curcio & Allen 1990, Figure 10, used for all meridians
    amacrine.density.supportMM.(cardinalMeridianNames{mm}) = [0, 0.11, 0.23, 0.4, 0.65, 0.86, 1.45, 2.46, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.6, 12.6, 13.6, 14.6, 15.6, 16.6, 17.7, 18.65, 19.65];
    amacrine.density.countsMMSq.(cardinalMeridianNames{mm}) = [73, 124, 317, 542, 689, 813, 1052, 1112, 1159, 1187, 1146, 1069, 992, 942, 888, 848, 798, 767, 731, 699, 677, 650, 642, 652, 676];

    % Convert from mm to visual degrees
    amacrine.density.supportDeg.(cardinalMeridianNames{mm}) = convert_mmRetina_to_degVisual(amacrine.density.supportMM.(cardinalMeridianNames{mm}),180);
    amacrine.density.countsDegSq.(cardinalMeridianNames{mm}) = amacrine.density.countsMMSq.(cardinalMeridianNames{mm}) .* calc_mmSqRetina_per_degSqVisual(amacrine.density.supportDeg.(cardinalMeridianNames{mm}), 180);
        
    % Obtain a spline fit to the amacrine densities
    amacrine.density.fitDegSq.(cardinalMeridianNames{mm}) = ...
        fit(amacrine.density.supportDeg.(cardinalMeridianNames{mm})', amacrine.density.countsDegSq.(cardinalMeridianNames{mm})', 'smoothingspline');
end

% Dacey measured the soma diameter in dopaminergic amacrine cells from 2-16
% mm eccentricity from the fovea in 112 M. nemestrina monkeys. The fit was
% extrapolated to 0mm eccentricity and out to just over 17mm. Soma diameter
% ranged from 11.2 to 16.6 microns. The data collected came from Dacey's
% best fit line which had a weakly significant positive correlation (p
% <.02).
% 
%   Dacey, Dennis M. "The dopaminergic amacrine cell." Journal of
%   Comparative Neurology 301.3 (1990): 461-489.
%
% Digitzed the line fit to the data in Figure 5B:
amacrine.diameter.supportMM = [0, 0.87, 1.68, 2.80, 3.76, 4.74, 5.86, 6.68, 7.36, 8.19, 8.71, 9.45, 9.97, 10.67, 11.45, 12.39, 13.47, 14.5, 15.62, 16.54, 17.82];
amacrine.diameter.supportDeg = convert_mmRetina_to_degVisual(amacrine.diameter.supportMM,180);
amacrine.diameter.sizeMM = [.01259, .01266, .01274, .01286, .01303, .0132, .01324, .01332, .01336, .01335, .01353, .01361, .01365, .01373, .01381, .01394, .01402, .01419, .01423, .01439, .01452];
amacrine.diameter.fitDeg = fit(amacrine.diameter.supportMM',amacrine.diameter.sizeMM','smoothingspline');

end

