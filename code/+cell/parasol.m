function parasol = parasol( cardinalMeridianAngles, cardinalMeridianNames, totalRGC, midget, bistratified )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



% %% Parasol RGCs
% % Average parasol cell densities across four meridians measured in
% % cells/square mm from six macaque retinas:
% % 
% %   Silveira, L. C. L., and V. H. Perry. "The topography of magnocellular
% %   projecting ganglion cells (M-ganglion cells) in the primate retina."
% %   Neuroscience 40.1 (1991): 217-237.
% 
% 
% % Parasol proportion of all retinal ganglion cells Data from Silviera et
% % al. 1991, Figure 17
% parasolMacaque.proportion.supportMM.nasal = [1.06, 2.07, 3.12, 4.12, 5.09, 6.17, 7.19, 8.22, 9.2, 10.24, 11.26, 12.24, 13.31, 14.32, 15.32, 16.42, 17.39];
% parasolMacaque.proportion.value.nasal = [.098, .0964, nan, nan, nan, .0758, .095, .1174, .1398, .1918, .2078, .1599, .1983, .2032, .2455, .2687, .1993];
% parasolMacaque.proportion.supportMM.temporal = [0 1.06, 2.07, 3.12, 4.12, 5.09, 6.17, 7.19, 8.22, 9.2, 10.24, 11.26, 12.24, 13.31, 14.32, 15.32, 16.42, 17.39];
% parasolMacaque.proportion.value.temporal = [.066, .0844, .0909, .0845, .0885, .1054, .1086, .0951, .0927, .0856, .0569, .0297, .0498, nan, nan, nan, nan];
% parasolMacaque.proportion.supportMM.superior = [.96, 1.99, 2.98, 3.93, 4.93, 5.99, 6.87, 7.94, 8.89, 9.88, 10.92, 11.91, 12.9, 13.89, 14.92];
% parasolMacaque.proportion.value.superior = [.058, .0471, .0448, .0568, .0791, .0808, .0857, .0938, .0923, .0853, .1044, .1148, .1332, .119, .0993];
% parasolMacaque.proportion.supportMM.inferior = [.96, 1.99, 2.98, 3.93, 4.93, 5.99, 6.87, 7.94, 8.89, 9.88, 10.92, 11.91, 12.9, 13.89, 14.92];
% parasolMacaque.proportion.value.inferior = [.0573, .0558, .0718, .056, .0696, .0872, .0881, .1001, .0828, .094, .0774, .0942, .0681, .0468, nan];
% 
% % Obtain a spline fit to the parasol densities
% for mm = 1:length(cardinalMeridianAngles)
%     nonNanSupportIdx = ~isnan(parasolMacaque.proportion.value.(cardinalMeridianNames{mm}));
%     tmpSupport = parasolMacaque.proportion.supportMM.(cardinalMeridianNames{mm})(nonNanSupportIdx)';
%     parasolMacaque.density.fitMMSq.(cardinalMeridianNames{mm}) = ...
%         fit( [0; tmpSupport], ...        
%         [0; totalRGC.density.fitMMSq.(cardinalMeridianNames{mm})(tmpSupport)' .* ...
%         parasolMacaque.proportion.value.(cardinalMeridianNames{mm})(nonNanSupportIdx)'], ...
%         'smoothingspline');
% end



%% Infer parasol densities
% We have solid measurements of total RGC density and a solid estimate of
% midget RGC density. The only parasol density measurements are from
% macaque. We assume that:
%   densityParasol = densityTotalRGC - densityMidget - densityBistratified
for mm = 1:length(cardinalMeridianAngles)
    parasol.density.fitDegSq.(cardinalMeridianNames{mm}) = @(x) totalRGC.density.fitDegSq.(cardinalMeridianNames{mm})(x')- ...
        midget.density.fitDegSq.(cardinalMeridianNames{mm})(x') - ...
        bistratified.density.fitDegSq.(cardinalMeridianNames{mm})(x');
end


% Parasol cell body sizes
% Data from Perry et al. 1984, Figure 6C:
%   Perry, V. H., R. Oehler, and A. Cowey. "Retinal ganglion cells that
%   project to the dorsal lateral geniculate nucleus in the macaque
%   monkey." Neuroscience 12.4 (1984): 1101-1123.
%
% These are the diameters observed as a function of eccentricity in the
% macauqe. Dacey 1993 reports that parasol cell bodies are 1.35 times
% larger in the human than in the macaque. We scale up the diameter values
% accordigly.
parasol.diameter.supportMM = [0.53, 1.07, 1.47, 1.96, 2.5, 3.1, 3.5, 4.0, 4.55, 5.0, 5.53, 6.0, 6.47, 7.0, 7.53, 8.05, 9.0, 9.65, 10.33, 11.4, 12.66, 13.63, 14.21];
parasol.diameter.supportDeg = convert_mmRetina_to_degVisual(parasol.diameter.supportMM, 180);
parasol.diameter.sizeMM = 1.35.*[0.0153, 0.017, 0.0186, 0.0203, 0.02187, 0.0196, 0.0246, 0.0252, 0.027, 0.0282, 0.0298, 0.0304, 0.0316, 0.0323, 0.0328, 0.0317, 0.037, 0.0351, 0.0327, 0.0321, 0.0306, 0.0302, 0.02713];
parasol.diameter.fitDeg = fit(parasol.diameter.supportDeg', parasol.diameter.sizeMM','poly1');


end

