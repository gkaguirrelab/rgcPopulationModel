function rgcThickness(varargin )
% Caclulates RGC layer thickness by reference to constituent cell classes
%
% Description:
%   We aim to create a model of the thickness of the retinal ganglion cell
%   layer and how it vaires with eccentricity, comparing our estimations to
%   reports from histology and OCT scanned images. We will calculate the
%   thickness of the RGC layer using data from histology on the population
%   of midget, parasol, and bistratified ganglion cells as well as
%   displaced amacrine cells. Specifically, the proportion of cell types,
%   their size, and their densities and how these three amounts vary with
%   eccentricity. The volume of space each cell type takes up will be
%   combined to create a coherent digital map of the human retinal ganglion
%   cell layer. This project utilizes largely generalized averages of cell
%   sizes and proportions, and does not take into account the different
%   subtypes of ganglion cells such as ON and OFF center midget and parasol
%   cells, widefield, smallfield, and diffuse bistratified cells, or the
%   many types of amacrine cells. The use of average proportions, sizes,
%   and densities, are noted throughout and still create an accurate
%   estimate of RGC thickness.
%
% Inputs:
%
% Optional key / value pairs:
%
% Outputs:
%
% Examples:
%

% Things for GKA to check out:
%   - See what x-axis scale factor would be required to bring our
%   thickness model into alignment with the Curcio 2011 data
%   - See if adding a higher packing density in the periphery fixes the
%   plateua


%% input parser
p = inputParser;

% Optional analysis params
p.addParameter('polarAngle',180,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal' 'superior' 'temporal' 'inferior'},@iscell);



% parse
p.parse(varargin{:})


%% Total RGC density function
% Create a function to return total RGC densities per square mm of retina
% at a given eccentricity position in mm retina.
% To do so, we call the rgcDisplacementMap functions which have a
% representation of the Curcio & Allen 1990 RGC density results. The
% rgcDisplacementMap toolbox is referenced in eccentricity units of degrees
% retina, and provides densities in square degrees retina. We convert those
% values here to mm and square mm.
for mm = 1:length(p.Results.cardinalMeridianAngles)
    tmpFit = getSplineFitToRGCDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
    totalRGC.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = @(posMMretina) tmpFit(convert_mmRetina_to_degRetina(posMMretina))'.*calc_degSqRetina_per_mmSqRetina();
end


%% Displaced amacrine cells
% Curcio & Allen 1990 distinguished displaced amacrine cells from retinal
% ganglion cells through imaging and evalutating their morphology, and
% determined their soma size and densities at eccentricities across the
% human retina. Amacrine cell densities are averages across four meridians.
%
%   Curcio, Christine A., and Kimberly A. Allen. "Topography of ganglion
%   cells in human retina." Journal of comparative Neurology 300.1 (1990):
%   5-25.

for mm = 1:length(p.Results.cardinalMeridianAngles)
    % Data from Curcio & Allen 1990, Figure 10, used for all meridians
    amacrine.density.supportMM.(p.Results.cardinalMeridianNames{mm}) = [0, 0.11, 0.23, 0.4, 0.65, 0.86, 1.45, 2.46, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.6, 12.6, 13.6, 14.6, 15.6, 16.6, 17.7, 18.65, 19.65];
    amacrine.density.countsMMSq.(p.Results.cardinalMeridianNames{mm}) = [73, 124, 317, 542, 689, 813, 1052, 1112, 1159, 1187, 1146, 1069, 992, 942, 888, 848, 798, 767, 731, 699, 677, 650, 642, 652, 676];
    % Obtain a spline fit to the amacrine densities
    amacrine.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = ...
        fit(amacrine.density.supportMM.(p.Results.cardinalMeridianNames{mm})', amacrine.density.countsMMSq.(p.Results.cardinalMeridianNames{mm})', 'smoothingspline');
end

% Dacey measured the soma diameter in dopaminergic amacrine cells from 2-16
% mm eccentricity from the fovea in 112 M. nemestrina monkeys. The fit was
% extrapolated to 0mm eccentricity and out to just over 17mm. Soma diameter
% ranged from 11.2 to 16.6 microns. The data collected came from Dacey's
% best fit line which had a weak positive significant correlation (p <.02).
% 
%   Dacey, Dennis M. "The dopaminergic amacrine cell." Journal of
%   Comparative Neurology 301.3 (1990): 461-489.
%
% Digitzed the line fit to the data in Figure:
amacrine.diameter.supportMM = [0, 0.87, 1.68, 2.80, 3.76, 4.74, 5.86, 6.68, 7.36, 8.19, 8.71, 9.45, 9.97, 10.67, 11.45, 12.39, 13.47, 14.5, 15.62, 16.54, 17.82];
amacrine.diameter.sizeMM = [.01259, .01266, .01274, .01286, .01303, .0132, .01324, .01332, .01336, .01335, .01353, .01361, .01365, .01373, .01381, .01394, .01402, .01419, .01423, .01439, .01452];
amacrine.diameter.fitMM = fit(amacrine.diameter.supportMM',amacrine.diameter.sizeMM','smoothingspline');


%% Parasol RGCs
% Average parasol cell densities across four meridians measured in
% cells/square mm from six macaque retinas:
% 
%   Silveira, L. C. L., and V. H. Perry. "The topography of magnocellular
%   projecting ganglion cells (M-ganglion cells) in the primate retina."
%   Neuroscience 40.1 (1991): 217-237.


% Parasol proportion of all retinal ganglion cells Data from Silviera et
% al. 1991, Figure 17
parasolMacaque.proportion.supportMM.nasal = [1.06, 2.07, 3.12, 4.12, 5.09, 6.17, 7.19, 8.22, 9.2, 10.24, 11.26, 12.24, 13.31, 14.32, 15.32, 16.42, 17.39];
parasolMacaque.proportion.value.nasal = [.098, .0964, nan, nan, nan, .0758, .095, .1174, .1398, .1918, .2078, .1599, .1983, .2032, .2455, .2687, .1993];
parasolMacaque.proportion.supportMM.temporal = [0 1.06, 2.07, 3.12, 4.12, 5.09, 6.17, 7.19, 8.22, 9.2, 10.24, 11.26, 12.24, 13.31, 14.32, 15.32, 16.42, 17.39];
parasolMacaque.proportion.value.temporal = [.066, .0844, .0909, .0845, .0885, .1054, .1086, .0951, .0927, .0856, .0569, .0297, .0498, nan, nan, nan, nan];
parasolMacaque.proportion.supportMM.superior = [.96, 1.99, 2.98, 3.93, 4.93, 5.99, 6.87, 7.94, 8.89, 9.88, 10.92, 11.91, 12.9, 13.89, 14.92];
parasolMacaque.proportion.value.superior = [.058, .0471, .0448, .0568, .0791, .0808, .0857, .0938, .0923, .0853, .1044, .1148, .1332, .119, .0993];
parasolMacaque.proportion.supportMM.inferior = [.96, 1.99, 2.98, 3.93, 4.93, 5.99, 6.87, 7.94, 8.89, 9.88, 10.92, 11.91, 12.9, 13.89, 14.92];
parasolMacaque.proportion.value.inferior = [.0573, .0558, .0718, .056, .0696, .0872, .0881, .1001, .0828, .094, .0774, .0942, .0681, .0468, nan];

% Obtain a spline fit to the parasol densities
for mm = 1:length(p.Results.cardinalMeridianAngles)
    nonNanSupportIdx = ~isnan(parasolMacaque.proportion.value.(p.Results.cardinalMeridianNames{mm}));
    tmpSupport = parasolMacaque.proportion.supportMM.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)';
    parasolMacaque.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = ...
        fit( [0; tmpSupport], ...        
        [0; totalRGC.density.fitMMSq.(p.Results.cardinalMeridianNames{mm})(tmpSupport)' .* ...
        parasolMacaque.proportion.value.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)'], ...
        'smoothingspline');
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
parasol.diameter.sizeMM = 1.35.*[0.0153, 0.017, 0.0186, 0.0203, 0.02187, 0.0196, 0.0246, 0.0252, 0.027, 0.0282, 0.0298, 0.0304, 0.0316, 0.0323, 0.0328, 0.0317, 0.037, 0.0351, 0.0327, 0.0321, 0.0306, 0.0302, 0.02713];
parasol.diameter.fitMM = fit(parasol.diameter.supportMM', parasol.diameter.sizeMM','poly1');


%% Bistratified RGCs

% Bistratified proportion
% Data from Dacey 1993, Figure 13b:
%   Dacey, Dennis M. "Morphology of a small-field bistratified ganglion
%   cell type in the macaque and human retina." Visual neuroscience 10.6
%   (1993): 1081-1098.
%
% Dacey provided values for the temporal retina; we assume these values
% hold for all meridians
for mm = 1:length(p.Results.cardinalMeridianAngles)
    bistratified.proportion.supportMM.(p.Results.cardinalMeridianNames{mm}) = [.97, 1.96, 2.91, 3.92, 4.91, 5.92, 6.89, 7.84, 8.85, 9.86, 10.89, 11.9, 12.91, 13.84, 14.91];
    bistratified.proportion.value.(p.Results.cardinalMeridianNames{mm}) = [.0135, .0168, .0202, .0241, .0284, .0324, .0364, .0403, .0447, .0485, .0538, .0573, .0603, .0641, .0662];
end

% Obtain a spline fit to the bistratified densities
for mm = 1:length(p.Results.cardinalMeridianAngles)
    nonNanSupportIdx = ~isnan(bistratified.proportion.value.(p.Results.cardinalMeridianNames{mm}));
    tmpSupport = bistratified.proportion.supportMM.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)';
    bistratified.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = ...
        fit( [0; tmpSupport], ...        
        [0; totalRGC.density.fitMMSq.(p.Results.cardinalMeridianNames{mm})(tmpSupport)' .* ...
        bistratified.proportion.value.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)'], ...
        'smoothingspline');
end


% Bistratified cell body sizes
% Data from Figure 3B:
%   Peterson, Beth B., and Dennis M. Dacey. "Morphology of wide-field
%   bistratified and diffuse human retinal ganglion cells." Visual
%   neuroscience 17.4 (2000): 567-578.
%
% and Table 1 of:
%
%   Dacey, Dennis M. "Morphology of a small-field bistratified ganglion
%   cell type in the macaque and human retina." Visual neuroscience 10.6
%   (1993): 1081-1098.
%
% We model the cell body diameter as constant as a function of
% eccentricity.
bistratified.diameter.fitMM = @(x) repmat(0.0189,size(x));


%% Midget RGCs

% Obtain a spline fit to RGC density values from Curcio & Allen
for mm = 1:length(p.Results.cardinalMeridianAngles)
    midget.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = createMidgetDensityFunc(p.Results.cardinalMeridianAngles(mm));
end

% Midget cell body sizes Liu and colleagues 2017 report midget soma sizes
% as a function of eccentricity as measured by AO-OCT:
%   Liu, Zhuolin, et al. "Imaging and quantifying ganglion cells and other
%   transparent neurons in the living human retina." Proceedings of the
%   National Academy of Sciences (2017): 201711734.
%
% Dacey 1993 (Table 1) reports a larget size of 17.9 microns.
%
%   Dacey, Dennis M. "Morphology of a small-field bistratified ganglion
%   cell type in the macaque and human retina." Visual neuroscience 10.6
%   (1993): 1081-1098.

% Support in the source data is in degrees of visual field along the
% temporal retina
midget.diameter.supportDeg = [(1.5+3)/2, (3+4.5)/2, (6+7.5)/2, (8+9.5)/2, (12+13.5)/2];

% We convert from visual degrees back to retinal mm.
midget.diameter.supportMM = convert_degVisual_to_mmRetina(midget.diameter.supportDeg, 180);
midget.diameter.sizeMM = [0.0115, 0.0113, 0.0114, 0.0118, 0.01315];

% Obtain a spline fit
midget.diameter.fitMM = fit(midget.diameter.supportMM', midget.diameter.sizeMM', 'poly1');


%% Infer parasol densities
% We have solid measurements of total RGC density and a solid estimate of
% midget RGC density. The only parasol density measurements are from
% macaque. We assume that:
%   densityParasol = densityTotalRGC - densityMidget - densityBistratified
for mm = 1:length(p.Results.cardinalMeridianAngles)
    parasol.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = @(x) totalRGC.density.fitMMSq.(p.Results.cardinalMeridianNames{mm})(x)- ...
        midget.density.fitMMSq.(p.Results.cardinalMeridianNames{mm})(x) - ...
        bistratified.density.fitMMSq.(p.Results.cardinalMeridianNames{mm})(x);
end


%% RGC and IPL thickness
% Determined by 18 histologically sectioned maculas. Data from Curcio 2011,
% Figure 7B and 7C
%
%   Curcio, Christine A., et al. "Human chorioretinal layer thicknesses 
%   measured in macula-wide, high-resolution histologic sections." Investigative 
%   ophthalmology & visual science 52.7 (2011): 3943-3954.
%
empiricalThickness.RGC.supportMM.temporal = [0, 0.09, 0.18, 0.26, 0.35, 0.4, 0.56, 0.69, 0.82, 0.91, 1.11, 1.26, 1.45, 1.7, 1.88, 2.06, 2.27, 2.48, 2.7, 2.94, 2.99];
empiricalThickness.RGC.thickMM.temporal = [.00455, .00455, .01212, .02273, .03485, .04848, .05303, .05758, .06061, .06061, .05909, .05606, .05303, .05, .04545, .03636, .03485, .0303, .02727, .02121, .02121];
empiricalThickness.RGC.supportMM.nasal = [3, 2.88, 2.72, 2.54, 2.35, 2.15, 1.99, 1.86, 1.68, 1.53, 1.33, 1.14, 1.03, 0.87, 0.7, 0.6, 0.48, 0.31, 0.19, 0.11, 0.05, 0];
empiricalThickness.RGC.thickMM.nasal = [.02121, .02121, .02273, .02576, .0303, .03485, .03788, .04394, .05, .05455, .07212, .06818, .07273, .07273, .07121, .0697, .05909, .05303, .03939, .0197, .00909, .00455];
empiricalThickness.IPL.supportMM.temporal = [0, 0.1, 0.23, 0.38, 0.59, 0.79, 1.04, 1.31, 1.48, 1.76, 1.98, 2.19, 2.45, 2.75, 3];
empiricalThickness.IPL.thickMM.temporal = [.00303, .00758, .01667, .02273, .02879, .03333, .03636, .03939, .03939, .04091, .04091, .03939, .03788, .03636, .03333];
empiricalThickness.IPL.supportMM.nasal = [3, 2.86, 2.69, 2.55, 2.34, 2.13, 1.92, 1.69, 1.56, 1.34, 1.15, 0.84, 0.55, 0.37, 0.13, 0.06, 0];
empiricalThickness.IPL.thickMM.nasal = [.0303, .03182, .03182, .03182, .03333, .03485, .03636, .03788, .03939, .03939, .03788, .03485, .02879, .02121, .00909, .00455, .00303];

% Obtain a spline fit to the thickness measurements
for mm = [1 3]
    tmpSupport = empiricalThickness.RGC.supportMM.(p.Results.cardinalMeridianNames{mm})';
    tmpVals = empiricalThickness.RGC.thickMM.(p.Results.cardinalMeridianNames{mm})';
    empiricalThickness.RGC.fitMM.(p.Results.cardinalMeridianNames{mm}) = ...
        fit( tmpSupport, tmpVals, 'smoothingspline');
    tmpSupport = empiricalThickness.IPL.supportMM.(p.Results.cardinalMeridianNames{mm})';
    tmpVals = empiricalThickness.IPL.thickMM.(p.Results.cardinalMeridianNames{mm})';
    empiricalThickness.IPL.fitMM.(p.Results.cardinalMeridianNames{mm}) = ...
        fit( tmpSupport, tmpVals, 'smoothingspline');
end



%% Figure prep
figure
supportMM = 0:0.01:6;

% Plot counts
subplot(1,3,1)
plot(supportMM, totalRGC.density.fitMMSq.temporal(supportMM))
hold on
plot(supportMM, midget.density.fitMMSq.temporal(supportMM))
plot(supportMM, parasol.density.fitMMSq.temporal(supportMM))
plot(supportMM, bistratified.density.fitMMSq.temporal(supportMM))
plot(supportMM, amacrine.density.fitMMSq.temporal(supportMM))
plot(supportMM, midget.density.fitMMSq.temporal(supportMM) + parasol.density.fitMMSq.temporal(supportMM) + bistratified.density.fitMMSq.temporal(supportMM) + amacrine.density.fitMMSq.temporal(supportMM),'xr')
xlabel('eccentricity [mm retina]');
ylabel('density [counts / sq mm]');
legend({'Curcio totalRGC','midget','parasol','bistratified','amacrine','model total all cells'});

% Plot cell volumes

% Volume of a sphere given diameter
sVol = @(d) 4/3*pi*(d./2).^3;

subplot(1,3,2)
plot(supportMM, sVol(midget.diameter.fitMM(supportMM)))
hold on
plot(supportMM, sVol(parasol.diameter.fitMM(supportMM)))
plot(supportMM, sVol(bistratified.diameter.fitMM(supportMM)))
plot(supportMM, sVol(amacrine.diameter.fitMM(supportMM)))
xlabel('eccentricity [mm retina]');
ylabel('individual cell volume [mm^3]');
legend({'midget','parasol','bistratified','amacrine'});


% Plot thickness

% Packing density of spheres
% We model the cells as spheres. For a single sphere size, the densest
% possible packing is Keppler's limit of ~0.74. In the RGC there are
% multiple cell classes. Except in the far periphery the two most numerous
% types of cells are the midget and parasol. If smaller cells are more
% numerous than larger cells, it is possible in principle to pack small
% cells into the spaces between larger cells, and achieve a density higher
% than Kepler's limit. When the ratio of the smaller to larger cell is
% greater than ~0.41, then the smaller is too large to pack into the
% larger. This is the situation we face in this model, so we assume
% Keppler's limit for sphere packing.
spherePackDensity = 0.74048048969;

% Need to assume a density higher than 1 to match Curcio's thickness
% measurements. This suggests that the Curcio measurements must have been
% subjected to some cell shrinkage.
spherePackDensity = 1.1;


subplot(1,3,3)
volumeProfile = (amacrine.density.fitMMSq.temporal(supportMM) .* sVol(amacrine.diameter.fitMM(supportMM)) + ...
    parasol.density.fitMMSq.temporal(supportMM) .* sVol(parasol.diameter.fitMM(supportMM)) + ...
    bistratified.density.fitMMSq.temporal(supportMM) .* sVol(bistratified.diameter.fitMM(supportMM)) + ...
    midget.density.fitMMSq.temporal(supportMM) .* sVol(midget.diameter.fitMM(supportMM))) ./ spherePackDensity;

plot(supportMM, volumeProfile);
hold on
plot(empiricalThickness.RGC.supportMM.temporal, empiricalThickness.RGC.thickMM.temporal, '*r');
legend({'Model thickness','Curcio RGC measure'});
xlabel('eccentricity [mm retina]');
ylabel('layer thickness [mm]]');




end % rgcThickness function


%% LOCAL FUNCTIONS

function midgetDensityFitMMSq = createMidgetDensityFunc(meridianAngle)

% The Barnett & Aguirre 2018 Linking Function parameters for each of the
% cardinal meridians (0, 90, 180, 270).
midgetFracLinkingParams = [...
   4.3286    1.6160; ...
   4.2433    1.5842; ...
   4.2433    1.5842; ...
   4.2433    1.6160];

% Define a support vector in retinal degrees
regularSupportPosDegRetina = 0:0.01:30;

% Obtain the spline fit function to total RGC density
fitRGCDensitySqDegRetina = getSplineFitToRGCDensitySqDegRetina(meridianAngle);

% Define a variable with RGC density over regular support and zero values
% at the optic disc positions
RGCDensityOverRegularSupport = ...
    zeroOpticDiscPoints(fitRGCDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngle);

% Define the mRGC density function using the Barnett & Aguirre 2018 linking
% function parameters
mRGCDensityOverRegularSupport = transformRGCToMidgetRGCDensityDacey(regularSupportPosDegRetina,RGCDensityOverRegularSupport,...
    'linkingFuncParams',midgetFracLinkingParams(meridianAngle/90+1,:));

% Perform a spline fit to the mRGC density expressed in mm and mm sq units
midgetDensityFitMMSq = fit(convert_degRetina_to_mmRetina(regularSupportPosDegRetina)', (mRGCDensityOverRegularSupport.*calc_degSqRetina_per_mmSqRetina())', 'smoothingspline');

end

function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegRetina, polarAngle)
opticDiscIndices = findOpticDiscPositions(regularSupportPosDegRetina, polarAngle);
vectorOut = vectorIn;
vectorOut(opticDiscIndices) = 0;
end

