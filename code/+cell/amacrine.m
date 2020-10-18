function amacrine = amacrine(cellSizeParams, showPlots)
% Size and count functions for the amacrine cell class
%
% Syntax:
%  amacrine = cell.amacrine(showPlots)
%
% Description:
%   Returns an array of structures, where each structure has a handle to a
%   function that returns cell counts (per degree squared) and cell
%   diameter (in mm) as a function of retinal eccentricity (in degrees).
%
% Examples:
%{
    amacrine = cell.amacrine(true);
%}

% Handle plotting
if nargin==0
    showPlots = false;
end

% The meridians over which the calculation is to be performed
cardinalMeridianAngles = [0 90 180 270];
cardinalMeridianNames = {'nasal' 'superior' 'temporal' 'inferior'};

%% Cell counts
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

% Loop over the specified meridians
for mm = 1:length(cardinalMeridianAngles)
    % Data from Curcio & Allen 1990, Figure 10, used for all meridians
    supportMM = [0, 0.11, 0.23, 0.4, 0.65, 0.86, 1.45, 2.46, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.6, 12.6, 13.6, 14.6, 15.6, 16.6, 17.7, 18.65, 19.65];
    countsMMSq = [73, 124, 317, 542, 689, 813, 1052, 1112, 1159, 1187, 1146, 1069, 992, 942, 888, 848, 798, 767, 731, 699, 677, 650, 642, 652, 676];

    % Hack the data to set the value at zero equal to zero, so that this
    % will accord with the thickness measures
    countsMMSq(supportMM==0)=0;
    
    % Convert retinal mm to visual degrees
    supportDeg = convert_mmRetina_to_degVisual(supportMM, cardinalMeridianAngles(mm));
    countsDegSq = countsMMSq .* calc_mmSqRetina_per_degSqVisual(supportDeg, cardinalMeridianAngles(mm));
    
    % Obtain a spline fit to the cell densities
    splineFit = fit(supportDeg', countsDegSq', 'cubicinterp');
    
    % Set up this meridian model element
    amacrine(mm).label = cardinalMeridianNames(mm);
    amacrine(mm).angle = cardinalMeridianAngles(mm);
    
    % Nan optic disc points and save the anonymous function
    amacrine(mm).countsDegSq = @(posDeg) ...
        nanOpticDiscPoints(splineFit(posDeg), posDeg, cardinalMeridianAngles(mm));
    
    if showPlots
        if mm == 1
            figure
        end
        plot(0:0.5:50,amacrine(mm).countsDegSq(0:0.5:50));
        hold on
        plot(supportDeg,countsDegSq,'*');
        title('amacrine cell density by meridian');
    end
    
end


%% Cell diameters
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
for mm = 1:length(cardinalMeridianAngles)

    % Data from Dacey 1990, Figure 5B, used for all meridians
	supportMM = [0, 0.87, 1.68, 2.80, 3.76, 4.74, 5.86, 6.68, 7.36, 8.19, 8.71, 9.45, 9.97, 10.67, 11.45, 12.39, 13.47, 14.5, 15.62, 16.54, 17.82];
	sizeMM = [.01259, .01266, .01274, .01286, .01303, .0132, .01324, .01332, .01336, .01335, .01353, .01361, .01365, .01373, .01381, .01394, .01402, .01419, .01423, .01439, .01452];
    meanSize = mean(sizeMM);

    % Convert retinal mm to visual degrees
    supportDeg = convert_mmRetina_to_degVisual(supportMM,cardinalMeridianAngles(mm));
    meanSupport = mean(supportDeg);

    % Model the size as mean with proportional growth slope
    myCellSize = @(x) (meanSize + meanSize.*((x-meanSupport).*cellSizeParams(1)+cellSizeParams(2)) )';

    % Obtain the fit and save
    amacrine(mm).diameter = myCellSize;
    
    if showPlots
        if mm == 1
            figure
        end
        plot(supportDeg,sizeMM,'*');
        hold on
        plot(0:0.5:70,amacrine(mm).diameter(0:0.5:70));
        title('amacrine soma size by eccentricity');
    end
    
end

