function totalRGC = totalRGC( cardinalMeridianAngles, cardinalMeridianNames )
%
%% Total RGC density function
%   Create a function to return total RGC densities per square mm of retina
%   at a given eccentricity position in mm retina. To do so, we call the
%   rgcDisplacementMap functions which have a representation of the Curcio
%   & Allen 1990 RGC density results. The rgcDisplacementMap toolbox is
%   referenced in eccentricity units of degrees visual field, and provides
%   densities in square degrees visual field.

for mm = 1:length(cardinalMeridianAngles)
    tmp = getSplineFitToRGCDensitySqDegVisual(cardinalMeridianAngles(mm));
    totalRGC.density.fitDegSq.(cardinalMeridianNames{mm}) = @(x) tmp(x)';
end



end

