function totalRGC = totalRGC( cardinalMeridianAngles, cardinalMeridianNames )
%
%% Total RGC density function
%   Create a function to return total RGC densities per square mm of retina
%   at a given eccentricity position in mm retina.
%   To do so, we call the rgcDisplacementMap functions which have a
%   representation of the Curcio & Allen 1990 RGC density results. The
%   rgcDisplacementMap toolbox is referenced in eccentricity units of degrees
%   retina, and provides densities in square degrees retina. We convert those
%   values here to mm and square mm.
%
for mm = 1:length(cardinalMeridianAngles)
    tmpFit = getSplineFitToRGCDensitySqDegRetina(cardinalMeridianAngles(mm));
    totalRGC.density.fitMMSq.(cardinalMeridianNames{mm}) = @(posMMretina) tmpFit(convert_mmRetina_to_degRetina(posMMretina))'.*calc_degSqRetina_per_mmSqRetina();
end



end

