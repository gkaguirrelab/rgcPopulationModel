function bistratified = bistratified( cardinalMeridianAngles, cardinalMeridianNames, totalRGC )
% Model of volume of bistratfiied ganglion cells in the RGC layer.
%   Digitized data from Dacey's publications concerning morphology of 
%   bistratified ganglion cells allows us to calculate how much the volume
%   of bistratified cells contributes to the overall RGC layer thickness
%   from data concerning bistratified cells density varying with eccentricity
%   across all meridians, and the cells' mean size.



%% Bistratified RGCs

% Bistratified proportion
% Data from Dacey 1993, Figure 13b:
%   Dacey, Dennis M. "Morphology of a small-field bistratified ganglion
%   cell type in the macaque and human retina." Visual neuroscience 10.6
%   (1993): 1081-1098.
%
% Dacey provided values for the temporal retina; we assume these values
% hold for all meridians
for mm = 1:length(cardinalMeridianAngles)
    bistratified.proportion.supportMM.(cardinalMeridianNames{mm}) = [.97, 1.96, 2.91, 3.92, 4.91, 5.92, 6.89, 7.84, 8.85, 9.86, 10.89, 11.9, 12.91, 13.84, 14.91];
    bistratified.proportion.value.(cardinalMeridianNames{mm}) = [.0135, .0168, .0202, .0241, .0284, .0324, .0364, .0403, .0447, .0485, .0538, .0573, .0603, .0641, .0662];
end

% Obtain a spline fit to the bistratified densities
for mm = 1:length(cardinalMeridianAngles)
    nonNanSupportIdx = ~isnan(bistratified.proportion.value.(cardinalMeridianNames{mm}));
    tmpSupport = bistratified.proportion.supportMM.(cardinalMeridianNames{mm})(nonNanSupportIdx)';
    bistratified.density.fitMMSq.(cardinalMeridianNames{mm}) = ...
        fit( [0; tmpSupport], ...        
        [0; totalRGC.density.fitMMSq.(cardinalMeridianNames{mm})(tmpSupport)' .* ...
        bistratified.proportion.value.(cardinalMeridianNames{mm})(nonNanSupportIdx)'], ...
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

end

