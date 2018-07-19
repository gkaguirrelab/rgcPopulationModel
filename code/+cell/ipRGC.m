function ipRGC = ipRGC( cardinalMeridianAngles, cardinalMeridianNames )

%% ipRGCs
%   Dacey, Dennis M., et al. "Melanopsin-expressing ganglion cells in
%   primate retina signal colour and irradiance and project to the LGN."
%   Nature 433.7027 (2005): 749.
%
% Data from staining instrinsically photosensitive retinal ganglion cells
% in human retinas. Note that 60% of all ipRGCs have outer stratifying
% dendritic fields, and 40% of these outer stratifying ipRGCs have their
% cell bodies lying within the inner plexiform layer. Thus, onlly 76% of
% the the total ipRGC count (1-0.6*0.4) should be considered as residing
% within the retinal ganglion cell layer

ipRGCpropInRGCLayer = 0.76;

% Density info
for mm = 1:length(cardinalMeridianAngles)
    % Data from Dacey 2005, Figure 1G, used for all meridians
    ipRGC.density.supportMM.(cardinalMeridianNames{mm}) = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17];
    ipRGC.density.countsMMSq.(cardinalMeridianNames{mm}) = ipRGCpropInRGCLayer.*[21.06, 15.36, 12.02 8.27, 5.95, 8.46, 5.87, 5.96, 5.31, 5.78, 4.99, 4.67, 7.55, 4.86, 6.21, 4.03];
    % Obtain a spline fit to the ipRGC densities
    ipRGC.density.fitMMSq.(cardinalMeridianNames{mm}) = ...
        fit(ipRGC.density.supportMM.(cardinalMeridianNames{mm})', ipRGC.density.countsMMSq.(cardinalMeridianNames{mm})', 'smoothingspline');
end

% Size info here
ipRGC.diameter.supportMM = [];
ipRGC.diameter.sizeMM = [];
ipRGC.diameter.fitMM = fit(ipRGC.diameter.supportMM',ipRGC.diameter.sizeMM','smoothingspline');

end

