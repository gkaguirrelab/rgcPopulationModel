function ipRGC = ipRGC( cardinalMeridianAngles, cardinalMeridianNames )

%% ipRGCs
% Citation info here

for mm = 1:length(cardinalMeridianAngles)
    % Data from Curcio & Allen 1990, Figure 10, used for all meridians
    ipRGC.density.supportMM.(cardinalMeridianNames{mm}) = [0, 0.11, 0.23, 0.4, 0.65, 0.86, 1.45, 2.46, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.6, 12.6, 13.6, 14.6, 15.6, 16.6, 17.7, 18.65, 19.65];
    ipRGC.density.countsMMSq.(cardinalMeridianNames{mm}) = [73, 124, 317, 542, 689, 813, 1052, 1112, 1159, 1187, 1146, 1069, 992, 942, 888, 848, 798, 767, 731, 699, 677, 650, 642, 652, 676];
    % Obtain a spline fit to the ipRGC densities
    ipRGC.density.fitMMSq.(cardinalMeridianNames{mm}) = ...
        fit(ipRGC.density.supportMM.(cardinalMeridianNames{mm})', ipRGC.density.countsMMSq.(cardinalMeridianNames{mm})', 'smoothingspline');
end

% Size info here
% Digitzed the line fit to the data in Figure:
ipRGC.diameter.supportMM = [0, 0.87, 1.68, 2.80, 3.76, 4.74, 5.86, 6.68, 7.36, 8.19, 8.71, 9.45, 9.97, 10.67, 11.45, 12.39, 13.47, 14.5, 15.62, 16.54, 17.82];
ipRGC.diameter.sizeMM = [.01259, .01266, .01274, .01286, .01303, .0132, .01324, .01332, .01336, .01335, .01353, .01361, .01365, .01373, .01381, .01394, .01402, .01419, .01423, .01439, .01452];
ipRGC.diameter.fitMM = fit(ipRGC.diameter.supportMM',ipRGC.diameter.sizeMM','smoothingspline');

end

