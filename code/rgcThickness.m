function rgcThickness(varargin )
% Caclulates RGC layer thickness by reference to constituent cell classes 
%
% Description:
%   Here goes some text.
%
% Inputs:
%
% Optional key / value pairs:
%
% Outputs:
%
% Examples:
%


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
fitRGCDensitySqDegRetina = getSplineFitToRGCDensitySqDegRetina(p.Results.polarAngle);
fitRGCDensitySqMMRetina = @(posMMretina) fitRGCDensitySqDegRetina(convert_mmRetina_to_degRetina(posMMretina)).*calc_degSqRetina_per_mmSqRetina();


%% Displaced amacrine cells
% Kara -- place here some text describing what we know about cell
% denisities for displaced amacrine cells, and the citations.
%
%   Curcio, Christine A., and Kimberly A. Allen. "Topography of ganglion
%   cells in human retina." Journal of comparative Neurology 300.1 (1990):
%   5-25.

% Data from Curcio & Allen 1990, Figure 10
amacrine.density.supportMM.temporal = [];
amacrine.density.countsMMSq.temporal = [];

% Amacrine cell body sizes
% Put comments here about where this info comes from
amacrine.diameter.supportMM = [];
amacrine.diameter.sizeMM = [];



%% Parasol RGCs
% Kara -- place here some text about parasol cells and the citation
% information for these

% Data from Silviera et al 1990, Figures 16 and 17
parasol.density.supportMM.temporal = 1:1:22;
parasol.density.countsMMSq.temporal = [21.5, 15.5, 9, 4.6, 3, 2.2, 1.7, 1.1, 0.7, 0.56, 0.4, 0.2, 0.1, nan, nan, nan, nan, nan, nan, nan, nan, nan];
parasol.density.supportMM.nasal = 1:1:22;
parasol.density.countsMMSq.nasal = [];

% Parasol cell body sizes
% Put comments here about where this info comes from
parasol.diameter.supportMM = [];
parasol.diameter.sizeMM = [];


%% Bistratified RGCs
bistratified.density.supportMM.temporal = [];
bistratified.density.countsMMSq.temporal = [];

% Parasol cell body sizes
% Put comments here about where this info comes from
bistratified.diameter.supportMM = [];
bistratified.diameter.sizeMM = [];


%% Midget RGCs

% These density functions will come from Geoff's code
midget.density.supportMM.temporal = [];
midget.density.countsMMSq.temporal = [];

% Parasol cell body sizes
% Put comments here about where this info comes from
midget.diameter.supportMM = [];
midget.diameter.sizeMM = [];


% Calculate thickness


end % rgcThickness function

