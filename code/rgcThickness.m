function rgcThickness( )
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
p.addParameter('sceneGeometryFileName','', @(x)(isempty(x) | ischar(x)));
p.addParameter('intrinsicCameraMatrix',[2600 0 320; 0 2600 240; 0 0 1],@isnumeric);
p.addParameter('sensorResolution',[640 480],@isnumeric);
p.addParameter('radialDistortionVector',[0 0],@isnumeric);
p.addParameter('cameraTranslation',[0; 0; 120],@isnumeric);
p.addParameter('cameraTorsion',0,@isnumeric);
p.addParameter('constraintTolerance',0.02,@isscalar);
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('medium','air',@ischar);
p.addParameter('spectralDomain','nir',@ischar);
p.addParameter('forceMATLABVirtualImageFunc',false,@islogical);

% parse
p.parse(varargin{:})





% Define cell proportions


% Define cell sizes


% Calculate thickness


end % rgcThickness function

