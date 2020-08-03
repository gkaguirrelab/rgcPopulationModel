function makeMMperDegPolyFit()

% This code should live in the data directory, so save the mat file here
filepath = fileparts(mfilename());
saveName = fullfile(filepath,'mmPerDegPolyFitEmmetrope.mat');

% Assemble the corneal curvature, spherical error, and axial length.
axialLength = 23.8;
SR = 0;

% Create the model eye
eye = modelEyeParameters('axialLength',axialLength,'sphericalAmetropia',SR,'calcLandmarkFovea',true);

% Define the visual field domain over which we will make the measure
horizVals = -30:15:30;
vertVals = -30:15:30;

% Define an empty matrix to hold the results
mmPerDeg = nan(length(horizVals),length(vertVals));

% Define the delta deg
deltaDegEuclidean = 1;
deltaAngles = [sqrt(deltaDegEuclidean/2) sqrt(deltaDegEuclidean/2)];

% Loop over horizontal and vertical field positions
for jj = 1:length(horizVals)
    for kk = 1:length(vertVals)
        % The position in the field relative to the optical axis of the
        % eye
        degField = [horizVals(jj) vertVals(kk)] + eye.landmarks.fovea.degField(1:2);
        % Obtain the retinal points that are delta degrees on either
        % side of the specified degree field position
        [~,X0] = calcRetinaFieldPoint( eye, degField - deltaAngles./2);
        [~,X1] = calcRetinaFieldPoint( eye, degField + deltaAngles./2);
        mmPerDeg(jj,kk) = norm(X0-X1) / norm(deltaAngles);
    end
end
% Fit a polynomial surface to the measure
[X,Y]=meshgrid(horizVals,vertVals);
mmPerDegPolyFit = fit([X(:),Y(:)],mmPerDeg(:),'poly33');

save(saveName,'mmPerDegPolyFit');

end