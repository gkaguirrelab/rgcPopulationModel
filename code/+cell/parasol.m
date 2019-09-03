function parasol = parasol( totalRGC, midget, bistratified, showPlots )
% Size and count functions for the parasol RGC class
%
% Syntax:
%  parasol = cell.parasol( totalRGC, midget, bistratified, showPlots )
%
% Description:
%   Returns an array of structures, where each structure has a handle to a
%   function that returns cell counts (per degree squared) and cell
%   diameter (in mm) as a function of retinal eccentricity (in degrees).
%
% Examples:
%{
    totalRGC = cell.totalRGC();
    midget = cell.midget( totalRGC );
    bistratified = cell.bistratified(totalRGC);
    parasol = cell.parasol( totalRGC, midget, bistratified, true );
%}

% Handle plotting and missing Params
if nargin==3
    showPlots = false;
end

% Define a maximum eccentricity of the model
maxEccenDeg = 50;


%% Cell counts
% We have solid information of total RGC density and an estimate of midget
% RGC density. The only parasol density measurements are from macaque. We
% therefore infer the parasol density by assuming that:
%   densityParasol = densityTotalRGC - densityMidget - densityBistratified

% Define a support vector in visual degrees
supportDeg = 0:0.01:maxEccenDeg;

% Loop over the specified meridians
for mm = 1:length(totalRGC)
    
    % Set up this meridian model element
    parasol(mm).label = totalRGC(mm).label;
    parasol(mm).angle = totalRGC(mm).angle;
    
    % Create a function for parasol density
    parasol(mm).countsDegSq = @(x) ...
        totalRGC(mm).countsDegSq(x')- ...
        midget(mm).countsDegSq(x') - ...
        bistratified(mm).countsDegSq(x');

    % Plot the fit
    if showPlots
        if mm == 1
            figure
        end
        subplot(ceil(length(totalRGC)/2),ceil(length(totalRGC)/2),mm)
        plot(supportDeg,totalRGC(mm).countsDegSq(supportDeg),'-k');
        hold on
        plot(supportDeg,parasol(mm).countsDegSq(supportDeg),'-r');
        title(totalRGC(mm).label);
        ylim([0 2500]);
    end

end


%% Cell diameters
% The first three values are taken from Figure 4B
%
%   Liu, Zhuolin, et al. "Imaging and quantifying ganglion cells and other
%   transparent neurons in the living human retina." Proceedings of the
%   National Academy of Sciences 114.48 (2017): 12803-12808.
%
% The 4th value is taken from Table 1 of:
%
%   Dacey, Dennis M. "Morphology of a small-field bistratified ganglion
%   cell type in the macaque and human retina." Visual neuroscience 10.6
%   (1993): 1081-1098.
%
% Lacking better data, I took the mean human parasol diameter value
% reported by Dacey and assigned this to 30 degrees of eccentricity.
%

% Loop over the specified meridians
for mm = 1:length(totalRGC)

    % Support in the source data is in degrees of visual field along the
    % temporal retina
    supportDeg = [6.75 8.75 12.25 30];
	sizeMM = [0.01863 0.01888 0.02095 0.0233];
    
    % Obtain the fit and save. We find that an exponential does best with
    % the impoverished set of measurements we have
    fx = @(a,b,c,x) (a.*x).^b+c;

    parasol(mm).diameter = fit(supportDeg', sizeMM',...
        fx,'StartPoint', [0.001 1 0.018], ...
        'Lower', [0 1 0]);

    % Plot the fit
    if showPlots
        if mm == 1
            figure
        end
        plot(supportDeg,sizeMM,'x');
        hold on
        plot(0:0.01:maxEccenDeg,parasol(mm).diameter(0:0.01:maxEccenDeg));
    end

end


end