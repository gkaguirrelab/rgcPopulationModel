

outDir = '~/Desktop/KaraCloud_VSS2019_figs';
mkdir(outDir);


%% Figure 1
% The RGC+IPL thickness map across subjects

% Find the data and load it
thisFuncLocation = mfilename('fullpath');
rgcFilePath = fullfile(fileparts(fileparts(thisFuncLocation)),'data','rgcIplThicknessMap.mat');
load(rgcFilePath,'rgcIplThicknessMap');

% Get the size and the support for the map
sizeX = size(rgcIplThicknessMap,1);
sizeY = size(rgcIplThicknessMap,2);
supportDeg = linspace(-15,15,sizeX);
rgciplThickHorizontalMeridian=rgcIplThicknessMap(sizeX/2,:);

% Save the figure twice, once with the "Painters" renderer to capture the
% PDF of the mesh
renderOptions = {'OpenGL','Painters'};
fileSuffix = {'.png','.pdf'};
for ii = 1:2
    fig1=figure(1);
    fig1.Renderer=renderOptions{ii};
    
    if strcmp(fileSuffix{ii},'.png')
        surf(supportDeg(1:sizeX),supportDeg(1:sizeY),rgcIplThicknessMap,'EdgeColor','none')
        alpha 0.95
    else
        surf(supportDeg(1:2),supportDeg(1:2),rgcIplThicknessMap(1:2,1:2),'EdgeColor','none')
        hold on
        plot3(supportDeg(1:sizeX),zeros(1,sizeX),rgciplThickHorizontalMeridian,'--k','LineWidth',2)
    end
    xlim([-15 15]);
    ylim([-15 15]);
    xlabel('temporal    visual angle [deg]      nasal');
    ylabel('superior    visual angle [deg]      inferior');
    zlabel('microns');
    pbaspect([.8 .8 .11])
    grid off
    set(gca,'Color','none')
    colorbar
    filename = fullfile(outDir,['rgcipl3D',fileSuffix{ii}]);
    switch fileSuffix{ii}
        case '.png'
            print(fig1,filename,'-dpng','-r600')
        case '.pdf'
            print(fig1,filename,'-dpdf')
    end
    close(fig1)
end

%% Figure 2
% Profile of RGC and IPL thickness on the horizontal meridian
fig2 = figure(2);

% Plot the RGC+IPL empirical thickness
plot(supportDeg(1:sizeX),rgciplThickHorizontalMeridian,'-k','LineWidth',2)
hold on
xlim([-15 15]);
ylim([0 125]);

% Add the IPL thickness, and calculated RGC thickness, for the temporal and
% nasal arms
ratioFuncByThickness = rgcLayerProportion;
temporalArmRGCProportion = fliplr(ratioFuncByThickness(linspace(0,15,sizeX/2),fliplr(rgciplThickHorizontalMeridian(1:sizeX/2))./1000));
temporalArmIPLthick = rgciplThickHorizontalMeridian(1:sizeX/2) - temporalArmRGCProportion.*rgciplThickHorizontalMeridian(1:sizeX/2);
plot(linspace(-15,0,sizeX/2),temporalArmIPLthick,'-r')
plot(linspace(-15,0,sizeX/2),rgciplThickHorizontalMeridian(1:sizeX/2)-temporalArmIPLthick,'-b');

nasalArmRGCProportion = ratioFuncByThickness(linspace(0,15,sizeX/2),rgciplThickHorizontalMeridian(sizeX/2+1:end)./1000);
nasalArmIPLthick = rgciplThickHorizontalMeridian(sizeX/2+1:end) - nasalArmRGCProportion.*rgciplThickHorizontalMeridian(sizeX/2+1:end);
plot(linspace(0,15,sizeX/2),nasalArmIPLthick,'-r')
plot(linspace(0,15,sizeX/2),rgciplThickHorizontalMeridian(sizeX/2+1:end)-nasalArmIPLthick,'-b');

% Format
xlabel('temporal    visual angle [deg]      nasal');
ylabel('microns');
pbaspect([3 1 1])
box off
grid off
filename = fullfile(outDir,['rgcIPL_2Dplot',fileSuffix{ii}]);
print(fig2,filename,'-dpdf')
close(fig2)


%% Figure 3
% Dacey and Drasdo midget fraction functions
fig3 = figure();
regularSupportPosDegVisual = 0:0.1:15;
totalRGC = cell.totalRGC([180], {'temporal'});
rgcDensitySqDegVisual = totalRGC.density.fitDegSq.temporal(regularSupportPosDegVisual');
[ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, rgcDensitySqDegVisual );
plot(regularSupportPosDegVisual,midgetFraction,'-k');
hold on
[ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDrasdo( regularSupportPosDegVisual, rgcDensitySqDegVisual );
plot(regularSupportPosDegVisual,midgetFraction,'--k');
xlabel('Eccentricity [deg]');
xlim([0 15]);
ylabel('Midget Fraction');
ylim([0 1]);
pbaspect([1 2 1])
filename = fullfile(outDir,['midgetFractionModels.pdf']);
print(fig3,filename,'-dpdf')
close(fig3)
