function sizeData = prepareMeanCellSizeData_Liu(meridianSetName)
%
%   Liu.
%

supportDeg = [2.25 12.75];
meanSizeMicrons = [11.4 13.9]';


% Set up an empty structure
sizeData = struct('label',{},'angle',{},'supportDeg',{},'diameter',{});


% Depending upon the setName, we may include vertical and / or horizontal
if ismember(meridianSetName,{'horiz','both'})
    
    sizeData(end+1).label = 'temporal';
    sizeData(end).angle = 180;
    sizeData(end).supportDeg = supportDeg;
    sizeData(end).diameter = meanSizeMicrons./1000;
    
    sizeData(end+1).label = 'nasal';
    sizeData(end).angle = 0;
    sizeData(end).supportDeg = supportDeg;
    sizeData(end).diameter = meanSizeMicrons./1000;
    
end

if ismember(meridianSetName,{'vert','both'})
    
    sizeData(end+1).label = 'inferior';
    sizeData(end).angle = 270;
    sizeData(end).supportDeg = supportDeg;
    sizeData(end).diameter = meanSizeMicrons./1000;
    
    sizeData(end+1).label = 'superior';
    sizeData(end).angle = 90;
    sizeData(end).supportDeg = supportDeg;
    sizeData(end).diameter = meanSizeMicrons./1000;
    
end

end