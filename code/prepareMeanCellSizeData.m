function sizeData = prepareMeanCellSizeData(meridianSetName)
%
%   Hebel, R., and H. Holländer. "Size and distribution of ganglion cells
%   in the human retina." Anatomy and embryology 168.1 (1983): 125-136.
%
% Table 1

supportDeg = [0 5 10 15 20 30 40 50 60 80];
retina1 = [10.36 10.35 11.32 12.36 13.13 13.69 13.50 13.37 13.01 12.31];
retina2 = [ 9.82 10.20 10.85 11.39 12.85 13.29 13.48 12.87 13.01 12.89];
retina3 = [ 9.47  9.58 10.45 11.12 12.59 13.03 13.08 12.39 12.41 11.27];

meanSizeMicrons = mean([retina1; retina2; retina3])';

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