function sizeData = prepareMeanCellSizeData()
%
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

sizeData(1).label = 'temporal';
sizeData(1).angle = 180;
sizeData(1).supportDeg = supportDeg;
sizeData(1).diameter = meanSizeMicrons./1000;

sizeData(2).label = 'nasal';
sizeData(2).angle = 180;
sizeData(2).supportDeg = supportDeg;
sizeData(2).diameter = meanSizeMicrons./1000;

end