function [squaredDifference, motherInt] = measuretriads(dataDir, nFiles, autoVal)
%MEAURETRIADS iterates througha  given directory of superSegger cell files
%to measure the squared difference sqDiff of the two daugher cells and the 
%intensity of the mother cell. 

% Get all of the files and load the fluorescence image. 
files = dir([dataDir 'cell*/*.mat']);
squaredDifference = 0;
motherInt = 0;
% Instantiate all of the storage vectors to be used. 
lineageFiles = {};
cellId = [];
matchedPairs = [];
cellIntensities = [];
cell1 = [];
cell2 = [];

counter = 0;  % Counter for unpredictable iteration. 
% Loop through all of the files and load the structure. 
for i=1:length(files)
    cellStruct = load([files(i).folder '/' files(i).name]);
    % Ensure that the cells 'died' on the last frame and were not present
    % in the original frame.
    if (cellStruct.death == nFiles) && (cellStruct.birth ~= 1)
        cellA = cellStruct.CellA{end};
        % Increase the counter and store this file as a valid cell.
        counter = counter + 1;
        lineageFiles{counter} = [files(i).folder '/' files(i).name];
        
        %Extract the intensity values and ID. 
        cellIntensities(counter) = cellA.fl2.sum - autoVal;
        cellId(counter) = cellStruct.ID;
    end
end

% Now loop through all of the files and determine which intensities are
% paired.
disp(length(lineageFiles))
for i=1:length(lineageFiles)
    cellStruct = load(lineageFiles{i});
    
    % Get the ID's. 
    sisterID = cellStruct.sisterID;
    ownID = cellStruct.ID;
    
    % Loop through each cell ID and make sure we aren't matching anything
    % up with itself or a pair we've already measured. 
    for j=1:length(cellId);
    if any(matchedPairs == ownID) || any(matchedPairs == sisterID)
        continue
    elseif cellId(j) == ownID;
        continue
    elseif cellId(j) == sisterID
        % Find the index of both and update the intensity values. 
        sisterIndex = find(cellId == sisterID);
        ownIndex = find(cellId == ownID);
        cell1 = [cell1 cellIntensities(ownIndex)];
        cell2 = [cell2 cellIntensities(sisterIndex)];
        
        % Add the daughter cell ID's to the matched pairs vector. 
        matchedPairs = [matchedPairs sisterID ownID];
    end
    end
    
% now compute the squqred differences an the total intensities. 
squaredDifference = (cell1 - cell2).^2;
motherInt = cell1 + cell2;
end

