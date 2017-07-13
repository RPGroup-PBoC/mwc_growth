%%Compute the autofluorescence value.
dataDir = '../../data/mbl_2017/images/supersegger_practice/20170608_fullset/concs/'

% Separate the CFP and mCherry files.
autoCFP = dir([dataDir '*auto*c2*.tif']);
autoRFP = dir([dataDir '*auto*c4*.tif']);
im = imread([autoCFP(1).folder '/' autoCFP(1).name]);
imNorm = mat2gray(im);
imFilt = medfilt2(imNorm);
imSeg = imFilt > 0.45;
imLarge = bwareaopen(imSeg, 100);
imshowpair(imFilt, imLarge)

intensities = [];
for i=1:length(autoRFP)
    % Load the image and segment. 
    im = imread([autoCFP(i).folder '/' autoCFP(i).name]);
    imRFP = imread([autoRFP(i).folder '/' autoRFP(i).name]);
    % Segment. 
    imSeg = simpleseg(im, 0.45);
    
    % Compute the properties.
    props = regionprops(imSeg, imRFP,...
        'Area', 'MeanIntensity');
    
    areas = [props.Area];
    meanInts = [props.MeanIntensity];
    totalInt = areas .* meanInts;
    intensities = [intensities totalInt];
end

meanAuto = mean(intensities);
% Loop through all images and compute the intensity. 

% 
%%

% Define the data directory. 
dataDir = '../../data/mbl_2017/images/supersegger_practice/20170608_fullset/x*';
files = dir(dataDir);
sqDiff = [];
itot = [];
for i=3:length(files)
    if files(i).isdir
    cellDir = [dataDir(1:end-2) files(i).name '/'];
    [x, y] = measuretriads(cellDir, 22, 0);
    sqDiff = [sqDiff x];
    itot = [itot y];
    end
end

%%
% Make it a table
squared_diff = sqDiff';
I_tot = itot';
data = table(squared_diff, I_tot);

% Sort the values. 
sortedVals = sortrows(data, 'I_tot');


% Binning strategy -- get 50 per bin.
alpha = [];
binned_Itot = [];
binned_sqdiff = [];
binSizeList = 1:100:2000;
for b=1:length(binSizeList)
for i=1:binSizeList(b):(length(sortedVals.squared_diff) - binSizeList(b))
    binned_Itot =[binned_Itot mean(sortedVals.I_tot(i:i+binSizeList(b)))];
    binned_sqdiff = [binned_sqdiff mean(sortedVals.squared_diff(i:i+binSizeList(b)))];
    residFunc = @(alphaVal) (log(binned_sqdiff) - alphaVal .* log(binned_Itot));
    % Fit the alpha parameter. 
    fitVal = lsqnonlin(residFunc, 100);
    alpha(b) = fitVal(1);
end

end

% Plot the cal factor as a function of bin size. 
plot(binSizeList(1:end -3), alpha, 'o')
xlabel('events per bin');
ylabel('calibration factor (a.u. / molecule)');

%% Bin with 50 points. 
binSize = 100;
binned_Itot = [];
binned_sqdiff = []
for i=1:binSize:length(sortedVals.squared_diff) - binSize
    binned_Itot =[binned_Itot mean(sortedVals.I_tot(i:i+binSize))];
    binned_sqdiff = [binned_sqdiff mean(sortedVals.squared_diff(i:i+binSize))];
end

residFunc = @(alphaVal) (log(binned_sqdiff) - alphaVal.*log(binned_Itot));
fitVal = lsqnonlin(residFunc, 1000);

loglog(itot, sqDiff, 'o');
hold on
loglog(binned_Itot, binned_sqdiff, 'ro', 'MarkerFaceColor', 'r');
itot_range = logspace(4, 6, 500);
y = exp(fitVal) * itot_range;
loglog(itot_range, y, 'r-');
hold off
xlabel('I_{tot}')
ylabel('(I_1 - I_2)^2')


%% Compute the auto copy number. 
autoR = intensities ./ fitVal;
hist(autoR, 100)

%% COmpute it for 10ng
%%Compute the autofluorescence value.
dataDir = '../../data/mbl_2017/images/supersegger_practice/20170608_fullset/concs/'

% Separate the CFP and mCherry files.
autoCFP = dir([dataDir '*10*c2*.tif']);
autoRFP = dir([dataDir '*10*c4*.tif']);
im = imread([autoCFP(1).folder '/' autoCFP(1).name]);
imNorm = mat2gray(im);
imFilt = medfilt2(imNorm);
imSeg = imFilt > 0.15;
imLarge = bwareaopen(imSeg, 100);
imshowpair(imFilt, imLarge)
%%
intensities = [];
for i=1:length(autoRFP)
    % Load the image and segment. 
    im = imread([autoCFP(i).folder '/' autoCFP(i).name]);
    imRFP = imread([autoRFP(i).folder '/' autoRFP(i).name]);
    % Segment. 
    imSeg = simpleseg(im, 0.15);
    
    % Compute the properties.
    props = regionprops(imSeg, imRFP,...
        'Area', 'MeanIntensity');
    
    areas = [props.Area];
    meanInts = [props.MeanIntensity];
    totalInt = areas .* meanInts;
    intensities = [intensities totalInt];
end

mean10 = mean(intensities - meanAuto);
% Loop through all images and compute the intensity. 