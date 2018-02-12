% Segment the autofluorescence and the original images using SuperSegger
% in MATLAB
addpath(genpath('../../../../SuperSegger'));
% Define the experiment parameters.
DATE = '20180209';
BASENAME = 'tenjin_37C_glucose_O2';
samples = {'snaps', 'growth_0', 'growth_1', 'growth_2'}

CONST = loadConstants('100XEc');
CONST.parallel.PARALLEL_FLAG = 1;
CONST.trackFoci.numSpots = 0;
CONST.align.ALIGN_FLAG = 1;
CONST.trackOpti.REMOVE_STRAY = 1;
cleanFlag = 0;

for i=1:length(samples)
    disp(samples{i})
    statement = ['Beginning segmentaton ', num2str(i), ' out of ',...
        num2str(length(samples))];
    disp(statement)

    % Define the data directory.
    directory = ['../../../data/images/', DATE, '_', BASENAME,...
        '_dilution/' samples{i}];

    % Perform the segmentation.
    BatchSuperSeggerOpti(directory, 1, cleanFlag, CONST);
end
disp('Finished!');
