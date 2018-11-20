% Segment the autofluorescence and the original images using SuperSegger
% in MATLAB
% Define the experiment parameters.
addpath(genpath('../../../../../SuperSegger'));
DATE = '20181120';
BASENAME = 'r1_37C_glycerol_O2';
samples = {'growth'};

% Get the snaps names.
dirpath = ['../../../../data/images/', DATE, '_', BASENAME, '_dilution/snaps/'];
snap_files = dir(dirpath);
snap_samples = {snap_files.name};
ignored = {'.', '..', '.DS_Store', 'snaps'};
  for i=1:length(snap_samples)
      valid = 0;
      for j=1:length(ignored)
          valid = valid + strcmp(snap_samples{i}, ignored{j});
      end
      if valid == 0
          samples{end+1} = ['snaps/', snap_samples{i}];
          end
  end

CONST = loadConstants('100XEc');
CONST.trackFoci.numSpots = 0;
CONST.align.ALIGN_FLAG = 1;
CONST.trackOpti.REMOVE_STRAY = 1;
cleanFlag = 0;
% parpool(2, 'IdleTimeout', 1000000)
for i=1:length(samples)
    disp(samples{i})
    statement = ['Beginning segmentaton ', num2str(i), ' out of ',...
        num2str(length(samples))];
    disp(statement)

    % Define the data directory.
    directory = ['../../../../data/images/', DATE, '_', BASENAME, '_dilution/', samples{i}]


    % Perform the segmentation.
    BatchSuperSeggerOpti(directory, 1, cleanFlag, CONST);
delete(gcp('nocreate'));
end
disp('Finished!');
exit();
