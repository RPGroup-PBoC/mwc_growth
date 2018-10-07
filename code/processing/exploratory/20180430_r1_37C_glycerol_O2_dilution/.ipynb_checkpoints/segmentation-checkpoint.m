% Segment the autofluorescence and the original images using SuperSegger
% in MATLAB
addpath(genpath('../../../../SuperSegger'));
% Define the experiment parameters.
DATE = '20180419';
BASENAME = '37C_glycerol_O2';
samples = {'growth_3'};

% Get the snaps names.
snap_files = dir(['../../../data/images/', DATE, '_', BASENAME,...
      '_dilution/','*snaps*']);
snap_samples = {snap_files.name};
ignored = {'.', '..', '.DS_Store', 'steady_state'};
  for i=1:length(snap_samples)

      valid = 0;
      for j=1:length(ignored)
          valid = valid + strcmp(snap_samples{i}, ignored{j});
      end
      if valid == 0
          samples{end+1} = snap_samples{i};
          end
  end

CONST = loadConstants('100XEc');
CONST.parallel.PARALLEL_FLAG = 1;
CONST.parallel.xy_parallel = 1;
CONST.trackFoci.numSpots = 0;
CONST.align.ALIGN_FLAG = 1;
CONST.trackOpti.REMOVE_STRAY = 1;
cleanFlag = 0;
samples = snap_samples;
for i=1:length(samples)
	parpool(24, 'IdleTimeout', 1000000)

    disp(samples{i})
    statement = ['Beginning segmentaton ', num2str(i), ' out of ',...
        num2str(length(samples))];
    disp(statement)

    % Define the data directory.
    directory = ['../../../data/images/', DATE, '_', BASENAME,...
        '_dilution/' samples{i}];

    % Perform the segmentation.
    BatchSuperSeggerOpti(directory, 1, cleanFlag, CONST);
delete(gcp('nocreate'));
end
disp('Finished!');
exit();
