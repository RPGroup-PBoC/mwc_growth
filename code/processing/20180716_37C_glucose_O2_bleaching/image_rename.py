# %%
import glob
import os
import shutil

DATE = 20180716
TEMP = 37
CARBON = 'glucose'
OPERATOR = 'O2'

data_dir = '../../../data/images/{}_{}C_{}_{}_bleaching'.format(DATE, TEMP, CARBON, OPERATOR)


# Get all bleaching files
bleaching_files = glob.glob('{}/photobleaching/*mCherry*.TIF'.format(data_dir))
unique_pos = []
for f in bleaching_files:
    _, _, pos, _, _ = f.split('/')[-1].split('_')
    pos = pos.split('ml')[1]
    if pos not in unique_pos:
        unique_pos.append(pos)

for i, pos in enumerate(unique_pos):
    # Grab particular files
    files = glob.glob('{}/photobleaching/*ml{}_*mCherry*.TIF'.format(data_dir, pos))

    # Get the information.
    for j, f in enumerate(files):
        _, _, conc, _, time = f.split('/')[-1].split('_') 
        time = int(time.split('t')[1].split('.')[0])
        new_name = '{}/photobleaching/{}_induction_8ngml_xy{:02d}_t{:02d}.tif'.format(data_dir, DATE, i, time)
        os.rename(f, new_name)