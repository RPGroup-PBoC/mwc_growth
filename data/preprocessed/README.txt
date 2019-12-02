# Segmentation Files and Cell Statistics
This folder contains the data extracted from snapshots and the time-lapse
images of growing bacterial colonies used in this work. 

## Naming Format
All folders have the same naming pattern to make them easier to iterate over
computationally. The pattern is as follows:
```
[DATE]_r[RUN NUMBER]_[TEMP]C_[CARBON SOURCE]_O2_dilution
```

* `DATE`: The date of the experiment in YYYYMMDD format
* `RUN NUMBER`: This is the experimental run of that day. This will usually be
   `1` or `2` indicating the first and second runs of the day, respectively.
* `TEMP`: Growth and imaging temperature for the experiment in C. Either `32`,
    `37`, or `42`.
* `CARBON SOURCE`: Identity of the carbon source. Either `glucose`, `glycerol`,
   or `acetate`.
* O2: Identity of the operator. Only operator O2 was used in this work. 
* dilution: Description of the experiment method. All files here are `dilution`
   indicating that they monitored the dilution of repressors through division.


## Internal Structure
Each folder is filled with folders which themselves are full of folders. This is
the core structure of the output from the SuperSegger software. Below, we
describe the hierarchy and explain what each contains. 

Level 1: `growth` or `snaps` - The `growth` folder contains the time-lapse
    imaging experiment data while `snaps` contains the snapshots of each individual 
    ATC induction condition and the relevant controls.

Level 2.a: `growth/xy##` - Each folder represents a single xy position imaged for
    the growth series. Within each `xy` position is a single file `clist.mat` 
    which contains all cell information for the experiment, including the
    details needed to piece together the family trees.

Level 2.b: `snaps/[STRAIN]_[##]ngml` - Each folder contains the data from the
processed snapshots for each ATC induction condition. The [STRAIN] field is
either `auto` (autofluorescence), `delta` (constitutive expression control), 
or `dilution` (the experimental strain). The field `##ngml` is the ATC induction
concentration in ng / mL. 

Level 3: `xy##/clist.mat` - Each folder represents a single position. Within 
    each folder is a single `clist.mat` file which contains information about
    all of the segmented cells at that position.