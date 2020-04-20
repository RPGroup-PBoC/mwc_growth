---
status: rejected
reason:
  This experiment was an abject failure. There were very few division events
  properly segmented courtesy of Z drift during acquisition of the growth movie.
---

# 2018-02-20 Dilution Experiment

## Purpose
This experiment served as another diagnostic test of using an LED illumination
source for the dilution experiment.

## Strain Information

| Location | Plasmid | Genotype | Host Strain | Shorthand |
| :------- | :------ | :------- | ----------: | --------: |
| dilution 1 - pos. 2 | None | `galK<>KD4+noKan, gspI<>4*5noO1v1-cfp` | HG105 | `autofluorescence` |
| dilution 1 - pos. 5 | `pZS3*PN25-tetR` | `galK<>25O2+11-YFP ybcN<>1-Wiggins2-lacI-mCherry, gspI<>4*5noO1v1-CFP` | HG105 | `dilution`|
| dilution 1 - pos. 9 | `pZS3*PN25-tetR`| `galK<>25O2+11-YFP, gspI<>4*5noO1v1-CFP` |  HG105 |`deltaLacI` |

## Titration Series

| Inducer | Concentration |
| :------ | ------------:
| Anhydrotetracycline HCl (ATC) | 0, 1, 2, 3, 4, 6, 10  [ng / mL] |
| Isopropylthiogalactopyranoside (IPTG) | 0 [mM] |

## Notes & Observations
* The cells were harvested for imaging at an OD<sub>600nm</sub> of ~0.4.
* There was considerable Z drift during the growth movie acquisition. This resulted
in many out-of-focus fields where the lineage tracing could not be performed.
* Position 4 of the growth movie was removed from consideration as it was a duplicated
field of Position 3.

## Analysis Files

**Calibration Factor Determination**
![](output/20180220_tenjin_37C_glucose_O2_calibration_factor)

**Fold-change**
Not Calculated

## Experimental Protocol

1. Cells as described in "Strain Information" were grown to saturation overnight in 3mL of LB Miller + chloramphenicol for the `dilution` strain. The cells were assumed to be saturated during this time.

2. Cells were diluted 1:1000 into 3mL of M9 + 0.5% Glucose (+ chloramphenicol for the `dilution` strain) in 14mL Falcon tubes. ATC was added from a 10µg/mL stock in EtOH to the appropriate concentration.

3. Tubes were placed in a rack and covered with a plastic box to protect from photocleavage of ATC. Cells were allowed to grow for ~ 8 hours at 37°C with shaking at ~ 220 RPM.

4. Once the cells reached an OD<sub>600nm</sub> between 0.2 - 0.4, the cells were removed from the warm room. A 200µL aliquot of the `dilution` samples from ATC concentrations of 1, 2, 3, 4, 6, and 8 ng/mL were combined in a 1.5mL eppendorf tube.

5. The `dilution`  was pelleted at 13000xg for 2 min. The supernatant was withdrawn and the pellet was resuspended in 1mL of ATC-free M9 + 0.5% glucose. This procedure was repeated twice more.

6. The washed `dilution` mixture was diluted 1:5 into ATC-free M9 + 0.5% glucose. Aliquots of 1µL were spotted onto 2% agarose pads made of M9 + 0.5% glucose.

7. Aliquots of 1µL from the other samples (`autofluorescence`, `deltaLacI`, and `dilution` for all ATC concentrations) were added to agarose pads.

8. Agarose pads spotted with cells were allowed to dry for 10 - 15 min and were then placed onto UV sterilized glass bottom dishes.

9. After mounting, the sample dish was affixed to the microscope using double stick tape. Between five and ten positions were marked per snapshot sample. Exposures were as follows:
    - Brightfield - 80ms
    - mCherry - 4000ms
    - YFP - 500ms

10. Between 20 and 30 positions were then marked on the `dilution` mixture pad. These positions were chosen requiring separation of cells and avoidance of debris.

11. These were positions were imaged every five minutes for two hours using only the Brightfield channel. After two hours, these positions were imaged once more using a Brightfield, mCherry, and YFP channels.

12. Between 15 and 20 positions were marked on a homogeneously fluorescent slide. These positions were imaged in both YFP and mCherry with the following exposures:
    - mCherry - 100ms
    - YFP - 5ms

13. Using the same positions marked on the  fluorescent slide, between 15 and 20 positions were imaged in mCherry and YFP with no illumination reaching the camera. These exposures were
    - mCherry - 4000ms
    - YFP - 500ms

14. The samples were discarded and the dataset was transferred to the storage server.


## Analysis Protocol
1. All image files were transferred to the local cluster (`delbrück`) for segmentation through SuperSegger.

2. The images were then renamed to match the required naming convention. These images were flattened using the images of the fluorescent slide and camera shot noise.

3. Blank images were created and saved for the growth series for all time points except the
terminal fluorescence image. This was performed such that SuperSegger would properly extract the fluorescence information from the final image. Fluorescence images are needed at all time steps for this to work.

4. After segmentation was complete, All `clist.mat` files were transferred to a local computer to analyze. These files were loaded and filtered based on size and aspect ratio. The fluctuations were computed and the best-fit for the calibration factor was determined via MCMC.