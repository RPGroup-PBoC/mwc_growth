---
status: Accepted 
reason:  All went as expected. 
---

# 2018-10-18 Microscopy Growth Measurement

## Purpose
This experiment was performed to measure the growth rate of our cells on M9 + 0.4% glycerol.

## Strain Information
| Location | Plasmid | Genotype | Host Strain | Shorthand |
| :------- | :------ | :------- | ----------: | --------: |
| Dilution 1 - 9| - | galK<>2*O2+11-YFP, ybcN<>1-Wiggins2-lacI, gspI<>NoO1v1-CFP| HG105 | delta |

## Titration Series

| Inducer | Concentration |
| :------ | ------------: |
| Anhydrotetracycline HCl (ATC) | 0 |
| Isopropylthiogalactopyranoside (IPTG) | 0|

## Notes & Observations
* The sample was diluted 10x from a culture that had been overgrown in the morning and had been diluted 5x for an induction experiment.

## Analysis Files

**Growth Movie**
![growth movie](output/growth_movie.gif)

**Growth Curves and Rate Determination**
[![growth curve](output/growth_rate.png)]

## Experimental Protocol

1. Cells as described in "Strain Information" were grown to saturation in 3mL of LB Miller.

2. Cells were diluted 1:1000 into 3mL of growth medium in a 14mL Falcon tube.

3. Cells were allowed to grow for ~12.5 hours at 37°C with shaking at ~ 220 RPM. The cells grew to saturation overnight and were further diluted 1:5 and allowed to grow for ~4.5 hours, at which point they were diluted 1:10 and allowed to grow for 4.5 more hours.

4. The cells were removed from the shaker and diluted 1:50.

5. A 2µL aliquot was spotted onto a 3% agarose pad made of the growth medium, allowed to dry, and mounted on a glass-bottom dish.

6. After mounting, the sample dish was affixed to the microscope using double stick tape. Exposures were as follows:
    - Brightfield - 100ms, gain 4, 12bit
    - YFP - 250ms, gain 4, 12bit

7. Between five and ten positions were marked and were imaged every 7 minutes for 8 hours.

8. The samples were discarded and the dataset was transferred to the storage server.