---
status: Rejected 
reason: experiment not yet completed
---

# 2018-11-07 37C Acetate O2 Microscopy Growth Measurement

## Purpose
This experiment aims to measure the growth rate of the *E. coli* strains of interest in acetate at 37°C.

## Strain Information
| Location | Plasmid | Genotype | Host Strain | Shorthand |
| :------- | :------ | :------- | ----------: | --------: |
| dilution 1 - pos. 9 | `pZS3*PN25-tetR`| `galK<>25O2+11-YFP, gspI<>4*5noO1v1-CFP` |  HG105 |`deltaLacI` |

## Titration Series

| Inducer | Concentration |
| :------ | ------------: |
| Anhydrotetracycline HCl (ATC) | 0 |
| Isopropylthiogalactopyranoside (IPTG) | 0|

## Notes & Observations
* `deltaLacI` grew very slowly.

## Analysis Files

**Growth Movie**
![growth movie](output/growth_movie.gif)

**Growth Curves and Rate Determination**
[![growth curve](output/growth_rate.png)]

## Experimental Protocol

1. Cells as described in "Strain Information" were grown to saturation in 3mL of LB Miller.

2. Cells were diluted 1:1000 into 3mL of growth medium in a 14mL Falcon tube.

3. Cells were allowed to grow for 32 hours at 37°C with shaking at ~ 220 RPM.

4. The cells were removed from the shaker and diluted 1:50.

5. A 2µL aliquot was spotted onto a 3% agarose pad made of the growth medium, allowed to dry, and mounted on a glass-bottom dish.

6. After mounting, the sample dish was affixed to the microscope using double stick tape. Exposures were as follows:
    - Brightfield - 100ms, gain 4, 12bit
    - YFP - 250ms, gain 4, 12bit

7. One or two positions were marked and were imaged every 10 minutes for 10 hours.

8. The samples were discarded and the dataset was transferred to the storage server.