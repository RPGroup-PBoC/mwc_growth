---
status: Accepted
reason: Need to examine the growth curves more carefully. 
---

# 2018-10-24 37C Glucose O2 Microscopy Growth Measurement

## Purpose
This experiment was the flagship growth experiment in M9 + 0.4% glucose.

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
* Cells were harvested at an OD_600nm_ of ~0.45. 

## Analysis Files

**Growth Movie**
![growth movie](output/growth_movie.gif)

**Growth Curves and Rate Determination**
[![growth curve](output/growth_rate.png)]

## Experimental Protocol

1. Cells as described in "Strain Information" were grown to saturation in 3mL of LB Miller.

2. Cells were diluted 1:1000 into 3mL of growth medium in a 14mL Falcon tube.

3. Cells were allowed to grow for 9 hours at 37°C with shaking at ~ 220 RPM. They were then diluted 1:2 and allowed to grow for ~1 hour.

4. The cells were removed from the shaker and diluted 1:50.

5. A 2µL aliquot was spotted onto a 3% agarose pad made of the growth medium, allowed to dry, and mounted on a glass-bottom dish.

6. After mounting, the sample dish was affixed to the microscope using double stick tape. Exposures were as follows:
    - Brightfield - 100ms, gain 4, 12bit
    - YFP - 250ms, gain 4, 12bit

7. One or two positions were marked and were imaged every 5 minutes for 6 hours.

8. The samples were discarded and the dataset was transferred to the storage server.