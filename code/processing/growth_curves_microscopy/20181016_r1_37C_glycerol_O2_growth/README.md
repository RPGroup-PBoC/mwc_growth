---
status: Accepted
reason: 
---

# 2018-10-18 Microscopy Growth Measurement

## Purpose
This experiment was performed to measure the growth rate of our cells on M9 + 0.4% glycerol.

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
* The sample was diluted 10x from a culture that had been overgrown in the morning and had been diluted 10x for an induction experiment.

## Analysis Files

**Growth Movie**
![](20181016_r1_37C_glycerol_O2_growth_xy7.gif)

**Growth Curves and Rate Determination**
![]()

## Experimental Protocol

1. Cells as described in "Strain Information" were grown to saturation in 3mL of LB Miller.

2. Cells were diluted 1:1000 into 3mL of growth medium in a 14mL Falcon tube.

3. Cells were allowed to grow for ~12 hours at 37°C with shaking at ~ 220 RPM. The cells grew to near saturation overnight (OD_600nm_ ~ 0.9) and were further diluted 1:10 and allowed to grow for ~9.5 hours, at which point they were diluted 1:10 again and allowed to grow for 3 more hours.

4. The cells were removed from the shaker and diluted 1:50.

5. A 1µL aliquot was spotted onto a 3% agarose pad made of the growth medium, allowed to dry, and mounted on a glass-bottom dish.

6. After mounting, the sample dish was affixed to the microscope using double stick tape. Exposures were as follows:
    - Brightfield - 100ms, gain 4, 12bit
    - YFP - 250ms, gain 4, 12bit

7. Between five and ten positions were marked and were imaged every 7 minutes for 8 hours.

8. The samples were discarded and the dataset was transferred to the storage server.