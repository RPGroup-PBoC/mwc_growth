---
status: Rejected 
reason: experiment not yet completed
---

# 2019-06-05 37C Glucose O2 Plate Reader Growth Measurement

## Purpose
This experiment aims to measure the growth rate of the *E. coli* strains of interest in M9 + 0.5% glucose with antibiotics at 37°C.

## Strain Information

| Location | Plasmid | Genotype | Host Strain | Shorthand |
| :------- | :------ | :------- | ----------: | --------: |
| dilution 1 - pos. 9 | `pZS3*PN25-tetR`| `galK<>25O2+11-YFP, gspI<>4*5noO1v1-CFP` |  HG105 |`deltaLacI` |

## Antibiotic Titration Series

| Antibiotic | Concentration |
| :------ | ------------: |
| Rifampicin | XX [ng / mL] |
| Apromycin | 0, 40, 100, 200, 400, 100, 2000, 4000, 1000, 20000, 40000 [ng / mL] |

## Notes & Observations
* Cells were harvested at OD_600nm .
* Rifamycin and streptomycin were titrated up to lethal concentrations (50µg/mL and 100µg/mL respectively).

## Analysis Files

**Whole Plate Growth Curves**
![plate layout](output/delta_glucoseA/gp_output_curves.png)

**Per Well Growth Rate Heatmap**
[![growth curves](output/delta_glucoseA/per_well_doubling_times_heatmap.png)]

## Experimental Protocol

1. Cells as described in "Strain Information" were grown to saturation in 1mL of LB Miller.

2. Cells were diluted 1:1000 into 50mL of growth media in a 250mL flask and grown to steady state.

3. 300µL of growth media were added to the first and last rows and columns of a square-welled, clear-bottomed 96 well plate, the total capacity of which was 630µL. Cells were diluted 100x and 300µL were added to the remaining wells.

4. Antibiotics were diluted in the arrangement depicted in 'output/growth_plate_layout.png'. Starting with column numbered 2 on the plate, the antibiotics were diluted 100x-1000x into the plate to achieve the concentrations listed in the "Antibiotic Titration Series" section, from the first non-zero concentration in column 2 and increasing to the right.

5. The plate was placed in a Biotek Gen5 plate reader and grown at 37C, shaking in a linear mode at the fastest speed. Measurements were taken every 7 minutes for approximately 18 hours.