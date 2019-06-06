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
| Rifamycin | 0, 5, 10, 50, 100, 500, 1000, 10000, 50000 [ng / mL] |
| Streptomycin | 0, 10, 20, 100, 200, 1000, 2000, 20000, 100000 [ng / mL] |

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

3. 300µL of growth media were added to the first and last rows and first two and last two columns of a square-welled, clear-bottomed 96 well plate, the total capacity of which was 630µL. 300µL of cells were added to the remaining wells, in the arrangement depicted in 'output/growth_plate_layout.png'.

4. The plate was placed in a Biotek Gen5 plate reader and grown at 37C, shaking in a linear mode at the fastest speed. Measurements were taken every 7 minutes for approximately 72 hours.