---
status: Accepted
reason: Looks great.
---

# 2019-02-04 37C Mixed Media Multistrain Plate Reader Growth Measurement

## Purpose
This experiment aims to measure the growth rate of the *E. coli* strains of interest in M9 + 0.5% glucose, M9 + 0.5% glycerol, and M9 + 0.5% acetate at 37°C.

## Notes & Observations
* Saturated LB cultures of `deltaLacI`, `ManuelDelta`, `GenevaLE`, and `GenevaBW` were diluted 1000x in M9 + 0.5% glucose about 4.5 hours prior to start of measurement.
* A saturated LB culture of `deltaLacI` was diluted 1000x into M9 + 0.5% glycerol and M9 + 0.5% acetate about 4.5 hours prior to start of measurement.
* The four glucose cultures were diluted 10x into the plate (from OD_600nm ~0.08-0.12). The glycerol and acetate cultures had OD_600nm <0.05 and were not diluted into the plate.
* `GenevaLE` and `GenevaBW` LB cultures (donated by the Van Valen Lab on the morning of 2019-01-29) were stored at 4˚C until dilution.

## Analysis Files

**Whole Plate Growth Curves**
![plate layout](output/delta_glucose/gp_output_curves.png)

**Per Well Growth Rate Heatmap**
[![growth curves](output/delta_glucose/per_well_doubling_times_heatmap.png)]

## Experimental Protocol

1. Cultures of `deltaLacI` and `ManuelDelta` were grown to saturation in 3mL of LB Miller. Saturated LB cultures of the strains `GenevaBW` and `GenevaLE` were supplied by Geneva in the Van Valen lab.

2. Saturated cultures of `deltaLacI`, `ManuelDelta`, `GenevaLE`, and `GenevaBW` were diluted 1000x in M9 + 0.5% glucose in 14mL Falcon tubes and allowed to grow for about 4.5 hours. The saturated culture of `deltaLacI` was also diluted 1000x into M9 + 0.5% glycerol and M9 + 0.5% acetate at the same time.

3. Cells were removed from the shaker and the four glucose cultures were diluted 1:10.

4. 100µL of water were added to the first and last row, and the first and last two columns of a round-welled, clear-bottomed 96 well plate, the total capacity of which was 250µL. 100µL of each culture were added to the remaining wells in the arrangement depicted in the file 'output/growth_plate_layout.png'.

5. Double-sided tape was placed between the edges of the plate and the lid to prevent them rubbing and creating plastic dust.

6. The plate was placed in a Biotek Gen5 plate reader and grown at 37C, shaking in a linear mode at the fastest speed. Measurements were taken every 7 minutes for approximately 44 hours.