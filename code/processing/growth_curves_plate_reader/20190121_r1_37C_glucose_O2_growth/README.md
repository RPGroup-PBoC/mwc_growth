---
status: Rejected 
reason: experiment not yet completed
---

# 2019-01-21 37C Glucose Multistrain Plate Reader Growth Measurement

## Purpose
This experiment aims to measure the growth rate of the *E. coli* strains of interest at 37°C.

## Notes & Observations
* `deltaLacI` (old) and `ManuelDelta` (potentially `ManuelDilution` box 1, pos 49 instead) were grown to saturation in LB + spec, then diluted 1000x in M9 + 0.5% glucose about 4.25 hours prior to start of measurement.
* Both samples were diluted 2x into the plate.
* Plastic dust was observed on the plate at the completion of the experiment, as well as noise in the growth curves of the center wells.

## Analysis Files

**Whole Plate Growth Curves**
![plate layout](output/delta_glucose/gp_output_curves.png)

**Per Well Growth Rate Heatmap**
[![growth curves](output/delta_glucose/per_well_doubling_times_heatmap.png)]

## Experimental Protocol

1. Cultures of `deltaLacI` and `ManuelDelta` were grown to saturation in 3mL of LB Miller + spec.

2. The cells were diluted 1000x into M9 + 0.5% glucose in 14mL Falcon tubes and allowed to grow for about 4.25 hours.

3. Cells were removed from the shaker and both samples were diluted 1:2.

4. 100µL of water were added to the first and last row and column of a round-welled, clear-bottomed 96 well plate, the total capacity of which was 250µL. 100µL of each culture were added to the remaining wells in the arrangement depicted in the output file 'growth_plate_layout.png'.

5. The plate was placed in a (Murray Lab) Biotek plate reader and grown at 37C, shaking in a linear mode at the fastest speed. Measurements were taken every 7 minutes for approximately 12 hours.