---
status: Rejected 
reason: experiment not yet completed 
---

# YEAR-MONTH-DATE TEMP CARBON OPERATOR Absorbance Growth Measurement (Run #)

## Purpose
This experiment aims to measure the growth rate of the *E. coli* strains of interest in XX at 37°C.
This file corresponds to the XX run, although both experiments were conducted concurrently.

## Strain information
| Location | Plasmid | Genotype | Host Strain | Shorthand |
| :------  | :------ | :------- | ----------: | --------: |
| dilution 1 - pos. 9 | `pZS3*PN25-tetR`| `galK<>25O2+11-YFP, gspI<>4*5noO1v1-CFP` |  HG105 |`deltaLacI` |

## Titration Series

| Inducer | Concentration |
| :-----  | ------------: |
| Anhydrotetracycline HCl (ATC) | 0 ng / mL |
| Isopropylthiogalactopyranoside (IPTG) | 0 mM |

## Notes & Observations
* Inoculum OD<sub>600nm</sub> was ~ X and was diluted 1:X into the growth medium.
* The run was started at XhXm.

## Analysis Files

![growth curve](output/growth_rate.png)

## Experimental Protocol

1. Cells as described in "Strain Information" from a X-hour saturated LB culture were diluted 1:X into 3 mL of growth medium in a 14 mL Falcon tube.

2. Cells grew overnight (Xh) at 37°C shaking at ~ 220 RPM.

3. The OD<sub>600nm</sub> of these cells were measured and X µL were diluted into 3 mL of prewarmed growth medium in a 14 mL Falcon tube and was thoroughly mixed.

4. This culture was allowed to grow at 37°C while shaking at ~ 220 RPM. Growth measurements were taken every half hourly unless otherwise noted.

5. Growth measurements were performed by directly measuring the 14 mL Falcon tube using a table-top spectrophotometer. At each time point, a blank was measured and used as a reference and tubes were held stationary for at least 30 seconds, followed by three measurements.

## Analysis Protocol

1. The elapsed time was calculated manually in a spreadsheet and was exported to a `.csv` file.

2. The data was read and was trimmed only to the region of absorbance between 0.1 and 0.8 a.u. I defined exponential growth to be in this region.

3. The initial absorbance (A<sub>0</sub>), growth constant (λ), and likelihood scale parameter (γ from a Cauchy likelihood) were inferred via MCMC.
