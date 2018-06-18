---
status: questionable
reason: 
    While the growth does appear to be exponential, the  doubling time
    estimate is still far above what we would expect for this growth medium.
---

# 2018-06-16 Growth Measurement


## Purpose
This experiment aims to measure the growth rate of the *E. coli*
strains of interest in M9 + 0.5% glucose at 37°C.

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
* Inoculum OD<sub>600nm</sub> was ~ 0.68 and was diluted 1:10 into the growth medium.
* The run was started at 10h50m.

## Analysis Files

![](output/20180221_r2_37C_glucose_O2_growth.png)

## Experimental Protocol

1. Cells from a 4-hour nearly saturated LB culture were diluted 1:10000 into 6 mL of M9 + 0.5% glucose in a 15 mL Falcon tube. The inoculum was stored at room temperature for several hours prior to the inoculation.

2. Cells grew overnight (14.5h) at 37°C shaking at 225 RPM.

3. The OD<sub>600nm</sub> of these cells were measured and were diluted into 50 mL of prewarmed M9 + 0.5% glucose and was thoroughly mixed.

4. Aliquots of 500 µL were serially pipetted to fill one 96-well plates (2 mL maximal volume). This plate was covered with a breathable nylon cover.

5. These cultures were allowed to grow at 37°C while shaking at 225 RPM. Growth measurements were taken hourly unless otherwise noted.

6. Growth measurements were performed by puncturing the nylon cover with a pipette and withdrawing 500 µL of culture. This aliquot was then transferred to a plastic cuvette (2 mL maximal volume) and measured using a table-top spectrophotometer. At each time point, an M9 blank was measured and used as a reference, followed by three measurements.

## Analysis Protocol

1. The elapsed time was calculated manually in the spreadsheet and was exported
to a `.csv` file.

2. The data was read and was trimmed only to the region of absorbance between
0.1 and 0.8 a.u. I defined exponential growth to be in this region.

3. The initial absorbance (A<sub>0</sub>), growth constant (λ), and likelihood
scale parameter (γ from a Cauchy likelihood) were inferred via MCMC.
