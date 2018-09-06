---
status: 
reason: 
    
---

# 2018-06-21 Growth Measurement (Run 1)


## Purpose
This experiment aims to measure the growth rate of the *E. coli*
strains of interest in M9 + 0.5% glycerol at 37°C. This file corresponds to the
first run, although both experimens were conducted concurrently.


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
* Inoculum OD<sub>600nm</sub> was ~ 0.83 and was diluted 1:10 into the growth medium.
* The run was started at 10h20m.

## Analysis Files

![](output/20180621_r1_37C_glycerol_O2_growth.png)

## Experimental Protocol

1. Cells from a 4.5-hour saturated LB culture were diluted 1:1000 into 6 mL of M9 + 0.5% glycerol in a 15 mL Falcon tube. The inoculum was stored at room temperature for several hours prior to the inoculation.

2. Cells grew overnight (13h) at 37°C shaking at 225 RPM.

3. The OD<sub>600nm</sub> of these cells were measured and 5mL were diluted into 50 mL of prewarmed M9 + 0.5% glycerol and was thoroughly mixed.

4. Aliquots of 500 µL were serially pipetted to fill one 96-well plates (2 mL maximal volume). This plate was covered with a breathable nylon cover.

5. These cultures were allowed to grow at 37°C while shaking at 225 RPM. Growth measurements were taken every half hourly unless otherwise noted.

6. Growth measurements were performed by puncturing the nylon cover with a pipette and withdrawing 500 µL of culture. This aliquot was then transferred to a plastic cuvette (2 mL maximal volume) and measured using a table-top spectrophotometer. At each time point, an M9 blank was measured and used as a reference, followed by three measurements.

## Analysis Protocol

1. The elapsed time was calculated manually in the spreadsheet and was exported
to a `.csv` file.

2. The data was read and was trimmed only to the region of absorbance between
0.1 and 0.8 a.u. I defined exponential growth to be in this region.

3. The initial absorbance (A<sub>0</sub>), growth constant (λ), and likelihood
scale parameter (γ from a Cauchy likelihood) were inferred via MCMC.
