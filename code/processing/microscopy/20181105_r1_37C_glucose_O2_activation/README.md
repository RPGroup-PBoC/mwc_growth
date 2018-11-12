---
status: Rejected 
reason: experiment not yet completed
---

# 2018-11-05 37C Glucose O2 Activation Snaps

## Purpose
This experiment was the first view of YFP and xapR-mCherry expression following integration of these genes for the activation construct.

## Strain Information

| Location | Plasmid | Genotype | Host Strain | Shorthand |
| :------- | :------ | :------- | ----------: | --------: |
| activation - pos. X |  |  |  | `autofluorescence` |
| activation - pos. X |  |  |  | `xapR-mCherry`|
| activation - pos. X |  |  |   |`28-yfp` |

## Titration Series

| Inducer | Concentration |
| :------ | ------------: |
| Anhydrotetracycline HCl (ATC) | None |
| Isopropylthiogalactopyranoside (IPTG) | None |

## Notes & Observations
* Cells were harvested at OD<sub>600nm</sub> of ~ 0.3.
* Cells waited on the dish for ~ 1.5 hours between plating them and imaging them, due to unexpected issues with the microscope.
* mCherry was very bright in the `xapR-mCherry` sample, and seemed to be aggregating into puncta at the poles.
* YFP was dim in the `28-yfp` sample, didn't seem to be expressing much, if at all.
 
## Analysis Files

**Calibration Factor Determination**
[![dilution summary](output/dilution_summary.png)](output/dilution_summary.html)

**Fold-change**
[![fold-change summary](output/foldchange_summary.png)](output/foldchange_summary.html)

## Experimental Protocol

1. Cells as described in "Strain Information" were grown to saturation in 3mL of LB Miller.

2. Cells were diluted 1:1000 into 3mL of growth media in 14mL Falcon tubes.

3. Tubes were placed in a rack and allowed to grow for 6 hours at 37°C with shaking at ~ 220 RPM.

4. Once the cells reached an OD<sub>600nm</sub> between ~ 0.3, the cells were removed from the warm room and harvested.

**Microscopy**

1. The samples were diluted 1:10 into growth medium. Aliquots of 1µL were added to agarose pads made of the growth medium.

2. Agarose pads spotted with cells were allowed to dry and were then placed onto a glass bottom dish.

3. After mounting, the sample dish was affixed to the microscope using double stick tape. Approximately 5 positions were marked per snapshot sample. Exposures were as follows:
    - Brightfield - 100ms, gain 4, 12bit
    - mCherry - 5000ms, gain 1, 12bit
    - YFP - 5000ms, gain 1, 12bit

4. The samples were discarded and the dataset was transferred to the storage server.