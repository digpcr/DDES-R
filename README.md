# DDES

The `DDES` R package converts digital polymerase chain reaction (dPCR) 
instrument data to the Digital PCR Data Essentials Standard (DDES) format.

For more information on the data standard, please refer to the preprint: 
https://www.biorxiv.org/content/10.64898/2026.04.10.717368v1

## Installing the package

If you do not yet have `pak` installed, install it first from CRAN:

```r
install.packages("pak")
```

Installing the `DDES` package directly from GitHub is then achieved with `pak`:

```r
pak::pak("digpcr/DDES-R")
```

## Usage

Usage of the package is detailed in the 
[QIAcuity](vignettes/qiacuity-to-ddes.Rmd) and
[QX200](vignettes/qx200-to-ddes.Rmd) vignettes.

## Version history

### **1.1.0** 
Added conversion capability for Bio-Rad QX200 QuantaSoft exports.

### **1.0.0** 
Original conversion package for QIAGEN QIAcuity exports.
