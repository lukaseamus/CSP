# Climate-driven shifts in kelp forest composition reduce carbon sequestration potential
This repository contains data and R code accompanying article 10.1111/gcb.16299 in *Global Change Biology*.

The repository is split into five folders. **Assimilation**, **Export**, **Decomposition** and **Sequestration** contain all files on carbon assimilation, export, remineralisation and potential sequestration. **Irradiance** contains files required to model the local light regime by depth, which is required in **Sequestration** to estimate cumulative detrital carbon assimilation. Below is a description of each file within those folders.

**Analysis**
1. `Environmental.csv`|`Environmental.R`: *In situ* light and temperature data.
    - *d* = day from start of experiment as an integer
    - *day* = day from start of experiment as a numeric
    - *lux* = irradiance in lumens per square metre
2. `Decomposition.csv`|`Decomposition.R`: Biomass loss data.
    - *species* = *Laminaria digitata* (d), *Laminaria hyperborea* (h) or *Laminaria ochroleuca* (o)
    - *bag* = ID of galvanised steel mesh bag
    - *age* = detrital age in days
    - *loss* = biomass loss per day (%)
3. `Extraction.csv`|`Extraction.R`: Comparison of pigment extraction from fresh and lyophilised tissue. All pigment concentrations are given in micrograms per gram.
    - *id* = sample ID
    - *pigment* = Chlorophll *a*, Chlorophyll *c* or Carotenoids
    - *wet* = fresh extract pigment concentration
    - *dry* = lyophilised extract pigment concentration
4. `Pigments.csv`|`Pigments.R`: Main pigment data. All pigment concentrations are given in micrograms per gram.
    - *id* = sample ID
    - *site* = Mount Batten (MB) or West Hoe (WH)
    - *species* = *Laminaria digitata* (d), *Laminaria hyperborea* (h) or *Laminaria ochroleuca* (o)
    - *bag* = ID of plant (P) or galvanised steel mesh bag (B)
    - *age* = detrital age in days
    - *bb.Car* = β,β-carotene
    - *Chl.a* = Chlorophll *a*
    - *Chl.c1* = Chlorophll *c*<sub>1</sub>
    - *Chl.c2* = Chlorophll *c*<sub>2</sub>
    - *Fuco* = Fucoxanthin
    - *Phe.a* = Pheophytin *a*
    - *Viola* = Violaxanthin
    - *Zea* = Zeaxanthin
5. `Map.R`: R code to plot the base map of Europe used to build Figure 1.

**Deconvolution**
1. `gaussian.peak.parameters.txt`: Gaussian peak parameters for 28 common algal pigments that determine the shape of the absorbance curve for each individual pigment (this file was left unchanged so see 10.1371/journal.pone.0137645 for details).
2. `specific.absorption.coefficients.txt`: Customised set of core pigments to be deconvoluted from *Laminaria* spp. pigment mixture absorbance spectra (details were obtained from the literature for acetone extractions since 10.1371/journal.pone.0137645 used ethanol)
    - *pigment* = pigment abbreviation as used in 10.1371/journal.pone.0137645
    - *nm* = absorbance maximum or the wavelength at which absorbance is highest given in nanometres
    - *L.g.cm* = absorption coefficient (sometimes also called the specific extinction coefficient) at the absorbance maximum given in litres per gram per centimetre
    - *solvent* = only values for acentone were used since this was the solvent we extracted our samples with
3. `pigment.function.R`: List of R functions required for spectral deconvolution. This file requires the above datasets.
4. `Wet.csv`|`wet.mass.csv`|`Wet.R`: Customised R code and raw spectrophotometric data from fresh tissue extractions required to obtain estimates of individual pigment concentrations.
    - *wavelength* = wavelength given in nanometres
    - *d1A* etc. = sample ID
    - *id* = sample ID
    - *g* = fresh mass of each sample in grams
5. `Dry.csv`|`Dry.R`: Same as above but for lyophilised tissue extraction. There is no equivalent of `wet.mass.csv` since sample weight was standardised to 100 ± 1 milligrams.
6. `Spectrum.R`: R code to visualise the deconvolution procedure on the basis of a single lyophilised *Laminaria digitata* sample (Figure S2). 

Luka Seamus Wright, 8 July 2022
