# Climate-driven shifts in kelp forest composition reduce carbon sequestration potential
This repository contains data and annotated R code accompanying article 10.1111/gcb.16299 in *Global Change Biology*.

The repository is split into five folders. **Assimilation**, **Export**, **Decomposition** and **Sequestration** contain all files on carbon assimilation, export, remineralisation and potential sequestration. **Irradiance** contains files required to model the local light regime by depth, a prerequisite for cumulative detrital carbon assimilation estimation in **Sequestration**. Below is a description of each file within those folders.

**Assimilation**
1. `Assimilation.csv`: Net and gross carbon assimilation data.
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **bag** = random factor (categorical variable) with levels plant (P) and mesh bag (B)
    - **age** = detrital age in days
    - **R** = respiration rate given in µmol oxygen per gram of buoyant mass per hour (buoyant mass is practically identical to wet mass)
    - **NPP** = net photosynthesis rate given in µmol oxygen per gram of buoyant mass per hour
    - **GPP** = gross photosynthesis rate given in µmol oxygen per gram of buoyant mass per hour (NPP + R)
    - **d:w** = dry to wet mass ratio
2. `Assimilation.R`: Code to analyse and visualise carbon assimilation.
    - **Input** = `Assimilation.csv`, Figure 5b from `Decomposition.R`
    - **Output** = Figure 5, carbon assimilation results

**Export**
1. `Export.csv`: Carbon export data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **time** = numerical expression of months
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **dw.export** = biomass export given in grams of dry mass per plant per day
    - **fw.export** = biomass export given in grams of wet mass per plant per day, converted from dry mass with plant-specific dry to wet mass ratios
    - **fw.export.avg** = biomass export given in grams of wet mass per plant per day, converted from dry mass with species- and month-specific dry to wet mass ratios
    - **C.export** = carbon export given in grams per plant per day, converted from dry mass with species- and month-specific carbon content (%)
2. `Carbon.csv`: Lamina carbon content data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **carbon** = carbon content (%)
3. `Mass.csv`: Sporophyte mass data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **mass** = whole plant wet mass given in grams
4. `DW.csv`: Dry to wet mass ratio data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **d.w** = dry to wet mass ratio
5. `Density.csv`: Sporophyte density data.
    - **month** = month and year given as MMM-YY
    - **season** = categorical variable with levels Spring, Summer, Autumn and Winter
    - **species** = categorical variable with levels *Laminaria digitata* (d), *Laminaria hyperborea* (h) and *Laminaria ochroleuca* (o)
    - **density** = number of plants per square metre. Note that *Laminaria hyperborea* and *Laminaria ochroleuca* share a quadrat while *Laminaria digitata* occurs spatially separated. 
6. `Export.R`: Code to analyse and visualise carbon export.
    - **Input** = `Export.csv`, `Carbon.csv`, `Mass.csv`, `DW.csv`, `Density.csv`
    - **Output** = Figure S2, `Constants.csv`, carbon export results

**Decomposition**

**Irradiance**

**Sequestration**


Luka Seamus Wright, 8 July 2022
