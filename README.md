# malaria-angumu-2024

Analysis code for a repeated cross-sectional evaluation of mass drug administration (MDA) combined with indoor residual spraying (IRS) for malaria control in Angumu health zone, Ituri province, Democratic Republic of the Congo (five population-based surveys, March 2022 to May 2024; Epicentre and MSF Operational Centre Geneva, with the Programme National de Lutte contre le Paludisme). Health areas were compared across three zones (Routine, IRS, MDA+IRS) for malaria prevalence, under-five and crude mortality, with entomological and routine facility surveillance as supporting data.

The scripts reproduce the tables, figures and in-text results of the manuscript and its supplementary files, in the order they appear:

| Script | Reproduces |
|---|---|
| `code/01_study_design.R` | Figure 1 (study area), Figure 2 (timeline and zone classification), S1 Table (baseline characteristics), S4 Table (clusters and households), S7 Table (participants by survey), S8 Table (MDA and IRS coverage) |
| `code/02_prevalence.R` | Table 1 (malaria prevalence by zone), Figure 3 (prevalence trend), Figure 4 (prevalence maps), S6 Table (Gengere), S9 Table and S4 Figure (prevalence by health area) |
| `code/03_prevalence_sensitivity.R` | S2 Table, panels A to G (ratio of ORs versus survey 1, Awasi-only comparator, health-area random intercept, area-level permutation test, MDA+IRS compared with IRS, E-values, net-adjusted prevalence), S12 Table (bed-net use) |
| `code/04_mortality.R` | Table 2 (under-five and crude mortality), S2 Figure (crude mortality rates), S3 Table, panels A to E, S10 Table (household-reported cause of death) |
| `code/05_entomology.R` | Table 3 (entomological indicators), Figure 5 (human biting rate by health area), S3 Figure (monthly entomological indicators), S5 Table (Additional file 2) |
| `code/06_facility_surveillance.R` | Figure 6 (facility test positivity rate against survey prevalence), S1 Figure (monthly facility test positivity rate), S11 Table |

`code/functions.R` holds helpers shared by the scripts. `code/prep.R` documents how the input files are built from the study team's internal data store; it cannot be run outside that store.

## Running

Run the scripts from the repository root, with R 4.5 or later, the packages loaded at the top of each script installed, and the input files in `data/` (see below). Outputs are written to `output/tables/` and `output/figures/`.

```
Rscript code/01_study_design.R
Rscript code/02_prevalence.R
Rscript code/03_prevalence_sensitivity.R
Rscript code/04_mortality.R
Rscript code/05_entomology.R
Rscript code/06_facility_surveillance.R
```

`03_prevalence_sensitivity.R` reads the model saved by `02_prevalence.R`, which takes about 50 minutes. For the net-adjusted prevalence models (S2 Table G) and the all-age mortality models, the optimiser was chosen by an `lme4::allFit()` search; the scripts fit the chosen optimiser directly. Figure 2 downloads the Inter font from Google Fonts.

## Inputs

Not included (see Data access below):

- `data/df_ind.rds`: one row per listed household member, surveys 1 to 5, with survey, cluster, household and health area codes, zone, age, sex, displacement status, survey weight, malaria rapid diagnostic test result, bed-net use, MDA and IRS receipt, death, cause of death and person-time (surveys 1, 3 and 5), and household location.
- `data/df_hh.rds`: one row per visited household, with consent and interview date.
- `data/entomology_trapnights.rds`, `data/entomology_mosquitoes.rds`: CDC light-trap nights and the Anopheles caught, with species and sporozoite results.
- `data/pnlp_health_area_monthly.rds`, `data/dhis2_monthly.rds`: monthly routine malaria indicators by health area (PNLP and DHIS2).
- `data/climate_weekly.rds`: weekly rainfall and temperature for the health zone (PNLP).
- `data/maps/health_areas.rds`: MSF health-area boundaries against which survey clusters were sampled, used for Figure 4.

Included in `data/maps/`:

| File | Source |
|---|---|
| `grid3_cod_health_areas_v9_angumu.geojson` | Health-area boundaries, GRID3 DRC health areas v9.0 (CC BY 4.0) |
| `drc_provinces_grid3_v9.gpkg` | DRC provinces, dissolved from GRID3 DRC health zones v9.0 (CC BY 4.0) |
| `geoBoundaries-COD-ADM0_simplified.geojson` | DRC national boundary, geoBoundaries (CC BY 4.0) |
| `ne_10m_lakes.rds` | Lakes, Natural Earth (public domain) |
| `lake_albert.rds` | Lake Albert, OpenStreetMap contributors (ODbL) |

## Data access

Data used in this study are subject to restrictions arising from participant informed consent and privacy protections for identifying and potentially stigmatising health information, and cannot be deposited in a public repository. Anonymised data are available on request to the [Epicentre Data Sharing Committee](https://epicentre.msf.org/en/ongoing-projects/study-data-access-request), which reviews all requests independently of the study authors in accordance with [Epicentre's Data Sharing Policy](https://epicentre.msf.org/sites/default/files/2020-02/Epicentre-Data%20Sharing%20Policy-EN-2020.pdf).

## Citation

If you use this code, please cite the manuscript (citation details will be updated on publication) and the code itself.
