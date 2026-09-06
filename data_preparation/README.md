# Data preparation

`build_elsa_longitudinal_analysis.R` builds `elsa_longitudinal_analysis.csv`
from the raw ELSA wave files (Waves 2-7 core data and derived-variables
files, in SPSS `.sav` format as distributed by the UK Data Service). It
constructs the loneliness, depression, transport-mobility, and mobility-
limitations variables used throughout the paper and reshapes the panel to
long (person-wave) format.

Portions of this preprocessing are adapted from Mayerl, H., Stolz, E., &
Freidl, W. (2023). *Lonely and depressed in older age: prospective
associations and common vulnerabilities.* Aging & Mental Health, 27(3),
640-645 (original code: https://osf.io/jrhq7/overview).

## Usage

1. Obtain ELSA Waves 2-7 access via the UK Data Service (see the root
   README's Data Availability section) and place the raw `.sav` files in
   this folder.
2. Run this script; it writes `elsa_longitudinal_analysis.csv`.
3. Copy that file into both `../BJCM/data/` and `../BMEOP/data/`.
