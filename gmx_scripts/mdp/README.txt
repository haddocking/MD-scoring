### This folder contains all .mdp files used for the preparations of the systems (protein-protein and protein-peptide) for MDruns.

### .mdp files for preparation (#01-05) originate from the MolMod tutorial(http://www.bonvinlab.org/education/molmod/simulation/) and are adapted to run two seeds.

### Unused file:
- 06_md_PME.mdp

### Template files:
- template_03_nvt_pr1000_PME.mdp
### for NVT conditions, but added thermalization steps of 50 K, 150 K and 300 K to let the system slowly adjust to the new conditions.
- 05_npt_NOpr_PME.mdp
- template_md_prod_.mdp
### Template file for the actual production of the simulation

### Files for system preparation:
- 01_em_vac_PME.mdp
- 02_em_sol_PME.mdp
- 03_nvt_pr1000_t300_PME.mdp
- 03_nvt_pr1000_t340_PME.mdp *
*for runs at 340 K (PRE5-PUP2, Barnase-barstar, E2A-HPr, Capri 11 and Capri 15
- 03_nvt_pr1000_t50_PME.mdp
- 03_nvt_pr1000_t150_PME.mdp
- 04_npt_pr_t300_PME.mdp
- 04_npt_pr_t340_PME.mdp
- 05_npt_NOpr_PME_t300.mdp
- 05_npt_NOpr_PME_t340.mdp
- md_prod_rosanna_t300.mdp
- md_prod_rosanna_t340.mdp

### Files for preparation of the second seed:
- 03_nvt_pr1000_t150_PME-secondseed.mdp
- 03_nvt_pr1000_t300_PME-secondseed.mdp
- 03_nvt_pr1000_t50_PME-secondseed.mdp
- md_prod_rosanna_t300-secondseed.mdp

### Files for preparation of the system taking interface A and B (and chain A and chain B) into account for writing energy terms to the .edr file during the production simulation run:
- md_prod_rosanna_t300_energygrps-chainA-B.mdp

### And the second seeds:
- md_prod_rosanna_t300_energygrps-chainA-B-secondseed.mdp
