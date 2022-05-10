## Stochastic simulation

These scripts (for Julia 1.6.1) perform stochastic simulation of the extended [Subbalakshmi *et al.* 2021](https://www.karger.com/Article/FullText/512520) models:

* With MMI2: `sim_jolly_mmi_cdh1.jl`
* With MMI2 and miR-101: `sim_jolly_mmis_cdh1.jl`

The model definitions `emt_*.jl` were generated from the deterministic Antimony models in Python scripts using `sbjulia.py`.
Example parameter sets and initial conditions are in the `parameterizations` directory.
