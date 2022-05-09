## Multistable parameter search

These scripts (for CPython 3.8.5) find example multistable parameter sets for the examined models:

* One-microRNA MMI4: `find_mmi4.py` (example parameter sets) and `count_mmi4.py` (count many parameterizations by number of attractors)
* Two-microRNA MMI4: `find_mmi4_22.py` (plausible example parameter sets) and `find_mmi4_22_uncorrelated.py` (same ranges as for one-microRNA model, for cooperativity distributions)
* Chained-MMI2: `find_chained_mmi2.py`
* Co-targeting-MMI2: `find_cotargeting_mmi2.py`
* Chained-MMI2 with transcriptional positive feedback loop: `find_tma.py` (find bistable transcriptional mutual activation parameters) and `find_combined_chained_mmi_tma.py` (search for many-attractor fusions)
* MMI2 with transcriptional repression: `bif_mmi2_mirep.py`
* Two-microRNA MMI4 with transcriptional repressions: `find_mmi4_22_mireps.py`
* [Subbalakshmi *et al.* 2021](https://www.karger.com/Article/FullText/512520) with MMI2: `find_jolly_mmi2.py` (find multistable parameters) and `plot_jolly_mmi2.py` (plot with E-cadherin readout)
* Subbalakshmi *et al.* 2021 with MMI2 and miR-101: `find_jolly_multimir.py` (find multistable parameters) and `plot_jolly_multimir.py` (plot with E-cadherin readout)

Parameterizations of other models can be plotted with `plotsystems.py`.
