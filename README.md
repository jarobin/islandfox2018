## Scripts for simulations in "Purging of strongly deleterious mutations explains long-term persistence and absence of inbreeding depression in island foxes" by Robinson et al. (2018).



**_SimCoalescent_SanNicolasPeaks_msprime_**

Scripts and files for coalescent simulations under putative models of San Nicolas island fox demographic history in msprime (Kelleher et al., 2016).

- `sim_SNIgenomes_msprime.py` runs the simulations in msprime.
- `SNI_ABC_posteriors_100.txt` contains parameters in the 100 best fitting San Nicolas demographic models from Robinson et al., 2016. The columns are Ne4, T, Ne3, Ne1.
- `Wong2010_dogRecRates_perbp.txt` contains per-site recombination rates in 1 Mb chunks of the dog genome, from Wong et al., 2010.



**_SimForward_IslandMainlandVariation_SLiM_**

Scripts for forward simulations of neutral and deleterious variation in island and mainland populations in SLiM (Haller and Messer, 2016).

- `run_slim_del_var.sh` executes the `make_slim_[model].job.sh` script given a set of parameters and a specific model, and submits the SLiM job.
- `make_slim_[model].job.sh` are scripts executed from within `run_slim_del_var.sh` to generate SLiM job scripts according to the specified parameter values.
