# Results

## Compare strategies
One of our first objective is to reproduce the results of ([Mondain-Monval et al., 2024](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14355)). Specifically, those of figure 3(a):

![figure3a mondain-monval et al., 2024](images/figure3_mondain-monval_2024.png)

In their work, ([Mondain-Monval et al., 2024](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14355)) use a dataset (simulated data) that is quite different to ours (STOC data). Before comparing our results, let us remind how much our methods differs:

| Parameter             | M.M.                                      | Ours                      |
| :--------:            | :-------------:                           | :-----:                   |
| Data type             | Simulated with `virtualspecies` package   | STOC annual sampling      |
| Number of species     | 50 species                                | 5 species                 |
| Number of variables   | 10 (randomly selected from 33 available)  | 6 (direct from STOC)      |
| Model type            | Ensemble model (GLM + GAM + RF)           | HMSC                      |
| Re-training strategy  | Re-sampling of 1%, 10% or 50% of dataset  | Addition of new samples   |
| Training size         | 2000 cells                                | 125 cells + new data      |


With this context in mind, we can now see our results (note: GF = fap-filling and SU = simplified-uncertainty): [click to see the file](../outputs/results/STOC-OED/dotwhisker_strategies_comparison.pdf)



**What can we see?** 
- Mondain-Monval show Delta MSE ranging up to 0.008, with means around 0.0015.
- Our results are consistent with theirs, with slightly better Delta MSE
- In our results, SU performs better than GF (significance not tested), but with a thinner difference than Mondain-Monval.

**What to conclude?** 

Given the methodoly used, our results are pretty comparable to those of Mondain-Monval. This is exciting because we are using a different model (HMSC) and a completely different dataset.

