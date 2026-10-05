# denovogear-bb (development snapshot)

Early development code (2012–2013) for a **beta-binomial extension of
[DeNovoGear](https://github.com/ultimatesource/denovogear)**, built to call
mutations from paired samples (for example tumor vs. matched normal) in
whole-exome sequencing data.

> **Looking for the usable version?** This work was cleaned up and released as
> [**DNG_BetaBinomial_Paired**](https://github.com/gatoravi/DNG_BetaBinomial_Paired),
> which has install instructions, usage and an output description. This
> repository is kept as a record of the development history.

## The idea

Standard callers model allele read counts as binomial. In exome data, capture
and amplification biases add extra variance, so counts are overdispersed. This
code:

- fits beta-binomial α and β parameters to read counts across sites by
  maximum likelihood (`AlphaBeta.cc`, with the optimisation in R through
  RInside), and
- uses beta-binomial genotype likelihoods inside DeNovoGear's mutation model
  to score paired-sample calls (`BBPaired.cc`, `DNGCaller.cc`).

## Contents

| File | Purpose |
| --- | --- |
| `denovogear.cc` | Command-line entry point |
| `BBPaired.cc` | Paired-sample beta-binomial calling |
| `AlphaBeta.cc`, `rinside.BBFit.cc` | Estimation of beta-binomial parameters (calls R through RInside) |
| `DNGCaller.cc` | Genotype configuration likelihoods and posterior probabilities |
| `makeLookup.cc`, `lookup.h` | Precomputed mutation and genotype lookup tables |
| `likCalc.java` | Prototype likelihood calculation |
| `callDNGBB.sh` | Example invocation |
| `newmat/` | Bundled newmat matrix library |

Files ending in `~` or wrapped in `#` are editor backups.

## Author

Avinash Ramu, Conrad Lab, Washington University in St. Louis.
