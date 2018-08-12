# The `buildmer` package

`buildmer` is an (experimental) `R` package written to simplify the process of testing whether the terms in your `lmer` (or equivalent) models make a significant contribution to the *-2 Log Likelihood*, AIC criterion, or BIC criterion. The aim of the package is to **fully automate model selection**, as is already possible for non-mixed regression models and lmerTest models using the `step` function.

In addition, the package (optionally) **determines the order of your predictors** by their contribution to the model fit. This is intended to mimic the results you might get from the stepwise functions available in, e.g., *SPSS* (although *SPSS* support stepwise elimination only for fixed-effect models...). In sum, the `buildmer` package aims to take the complex and time-consuming parts of the model fitting procedure out of your hands -- all you need to do is specify your intended maximal model and your dataset, and the package will take care of the rest.

**Nonconvergence** of models is handled properly by removing a random slope if the maximal model you have specified turns out to not converge. **p-values** are calculated for you using Wald z-scores. Better alternatives (Kenward-Roger, Satterthwaite) are available at the user's option (if the `lmerTest` package is installed).

The package supports the fitting of a wide variety of models, if the relevant packages are available. The following `buildmer` functions make it possible to fit the following types of models:
 * *buildmer*: `lm`, `glm`, `lmer`, `glmer`, `gamm4` (package `gamm4`)
 * *buildgam*, *buildbam*: `gam`, `bam` (package `mgcv`)
 * *buildglmmTMB*: `glmmTMB` (package `glmmTMB`)
 * *buildjulia*: uses `RCall` to drive a Julia installation to fit models using Julia package `MixedModels`
No support is available (or planned) for fitting generalized additive mixed models via `mgcv`'s `gamm` function; use `gamm4`s eponymous function instead.

Automatic elimination of fixed, random, and/or smooth terms, is possible and enabled by default using the backward (default) or forward stepwise method. Bi-directional elimination is also possible, by passing e.g. `direction=c('forward','backward','forward')`, although I would not want to recommend doing this.

The intention of the `buildmer` package is to make your life as simple as:

```
library(buildmer)
model = buildmer(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)
```

In other words, the intention of the package is for you to specify the maximal model that you *would like* to fit, and let the package worry about whether or not your model converges, and whether or not your effect structure before, during, or after non-significant term removal needs to be passed to `gamm4`, `lmer`, or `lm` (or `glm`(`er`) or `{g,b}am`).

Many more options are available; please see the documentation for details. For users familiar with `SPSS`, a less-convoluted alternative to the `buildmer` command is available called `stepwise`. That command only requires three arguments: your model formula, your dataset, and (optionally) the error family you are fitting (leave empty if you're doing linear regression, specify `family=binomial` for logistic regression, etc).

## Installation

```
library(devtools)
install_github('cvoeten/buildmer')
```

## Known issues

This package is work-in-progress. **You should always check the output `buildmer` provides during the fitting process, and verify that it is making the correct decisions in terms of which terms it chooses to add or remove from the model**. If you notice a mistake, you can always correct it and instruct `buildmer` to leave the problematic parts alone, e.g.: `buildmer(...,reduce.random=FALSE)`. Please file a bug report if you ever need to do this!

Bugs will be fixed as they are uncovered. If you notice a bug, please file an issue for it and I will look into it if I have time.
