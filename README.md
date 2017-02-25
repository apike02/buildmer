# The `buildmer` package

`buildmer` is an (experimental) `R` package written to simplify the process of testing whether the terms in your `lmer` (or equivalent) models make a significant contribution to the deviance (or *-2 Log Likelihood*). The aim of the package is **to fully automate model selection**, as is already possible for non-mixed regression models using the `step` function.

In addition, the package (optionally) **determines the order of your predictors** by their contribution to the deviance. This is intended to mimic the results you might get from the stepwise functions available in, e.g., *SPSS* (although they only support fixed-effect models only). In sum, the `buildmer` package aims to take the complex and time-consuming parts of the model fitting procedure out of your hands -- all you need to do is specify your intended maximal model and your dataset, and the package will take care of the rest.

**Nonconvergence** of models is handled properly by removing a random slope if the maximal model you have specified turns out to not converge. **p-values** are calculated for you using Wald z-scores. Better alternatives (Kenward-Roger, Satterthwaite) are available at the user's option (if the `lmerTest` package is installed).

The package supports the fitting of (`g`)lm, (`g`)`lmer`, and `gamm` (from the `gamm4` package) models, if the relevant packages are available. Forward and backward elimination of fixed and random effects are both possible. The intention is for your life to be as simple as:

```
library(buildmer)
model = buildmer(Reaction ~ Days + (Days|Subject),sleepstudy)
```

Many more options are available; please see the documentation (`?buildmer`, after installing the package) for details.

## Installation

```
library(devtools)
install_github('cvoeten/buildmer')
```

## Known issues

This package is work-in-progress. **You should always check the output `buildmer` provides during the fitting process, and verify that it is making the correct decisions in terms of which terms it chooses to add or remove from the model**. Bugs will be fixed as they are uncovered.

The correct evaluation of marginality restrictions is currently a slightly sore spot. It is surely broken in models with diagonal random-effect covariance structures (`y ~ ... + (1|g) + (0+x|g) + (0+y|g) + (0+x:y|g)` will happily consider `x|g` and `y|g` for exclusion, even though this is wrong). In other cases it also appears to sometimes fail, though only with the random effects. (If you notice a mistake, you can always correct it and instruct `buildmer` to leave the problematic parts alone, e.g.: `buildmer(...,reduce.random=FALSE)`.)

If you notice a bug, please file an issue for it and I may or may not look into it if I have time.

The design of buildmer is currently suboptimal. Most important functions are nested within the main `buildmer` function, which makes bookkeeping really easy, but is bad practice. Ideas on alternative approaches would be welcomed.
