--- 
title: "*Doing Bayesian Data Analysis* in brms and the tidyverse"
subtitle: "version 1.1.0"
author: "A Solomon Kurz"
date: "2023-01-25"
site: bookdown::bookdown_site
output: 
  bookdown::gitbook
documentclass: book
bibliography: bib_zotero.bib
biblio-style: apalike
csl: apa.csl
link-citations: yes
geometry:
  margin = 0.5in
urlcolor: blue
highlight: tango
header-includes:
  \usepackage{underscore}
  \usepackage[T1]{fontenc}
github-repo: ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse
twitter-handle: SolomonKurz
description: "This project is an attempt to re-express the code in Kruschke's (2015) textbook. His models are re-fit in brms, plots are redone with ggplot2, and the general data wrangling code predominantly follows the tidyverse style."
---

# What and why {-}

Kruschke began his text with "This book explains how to actually do Bayesian data analysis, by real people (like you), for realistic data (like yours)." In the same way, this project is designed to help those real people do Bayesian data analysis. My contribution is converting Kruschke's JAGS and Stan code for use in Bürkner's [**brms** package](https://github.com/paul-buerkner/brms) [@R-brms; @burknerBrmsPackageBayesian2017; @burknerAdvancedBayesianMultilevel2018], which makes it easier to fit Bayesian regression models in **R** [@R-base] using Hamiltonian Monte Carlo. I also prefer plotting and data wrangling with the packages from the [**tidyverse**](http://style.tidyverse.org) [@R-tidyverse; @wickhamWelcomeTidyverse2019]. So we'll be using those methods, too.

This ebook is not meant to stand alone. It's a supplement to the second edition of Kruschke's [-@kruschkeDoingBayesianData2015] [*Doing Bayesian data analysis: A tutorial with R, JAGS, and Stan*](https://sites.google.com/site/doingbayesiandataanalysis/). Please give the source material some love.

## **R** setup {-}

To get the full benefit from this ebook, you'll need some software. Happily, everything will be free (provided you have access to a decent personal computer and an good internet connection).

First, you'll need to install **R**, which you can learn about at [https://cran.r-project.org/](https://cran.r-project.org/).

Though not necessary, your **R** experience might be more enjoyable if done through the free RStudio interface, which you can learn about at [https://rstudio.com/products/rstudio/](https://rstudio.com/products/rstudio/).

Once you have installed **R**, execute the following to install the bulk of the add-on packages. This will probably take a few minutes to finish. Go make yourself a coffee.


```r
packages <- c("bayesplot", "brms", "coda", "cowplot", "cubelyr", "devtools", "fishualize", "GGally", "ggdist", "ggExtra", "ggforce", "ggmcmc", "ggridges", "ggthemes", "janitor", "lisa", "loo", "palettetown", "patchwork", "psych", "remotes", "rstan", "santoku", "scico", "tidybayes", "tidyverse")

install.packages(packages, dependencies = T)
```

A few of the other packages are not officially available via the Comprehensive R Archive Network (CRAN; https://cran.r-project.org/). You can download them directly from GitHub by executing the following.


```r
remotes::install_github("clauswilke/colorblindr")
devtools::install_github("dill/beyonce")
devtools::install_github("ropenscilabs/ochRe")
```

It's possible you'll have problems installing some of these packages. Here are some likely suspects and where you can find help:

* for difficulties installing **brms**, go to [https://github.com/paul-buerkner/brms#how-do-i-install-brms](https://github.com/paul-buerkner/brms#how-do-i-install-brms) or search around in the [**brms** section of the Stan forums](https://discourse.mc-stan.org/c/interfaces/brms/36); and
* for difficulties installing **rstan**, go to [https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

## We have updates {-}

For a brief rundown of the version history, we have:

### Version 0.1.0. {-}

I released the 0.1.0 version of this project in February 17, 2020. It was the first [fairly] complete draft including material from all the chapters in Kruschke's text. The supermajority of Kruschke's JAGS and Stan models were fit **brms** 2.11.5. The results were saved in the [`fits` folder on GitHub](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/tree/master/fits) and most of the results are quite comparable to those in the original text. We also reproduced most of the data-related figures and tables and little subpoints and examples sprinkled throughout Kruschke's prose.

### Version 0.2.0. {-}

The 0.2.0 update came in May 19, 2020. Noteworthy changes included:

* reproducing the simulation necessary for Figure 7.3 (see [GitHub issue #14](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/14)) with help from Cardy Moten III ([\@cmoten](https://github.com/cmoten));
* with guidance from Bjørn Peare Bartholdy ([\@bbartholdy](https://github.com/bbartholdy)), Mladen Jovanović ([\@mladenjovanovic](https://github.com/mladenjovanovic)), Cory Whitney ([\@CWWhitney](https://github.com/CWWhitney)), and Brenton M. Wiernik ([\@bwiernik](https://github.com/bwiernik)), we improved in-text citations and reference sections using [BibTex](http://www.bibtex.org/) [@BibTeX2020], [Better BibTeX](https://github.com/retorquere/zotero-better-bibtex) [@BetterBibTeXZotero2020], and [zotero](https://www.zotero.org/) [@ZoteroYourPersonal2020];
* the plot resolution increased with `fig.retina = 2.5`; and
* small code, hyperlink, and typo corrections.

### Version 0.3.0. {-}

The 0.3.0 update came in September 22, 2020. Noteworthy changes included:

* adding the [Kruschke-style model diagrams](https://solomonkurz.netlify.app/post/make-model-diagrams-kruschke-style/) throughout the text (e.g., [Figure 8.5][Example: Difference of biases]);
* adding chapter-specific plotting schemes with help from the [**cowplot** package](https://wilkelab.org/cowplot) [@R-cowplot], Wilke's [-@wilkeFundamentalsDataVisualization2019] [*Fundamentals of data visualization*](https://clauswilke.com/dataviz/), and many other great color-scheme packages; 
* an overhaul to the plotting workflow in [Section 6.4.1][Prior knowledge expressed as a beta distribution.]; and
* updating all model fits with **brms** version 2.13.5.

### Version 0.4.0. {-}

The 0.4.0 update came in May 6, 2021. Noteworthy changes included:

* using the Metropolis algorithm to fit the bivariate Bernoulli model for Figure 7.6 ([Section 7.4.3][The posterior via the Metropolis algorithm.]), thanks to help from [Omid Ghasemi](https://github.com/OmidGhasemi21);
* corrections to mistakes around the `lag()` and `lead()` functions in [Section 7.5.2][MCMC accuracy.];
* an added bonus section clarifying the pooled standard deviation for standardized mean differences ([Section 16.3.0.1][Bonus: Pooled standard deviation.]);
* refining the custom `stat_wilke()` plotting function in [Chapter 18][Metric Predicted Variable with Multiple Metric Predictors];
* an overhaul of the bonus section covering effect sizes ([Section 19.6][~~Exercises~~ Walk out an effect size]);
* refining/correcting the threshold workflow for univariable logistic regression models ([Chapter 21][Dichotomous Predicted Variable]);
* fixing the divergent transitions issue for the robust logistic regression model by adding boundaries on the prior ([Section 21.3][Robust logistic regression]);
* corrections to a few incorrectly computed effect sizes in [Chapter 23][Ordinal Predicted Variable];
* the addition of a new bonus section ([Section 22.3.3.1.1][Bonus: Consider the interceps-only softmax model.]) highlighting the benefits of the intercepts-only softmax model;
* expansions to the material on censored data ([Section 25.4][Censored Data in ~~JAGS~~ brms]) and the addition of a brief introduction to truncated data ([Section 25.4.4][Bonus: Truncation.]); and
* updating all HMC fits to the current version of **brms** (2.15.0).

### Version 1.0.0. {-}

The big 1.0.0 update came in May 4, 2022. This was the first full draft in the sense that it contained **brms** versions of all of Kruschke's JAGS and Stan models, excluding examples not currently possible with the **brms** paradigm (e.g., [Section 10.3.2][Hierarchical MCMC computation ~~of relative model probability~~ is not available in brms: We'll cover information criteria instead.]). Noteworthy changes included:

* two new solutions for the conditional logistic models of [Chapter 22][Nominal Predicted Variable] thanks to the generous efforts by [Henrik Singmann](https://github.com/singmann) and [Mattan Ben-Shachar](https://github.com/mattansb/);
* replacing the depreciated `posterior_samples()` function with the new `posterior::as_draws_df()`-based workflow;
* adding a new solution for the multivariate Bernoulli model with different trial numbers via the `resp_subset()` function in [Section 7.4.4.1][Uneven multivariate Benoulli via the `resp_subset()` approach.];
* improving the efficiency of the intercept-only Bernoulli models with the new `lb` and `ub` arguments for priors of `class = Intercept` (e.g., `fit7.1a`, `fit7.1b`, and `fit7.2` in [Chapter 7][Markov Chain Monte Carlo]);
* updating all model fits with **brms** version 2.17.0; and
* various minor code, hyperlink, and typo corrections.

### Version 1.1.0. {-}

Welcome to version 1.1.0! Noteworthy changes include:

* replacing my incorrect use of `tidyr::expand()` with a more appropriate `tidyr::expand_grid()` workflow, thanks to insights from [Desislava Petkova](https://github.com/dipetkov);
* adopting the new `linewidth` argument for several **ggplot2** geoms (see [here](https://www.tidyverse.org/blog/2022/08/ggplot2-3-4-0-size-to-linewidth/)); and
* minor prose and code edits throughout.

### Some minor issues remain. {-}

There are some minor improvements I'd like to add in future versions. Most importantly, a few simulations, figures, and models are beyond my current skill set. I've opened separate GitHub issues for the most important ones and they are as follows:

* the effective-sample-size simulations in Section 7.5.2 and the corresponding plots in Figures 7.13 and 7.14 ([issue #15](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/15)),
* several of the simulations in Sections 11.1.4, 11.3.1, and 11.3.2 and their corresponding figures (issues [#16](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/16), [#17](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/17), [#18](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/18), and [#19](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/19)),
* the stopping-rule simulations in Section 13.3.2 and their corresponding figures ([issue #20](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/20)),
* the data necessary to properly reproduce the HMC proposal schematic presented in Section 14.1 and Figures 14.1 through 14.3 ([issue #21](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/21)), and

If you know how to conquer any of these unresolved challenges, I'd love to hear all about it. In addition, please feel free to open a new [GitHub issue](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues) if you find any flaws in the other sections of the ebook.

## Thank-you's are in order {-}

Before we enter the primary text, I'd like to thank the following for their helpful contributions:

* Bjørn Peare Bartholdy ([\@bbartholdy](https://github.com/bbartholdy)),
* David Baumeister ([\@xdavebx](https://github.com/xdavebx)),
* Mattan Ben-Shachar ([\@mattansb](https://github.com/mattansb/)),
* Paul-Christian Bürkner ([\@paul-buerkner](https://github.com/paul-buerkner)),
* Christopher Flach ([\@flachboard](https://github.com/flachboard)),
* Andrew Gelman ([\@andrewgelman](https://github.com/andrewgelman)),
* Omid Ghasemi ([\@OmidGhasemi21](https://github.com/OmidGhasemi21)),
* Mladen Jovanović ([\@mladenjovanovic](https://github.com/mladenjovanovic)),
* Matthew Kay ([\@mjskay](https://github.com/mjskay)),
* TJ Mahr ([\@tjmahr](https://github.com/tjmahr)),
* Cardy Moten III ([\@cmoten](https://github.com/cmoten)),
* Lukas Neugebauer ([\@LukasNeugebauer](https://github.com/LukasNeugebauer)),
* Demetri Pananos ([\@Dpananos](https://github.com/dpananos)),
* Desislava Petkova ([\@dipetkov](https://github.com/dipetkov)),
* Peter Ralph ([\@petrelharp](https://github.com/petrelharp)),
* Henrik Singmann ([\@singmann](https://github.com/singmann)),
* Aki Vehtari ([\@avehtari](https://github.com/avehtari)),
* Matti Vuorre ([\@mvuorre](https://github.com/mvuorre)),
* Cory Whitney ([\@CWWhitney](https://github.com/CWWhitney)), and
* Brenton M. Wiernik ([\@bwiernik](https://github.com/bwiernik)).

## License and citation {-}

This book is licensed under the Creative Commons Zero v1.0 Universal license. You can learn the details, [here](https://github.com/ASKurz/Applied-Longitudinal-Data-Analysis-with-brms-and-the-tidyverse/blob/master/LICENSE). In short, you can use my work. Just please give me the appropriate credit the same way you would for any other scholarly resource. Here's the citation information:


```r
@book{kurzDoingBayesianDataAnalysis2023,
  title = {Doing {{Bayesian}} data analysis in brms and the tidyverse},
  author = {Kurz, A. Solomon},
  year = {2023},
  month = {1},
  edition = {Version 1.1.0},
  url = {https://bookdown.org/content/3686/}
}
```

