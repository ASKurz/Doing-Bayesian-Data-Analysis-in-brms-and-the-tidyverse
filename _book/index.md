--- 
title: "*Doing Bayesian Data Analysis* in brms and the tidyverse"
subtitle: "version 0.0.2"
author: ["A Solomon Kurz"]
date: "2019-10-27"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
geometry:
  margin = 0.5in
urlcolor: blue
highlight: tango
header-includes:
  \usepackage{underscore}
  \usepackage[T1]{fontenc}
link-citations: yes
github-repo: ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse
twitter-handle: SolomonKurz
description: "This project is an attempt to re-express the code in Kruschke's (2014) textbook. His models are re-fit in brms, plots are redone with ggplot2, and the general data wrangling code predominantly follows the tidyverse style."
---

# What and why {-}

Kruschke began his text with "This book explains how to actually do Bayesian data analysis, by real people (like you), for realistic data (like yours)." In the same way, this project is designed to help those real people do Bayesian data analysis. My contribution is converting Kruschke's JAGS code for use in Bürkner's [**brms** package](https://github.com/paul-buerkner/brms), which makes it easier to fit Bayesian regression models in **R** using Hamiltonian Monte Carlo (HMC). I also prefer plotting and data wrangling with the packages from the [**tidyverse**](http://style.tidyverse.org). So we'll be using those methods, too.

This project is not meant to stand alone. It's a supplement to the second edition of [Kruschke’s *Doing Bayesian Data Analysis*](https://sites.google.com/site/doingbayesiandataanalysis/). Please give the source material some love.

## Caution: Work in progress {-}

The first release of this project only contains Chapters 1 through 5. Welcome to version 0.0.2! The notable updates are:

* the addition of Chapters 6 through 10 and
* minor typo fixes to Chapter 5.

Version 0.0.2 is also noteworthy in that it's the first version containing an incomplete chapter. At present, I do not know how to perform the simulation for Figure 7.3, nor do I understand how to run Kruschke's effective sample size simulations from subsection 7.5.2. Fans of this project are welcome to share solutions to those sections in the [Issues section of the GitHub repository for this project](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues).

Most of the remaining chapters have completed drafts and just need another round of edits. I'll add them, soon. In addition to checking in here, you can follow my GitHub [progress log](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/1), in which I will point out other figures or sections I'm having trouble with.



