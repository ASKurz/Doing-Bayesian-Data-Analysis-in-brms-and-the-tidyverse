--- 
title: "*Doing Bayesian Data Analysis* in brms and the tidyverse"
subtitle: "version 0.0.6"
author: ["A Solomon Kurz"]
date: "2020-01-25"
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
description: "This project is an attempt to re-express the code in Kruschke's (2015) textbook. His models are re-fit in brms, plots are redone with ggplot2, and the general data wrangling code predominantly follows the tidyverse style."
---

# What and why {-}

Kruschke began his text with "This book explains how to actually do Bayesian data analysis, by real people (like you), for realistic data (like yours)." In the same way, this project is designed to help those real people do Bayesian data analysis. My contribution is converting Kruschke's JAGS code for use in Bürkner's [**brms** package](https://github.com/paul-buerkner/brms), which makes it easier to fit Bayesian regression models in **R** using Hamiltonian Monte Carlo (HMC). I also prefer plotting and data wrangling with the packages from the [**tidyverse**](http://style.tidyverse.org). So we'll be using those methods, too.

This project is not meant to stand alone. It's a supplement to the second edition of [Kruschke's *Doing Bayesian Data Analysis*](https://sites.google.com/site/doingbayesiandataanalysis/). Please give the source material some love.

## Caution: Work in progress {-}

The first release of this project only contained Chapters 1 through 5. Version 0.0.2 added Chapters 6 through 10. Version 0.0.3 added Chapters 11, 12, and 15 through 18, with placeholders for Chapters 13 and 14. Version 0.0.4 added Chapters 19 through 21. Version 0.0.5 added Chapter 22 and a bit of prose to Chapters 13 and 14. Welcome to version 0.0.6! The notable updates are:

* the addition of Chapters 23 and 24 and 
* minor typo fixes throughout.

Several sections still remain in Chapters 7 and 11 that I have yet to master. They and other difficulties are detailed in the GitHub [progress log](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues/1). If you would like to share your solutions to those sections or any other section in the project, please open an [issue in this project's GitHub repository](https://github.com/ASKurz/Doing-Bayesian-Data-Analysis-in-brms-and-the-tidyverse/issues). It would be a tremendous help.

Content will continue trickling in as I get to it.

## Thank-you's are in order {-}

Before we enter the primary text, I'd like to thank the following for their helpful contributions:

* Paul-Christian Bürkner ([\@paul-buerkner](https://github.com/paul-buerkner)), 
* Andrew Gelman ([\@andrewgelman](https://github.com/andrewgelman)), 
* Matthew Kay ([\@mjskay](https://github.com/mjskay)), 
* TJ Mahr ([\@tjmahr](https://github.com/tjmahr)), 
* Lukas Neugebauer  ([\@LukasNeugebauer](https://github.com/LukasNeugebauer)), 
* Demetri Pananos ([\@Dpananos](https://github.com/dpananos)), 
* Aki Vehtari ([\@avehtari](https://github.com/avehtari)), 
* Matti Vuorre ([\@mvuorre](https://github.com/mvuorre)), and
* Brenton M. Wiernik ([\@bwiernik](https://github.com/bwiernik)).



