---
title: "Area of applicability of spatial prediction models"
author: "Hanna Meyer"
date: "2020-06-19"
output:
  rmarkdown::html_document:
    toc: true
    theme: united
vignette: >
  %\VignetteIndexEntry{Area of applicability of spatial prediction models}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



---


# Introduction
In spatial predictive mapping, models are often applied to make predictions far beyond sampling locations (i.e. field observarions used to map a variable even on a global scale), where new locations might considerably differ in their environmental properties. However, areas in the predictor space without support of training data are problematic. The model has no knowledge about these environments and predictions for such areas have to be considered highly uncertain. 

Here we implement the methodology described in [Meyer\&Pebesma (2020)](http://arxiv.org/abs/2005.07939) to estimate the "area of applicability" (AOA) of spatial prediction models. The AOA is defined as the area for which, in average, the cross-validation error of a trained model applies. To delineate the AOA, first an dissimilarity index (DI) is calculated that is based on distances to the training data in the multidimensional predictor variable space. To account for relevance of predictor variables responsible for prediction patterns we weight variables by the model-derived importance scores prior to distance calculation. The AOA is then derived by applying a threshold based on the DI observed in the training data.
 
This tutorial shows an example of how to estimate the area of applicability of spatial prediction models. 

For further information see: Meyer, H., Pebesma, E. (2020): Predicting into unknown space? Estimating the area of applicability of spatial prediction models. https://arxiv.org/abs/2005.07939

### Getting started














































