# MT LIP GLM

This Repo provides MATLAB code for the basic GLM analyses in Yates et al. (2017)

## Getting started / installation

### Installing the code

Create a local copy of mtlipglm by cloning or downloading a zip

```shell
	git clone https://github.com/jcbyts/mtlipglm
```

You will also need the `classy` branch of neuroGLM code

```shell
	git clone https://github.com/jcbyts/neuroGLM
	cd neuroGLM
	git checkout classy
```

### Getting the data

We plan to release the data publically, but there's no telling when that will happen. Email me for a link to download sample data.
<!-- The code will not run without local copies of the data. 

To quickly reproduce figures from the manuscript requires both the data and the fitted GLMs. Because all the data and fits takes up a fair amount of disk space, there are a number of options.

1. [neurons, stimulus (with eye position)](https://www.dropbox.com/s/ix7vtrid2rlfuyx/mtlipglm_data_full.zip?dl=0)
2. [neurons, stimulus (no eye position data)](https://www.dropbox.com/s/0myzjnh5xy014pi/mtlipglm_data_small.zip?dl=0)
3. [LIP fits](https://www.dropbox.com/s/xc5n2wh02wsjhzg/mtlipglm_lip_fits.zip?dl=0)

You can download the data through #1 or #2 and then refit yourself, or additionally download #3 as well for a quick reproduction of figure 5.

If you downloaded the fits, move `main_fits` and `lip_trunc_fits` to the directory where neurons and stim are. The `mtlipglm` code will assume that they are all in the same base directory. -->

### Before running
Note: this code has only been tested on MATLAB 2015b and newer.


### setup paths
Open Matlab. Change directories to mtlipglm. You need to set up the specific paths for your machine.

```matlab
edit addMTLIPpaths
```

Find the lines

```matlab
% this is where your data live
dataPath = '';
```

and

```matlab
% neuroGLM - must be downloaded from https://github.com/jcbyts/neuroGLM 
neuroGLMpath = '';
```

and enter the appropriate path locations for the base data directory and the neuroGLM repository.

run `addMTLIPpaths` and if setup correctly, then 

```matlab
ls(getpref('mtlipglm', 'dataPath'))
```
should show you the `neurons`, `stim`, `main_fits`, and `lip_trunc_fits` directories.


### Getting started
To begin exploring the data, the script `cmdExamplePTA.m` walks through the pulse-triggered average analyses from figure 2. The script works as a walkthrough and should be self-explanatory.
```matlab
edit cmdExamplePTA.m
```

### Example GLM
To explore the data / try alternate GLM parameterizations,
```matlab
edit cmdExample.m
```

### A note on fitting procedures and hyperparameters
There are many parameters in the GLM and fitting usually requires some regularization. The GLMs here are regularized in two ways: The first is that the weights are parameterized on user-specified basis functions. The basis functions that are chosen specify the resolution of the fit and impose some smoothness. The second is that the weights on the basis functions are fit with ridge regression. This imposes shrinkage and biases the weights to be small. These hyperparameters can change the exact values of the weights and the resulting model fits by a non-trivial amount, however the fits are usually very robust to many settings of the basis functions and hyperparameters. To get a feel for how the basis and the ridge paramter changes the fitted models, explore the code in cmdExample. In the fits that are provided in the download link above, the basis functions were set to be identical for all units from each area. The ridge parameter was learned by searching over a grid on witheld data from the training set. The fitting code provided here fits the GLM with a fixed ridge parameter. The resulting fits will be subtly different than the ones provided in the link above.

## Basic Overview

There are a few classes that govern the analyses in mtlipglm. 

Try running
```matlab
dataPath = getpref('mtlipglm', 'dataPath');
neurons = getNeurons('p20140304', dataPath)
```

You'll notice that `neurons` is an array of `neuro` objects. This is the most memory efficient way to access the neuron data. Each neuron file is an [hdf5](https://www.google.com/search?q=hdf5) file that includes many attributes about the unit under study (e.g., receptive field maps, waveform properties). Rather than load all of these files and atributes into memory, the `neuro` object can access each attribute when needed. To get the full struct for a neuron, just run

```matlab
n = neurons(1).getStruct;
```
or load the file directly. However, as mentioned, the code is going to use the `neuro` class, both because of its lazy access properties, but also because it has a number of useful methods for getting more info about the neuron. The `neuro` class has a property `stim` that similarly accesses the stimulus properties when needed.

The rest of the code revolves around two classes: `mtlipglm` and `glmspike` 

`glmspike` is an inherits `neuroGLM` and can be used the same way. It contains additional methods for setting up particular parameterizations of the basis functions, compiling the design matrix, and performing the fitting.

`mtlipglm` is the workhorse. It creates a `trial` structure that is required by `neuroGLM`. It sets up the specific `glmspike` objects for each GLM and fits them. It also has model comparison code to analyze the data.

## Figure code

This repository can reproduce figure 5. It has all the fitting code for fixed hyperparameters. It also has an example script for exploring the data and creating new GLMs.

### Figure 5 b-f
If you have downloaded the fits, either run `figure05` or open it and run cell by cell. This will reproduce figure 5, panels b-f from the manuscript.

If you have not downloaded the fits, you will need to run `cmdFitLIP`. This will probably take some time (30min+) to refit all the required models with a fixed hyperparameter. Once this is completed running, you can run `figure05` to plot the results.

### Figure 5 g,h
To reproduce the effects of elongating the acausal portion of the pre-saccadic covariates, run `cmdFitLIPtrunc`. If you downloaded the fits, this will re-analyze them and produce the required data for plotting results. Otherwise, it will refit each model. You can adjust the number of truncation steps to make this faster.









