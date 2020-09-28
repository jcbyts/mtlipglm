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
An example session can be found [here](https://drive.google.com/file/d/1RhoC3VIxQ_dBem1oPngY-m8f-_vfi8cc/view)

### Before running
Note: this code has only been tested on MATLAB versions between 2015b and 2018a.


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
% neuroGLM - must be downloaded from https://github.com/jcbyts/neuroGLM -- you want the classy branch for this to work properly 
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
neurons = getNeurons('p20140307', dataPath)
```

You'll notice that `neurons` is an array of `neuro` objects. This is the most memory efficient way to access the neuron data. Each neuron file is an [hdf5](https://www.google.com/search?q=hdf5) file that includes many attributes about the unit under study (e.g., receptive field maps, waveform properties). Rather than load all of these files and atributes into memory, the `neuro` object can access each attribute when needed. To get the full struct for a neuron, just run

```matlab
n = neurons(1).getStruct;
```
or load the file directly. However, as mentioned, the code is going to use the `neuro` class, both because of its lazy access properties, but also because it has a number of useful methods for getting more info about the neuron. The `neuro` class has a property `stim` that similarly accesses the stimulus properties when needed.

The rest of the code revolves around two classes: `mtlipglm` and `glmspike` 

`glmspike` is an inherits `neuroGLM` and can be used the same way. It contains additional methods for setting up particular parameterizations of the basis functions, compiling the design matrix, and performing the fitting.

`mtlipglm` is the workhorse. It creates a `trial` structure that is required by `neuroGLM`. It sets up the specific `glmspike` objects for each GLM and fits them. It also has model comparison code to analyze the data.








