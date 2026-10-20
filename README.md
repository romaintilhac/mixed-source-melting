# Mixed-Source Melting model (MSM3)

MATLAB implementation of the Mixed-Source Melting model (MSM3) to investigate the impact of recycled oceanic crust and its geochemical signatures in mid-ocean ridge basalts.

The package includes:

- The main script (*mainMixedSourceMelting.m*), which runs forward simulations or adaptive Monte Carlo simulations based on the MSM3 model.
- A post-processing script (*processMonteCarlo.m*) to process Monte Carlo simulations and plot their outputs.
- A user-defined model components file (*userModelComponents.m*) to specify the minerals, major elements, trace elements and isotopes included in the simulations.
- A user configuration file (*userFolderPaths.m*) to set the paths to the required model directories.
- A *setups* folder containing the configuration files used to initialize the model, including source compositions, model parameters and target geochemical constraints.
- A *resources* folder containing the databases and lookup tables required for the geochemical calculations.
- A *utils* folder containing the supporting functions required to run the model.

## Instructions & user changes

### 1. Prerequisites

- **[Alphamelts 1.9](https://magmasource.caltech.edu/alphamelts/)** is required to run pMELTS as part of this code and you need to make the following changes in _alphamelts_default_env.txt: *ALPHAMELTS_MINP* +0 and *ALPHAMELTS_SAVE_ALL* true.

- MATLAB with the **Statistics and Machine Learning Toolbox** is required to run the adaptive Monte Carlo simulations and plotting: _randsample_ performs weighted resampling, and _histogram_ is used to plot input distributions. The code has been tested in Matlab R2024b.

- Before running the model, make sure all paths in *userFolderPaths.m* are consistent with your installation, particularly the AlphaMELTS directory (*Alpha_install*) and executable (*Alpha_link*).

### 2. Setup

The main user-adjustable setups are defined as follows in subfolders of the *setups* folder:

- **Model options** (*options.txt*): control the modelling workflow, including the activation of Monte Carlo simulations and of each modelling stages, as well as the main model assumptions. When both stages are activated, options selected outputs from stage 1 can be passed to stage 2. The available options determine which compositions and parameters are transferred between stages and which are instead taken from the input files.
- **Model parameters** (*parameters.txt*): define the fixed parameters used in the simulations.
- **Source compositions** (*sources.txt*): specify the starting compositions and their references for each lithology, as well as modal, melting and partition coefficients. These value are read based on their name-reference in a user-adjustable lookup table (*input.xlsx*).
- **Monte Carlo ranges** (*montecarlo_ranges.txt*): define the ranges of input parameters and compositions sampled during Monte Carlo simulations, including which variables are adaptively resampled.
- **Monte Carlo targets** (*montecarlo_targets.txt*): define the target geochemical values and their associated uncertainties used to evaluate model fits.
- **Computational settings** (*settings.txt*): define the main parameters controlling the adaptive Monte Carlo algorithm (*e.g.,* number of iterations per pass, the maximum number of passes, the target goodness-of-fit and the stagnation criteria) as well as rounding and interpolation used in the model.

You may create new subfolders. The model setup that will be run is the one specified by *setup_name* in *userFolderPaths.m*. You may also want to define an output folder name (*folder_name*) accordingly.

***Note:*** Parameters and compositions may be overwritten by randomized inputs if also defined as Monte Carlo ranges.

The model components included in the calculations can also be modified with caution in *userModelComponents.m*. This file allows the user to specify the minerals, lithophile major elements and trace elements, chalcophile major and trace elements, and isotopes to be considered. When selecting an isotope, ensure that its parent element is included in the corresponding element list. If adding a new isotope, *buildTableIsotope.m* must also be updated accordingly for proper element indexing and plot labelling.

### 3. Running simulations

The model can then be directly launched by running *mainMixedSourceMelting.m* from the repository's main directory.

- For adaptive Monte Carlo simulations, set *MonteCarlo = 1* in *options.txt*. The model evaluates the simulated outputs against the target values and uncertainties specified in *montecarlo_targets.txt*. The resulting goodness-of-fit (GOF) values are used to update the sampling distributions between successive passes for variables marked for adaptation in *montecarlo_ranges.txt*. The number of iterations, maximum number of passes, target GOF and stagnation criteria are controlled through *settings.txt*.

- For a single forward simulation with fixed inputs, set *MonteCarlo = 0*. In this mode, the model runs once using the specified source compositions and model parameters.

The *resources* folder contains the geochemical databases and lookup tables required by the model. The *utils* folder contains the supporting functions used to initialize and run the simulations, perform the geochemical calculations and process the results.

***Note:*** the model relies on AlphaMELTS and Melt-PX calculations where enabled. When the corresponding options allow it, existing calculation outputs may be reused to reduce computation time. Changing the model configuration may require these calculations to be regenerated. Any change must be made consistent with the name-reference available in *input.xlsx* if relevant.

### 4. Processing and plotting

The post-processing of adaptive Monte Carlo simulations can be performed using *processMonteCarlo.m* after the simulations have finished. This script processes the saved outputs and plot the relevant figures.

**During adaptive Monte Carlo simulations** (MonteCarlo = 1), the model generates plots of the randomized input distributions and simulated outputs for each pass. These plots can be used to monitor the evolution of the sampling distributions and assess the quality of the model fits against the target geochemical constraints. The corresponding plotting scripts are *plotMCinputs.m* and *plotMCoutputs.m*. Compiled outputs and tabulated results are saved in the output directory specified by *folder_name* in *userFolderPaths.m*.

**During single forward simulations**(MonteCarlo = 0), the model automatically generates detailed plots of the model inputs and outputs. The corresponding plotting scripts are *plotStage2Inputs.m* and *plotStage2Outputs.m*.

## Reproducibility note



## Credits

By Romain Tilhac, CNRS researcher at Géosciences Montpellier, France. <br />
Updated on October 20th, 2026 <br />
Contact: romain.tilhac@umontpellier.fr

## References

König, S., Varas-Reus, M.I., Tilhac, R., Garrido, C. J. Selenium Isotopes Uncover the Mechanism of the Deep Plate Tectonic Cycle in the MORB Mantle. Science Advances. Manuscript no. aec9625 (preliminary acceptance).

