# Mixed-Source Melting model (MSM3)

MATLAB implementation of the Mixed-Source Melting model (MSM3) to investigate the impact of recycled oceanic crust and its geochemical signatures in mid-ocean ridge basalts.

The package includes:
- The main script (*mainMixedSourceMelting.m*), which runs forward simulations or adaptive Monte Carlo simulations based on the MSM3 model.
- A user-defined model components file (*userModelComponents.m*) to specify the minerals, major elements, trace elements and isotopes included in the simulations.
- A user configuration file (*userFolderPaths.m*) to set the paths to the required model directories.
- A *setups* folder containing the configuration files used to initialize the model, including source compositions, model parameters and target geochemical constraints.
- A *resources* folder containing the databases and lookup tables required for the geochemical calculations.
- A *utils* folder containing the supporting functions required to run the model.
- A post-processing script (*processMonteCarlo.m*) to process Monte Carlo simulations and plot their outputs.

The model is applied to investigate the origin and evolution of geochemical signatures in the mid-ocean ridge basalt (MORB) mantle, with particular emphasis on selenium isotopes and the mechanisms involved in the deep recycling of oceanic plates, as presented in Koenig et al., *Selenium Isotopes Uncover the Mechanism of the Deep Plate Tectonic Cycle in the MORB Mantle* (*Science Advances*, manuscript no. aec9625).

## Instructions & user changes

- [Alphamelts 1.9](https://magmasource.caltech.edu/alphamelts/) is required to run pMELTS as part of this code and you need to make the following changes in _alphamelts_default_env.txt:
    %   ALPHAMELTS_MINP +0
    %   ALPHAMELTS_SAVE_ALL true

- MATLAB with the **Statistics and Machine Learning Toolbox** is required to run the adaptive Monte Carlo simulations and plotting: _randsample_ performs weighted resampling, and _histogram_ is used to plot sampled input distributions. The code has been tested in Matlab R2024b.

- Before running the model, make sure all paths in *userFolderPaths.m* are consistent with your installation, particularly the AlphaMELTS directory (*Alpha_install*) and executable (*Alpha_link*).

The user can also select the model setup (*setup_name*) and output folder (*folder_name*).

The model can be directly launched by running *mainMixedSourceMelting.m* from the repository's main directory.

The main user-adjustable settings are defined in the *setups* folder, in the subfolder specified by *setup_name* (default: *default*). The following files control the model configuration:

- **Model options** (*options.txt*): control the modelling workflow, including the activation of Monte Carlo simulations, the two modelling stages, melting calculations, refertilization and sulfur-related calculations.
- **Model parameters** (*parameters.txt*): define the fixed parameters used in the simulations.
- **Source compositions** (*sources.txt*): specify the starting compositions and their references for peridotite, eclogite and pyroxenite.
- **Monte Carlo ranges** (*montecarlo_ranges.txt*): define the ranges of input parameters and compositions sampled during Monte Carlo simulations, including which variables are adaptively resampled.
- **Monte Carlo targets** (*montecarlo_targets.txt*): define the target geochemical values and their associated uncertainties used to evaluate model fits.

The model components included in the calculations can be modified in *userModelComponents.m*. This file allows the user to specify the minerals, lithophile major elements and trace elements, chalcophile major and trace elements, and isotopes to be considered. When selecting an isotope, ensure that its parent element is also included in the corresponding element list.

The main workflow is controlled through *options.txt*. The two modelling stages can be activated independently: stage 1 addresses source mixing and, depending on the selected options, melting; stage 2 models subsequent melting and melt–mantle interactions. When both stages are activated, selected outputs from stage 1 can be passed to stage 2. The available options determine which compositions and parameters are transferred between stages and which are instead taken from the input files.

For adaptive Monte Carlo simulations, set *MonteCarlo = 1* in *options.txt*. The model evaluates the simulated outputs against the target values and uncertainties specified in *montecarlo_targets.txt*. The resulting goodness-of-fit values are used to update the sampling distributions between successive passes for variables marked for adaptation in *montecarlo_ranges.txt*. The number of iterations, maximum number of passes, target goodness of fit and stagnation criteria are controlled through *settings.txt*.

For a single forward simulation with fixed inputs, set *MonteCarlo = 0*. In this mode, the model runs once using the specified source compositions and model parameters.

The *resources* folder contains the geochemical databases and lookup tables required by the model. The *utils* folder contains the supporting functions used to initialize and run the simulations, perform the geochemical calculations and process the results.

Note that the model relies on external AlphaMELTS and pMELTS calculations, as well as Melt-PX calculations where enabled. Existing calculation outputs may be reused when the corresponding options allow it, reducing computation time. Changing the model configuration may require these calculations to be regenerated.

## Reproducibility note

## Credits

By Romain Tilhac, CNRS researcher at Géosciences Montpellier, France. <br />
Updated on October 20th, 2026 <br />
Contact: romain.tilhac@umontpellier.fr

## References

König, S., Varas-Reus, M.I., Tilhac, R., Garrido, C. J. Selenium Isotopes Uncover the Mechanism of the Deep Plate Tectonic Cycle in the MORB Mantle. Science Advances. Manuscript no. aec9625 (preliminary acceptance).

