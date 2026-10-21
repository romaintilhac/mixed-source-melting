# Mixed-Source Melting model (MSM3) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23533320.svg)](https://doi.org/10.5281/zenodo.23533320)

MATLAB implementation of the Mixed-Source Melting model (MSM3) to investigate the geochemical signature of recycled oceanic crust in the mantle source of mid-ocean ridge basalts.

The package includes:

- The main script (*mainMixedSourceMelting.m*) to run single forward simulations or adaptive Monte Carlo simulations based on the MSM3 model.
- A post-processing script (*processMonteCarlo.m*) to process and plot Monte Carlo simulation outputs.
- A user-defined model components file (*userModelComponents.m*) to specify the minerals, major, trace elements and isotopes included in the simulations.
- A user configuration file (*userFolderPaths.m*) to set the paths to the required model directories and define the model setup and output folder.
- A *setups* folder containing the configuration files used to initialize the model, including source compositions, model parameters and target geochemical constraints.
- A *utils* folder containing all the supporting scripts and functions, including the two main model stages (*Stage1model.m* and *Stage2model.m*), with an explicit treatment of sulfur sytematics based on Sulfur Concentration at Sulfide Saturation (SCSS) calculations in *calcSCSS.m*.
- A *resources* folder containing the data required to run the model.
- An *outputs* folder containing an example of pre-run model results.

## Instructions & user changes

### Prerequisites

- MATLAB is required with the **[Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)** (_randsample_ and _histogram_ are need for the adaptive Monte Carlo simulations). The code has been tested in MATLAB R2024b.

- **[Alphamelts 1.9](https://magmasource.caltech.edu/alphamelts/)** is required to run pMELTS as part of this code. The following changes are needed in *alphamelts_default_env.txt*: *ALPHAMELTS_MINP* +0 and *ALPHAMELTS_SAVE_ALL* true.

### Configuration

Before running the model, make sure all paths in *userFolderPaths.m* are consistent with your installation, particularly the AlphaMELTS directory (_Alpha_install_) and executable (_Alpha_link_).

The code has been designed to be flexibly used with multiple geochemical proxies, including lithophile and chalcophile major, trace elements and isotopes. To that end, *userModelComponents.m* allows the user to specify the minerals, elements and isotopes to be considered. When selecting an isotope, ensure that its corresponding element is included in the relevant element list.

***Note:*** If adding a new isotope, *buildTableIsotope.m* must be updated accordingly for element indexing and labelling. Similarly, if new a mineral is added, *buildTableMineral.m* must be updated.

### Model setup

The main user-adjustable setups are defined as follows in subfolders of the *setups* folder:

- **Model options** (*options.txt*): control the modelling workflow, including the activation of Monte Carlo simulations and of each modelling stages, as well as the main model assumptions. When both stages are activated, specific options determine which compositions and parameters are passed from stage 1 to stage 2 or taken from the input files.
- **Model parameters** (*parameters.txt*): define the fixed parameters used in the simulations.
- **Source compositions** (*sources.txt*): specify the starting compositions and their references for each lithology, as well as modal, melting and partition coefficients. These value are read based on their names-references in a user-adjustable lookup table (*input.xlsx*).
- **Monte Carlo ranges** (*montecarlo_ranges.txt*): define the ranges of input parameters and compositions sampled during Monte Carlo simulations, including which variables are adaptively resampled.
- **Monte Carlo targets** (*montecarlo_targets.txt*): define the target geochemical values and their associated uncertainties used to evaluate model fits.
- **Computational settings** (*settings.txt*): define the main parameters controlling the adaptive Monte Carlo algorithm (*e.g.,* number of iterations per pass, the maximum number of passes, the target goodness-of-fit and the stagnation criteria) as well as rounding and interpolation used in the model.

You may create new subfolders. The model will run the setup specified by _setup_name_ in *userFolderPaths.m*. You may also want to define the output folder name (_folder_name_) accordingly. Any change must be made consistent with the names-references available in *input.xlsx* if relevant.

***Note:*** Parameters and compositions will be overwritten by randomized inputs if also defined as Monte Carlo ranges.

### Running simulations

The model can be directly launched by running *mainMixedSourceMelting.m* from the repository's main directory.

- **For adaptive Monte Carlo simulations**, set _MonteCarlo_ = 1 in *options.txt*. The model generates randomized input distributions as defined in *montecarlo_ranges.txt* and run iteratively. Simulated outputs are then evaluated against the target values and uncertainties specified in *montecarlo_targets.txt*. The resulting goodness-of-fit (GOF) values are used to update the sampling distributions between successive passes for variables marked for adaptation in *montecarlo_ranges.txt*. The number of iterations, maximum number of passes, target GOF and stagnation criteria are controlled through *settings.txt*.

- **For a single forward simulation**, set _MonteCarlo_ = 0. In this mode, the model runs once with fixed source compositions (*sources.txt*) and model parameters (*parameters.txt*).

***Note:*** the model relies on AlphaMELTS and Melt-PX calculations where enabled. Existing calculation outputs are stored in a *preload* folder and reused to reduce computation time. Changing the model configuration may require these calculations to be regenerated.

### Processing and plotting

The post-processing of adaptive Monte Carlo simulations can be performed using *processMonteCarlo.m* after the simulations have finished. This script processes the saved outputs and plot the relevant figures.

**During adaptive Monte Carlo simulations** (_MonteCarlo_ = 1), the model generates plots of the randomized input distributions and simulated outputs for each pass. These plots can be used to monitor the evolution of the sampling distributions and assess the quality of the model fits against the target geochemical constraints. The corresponding plotting scripts are *plotMCinputs.m* and *plotMCoutputs.m*. Compiled outputs and tabulated results are saved in the output directory specified by _folder_name_ in *userFolderPaths.m*.

**During single forward simulations** (_MonteCarlo_ = 0), the model automatically generates detailed plots of the model inputs and outputs. The corresponding plotting scripts are *plotStage2Inputs.m* and *plotStage2Outputs.m*.

## Reproducibility note

The model was applied to investigate the geochemical signatures of oceanic recycling in mid-ocean ridge basalts (MORB), based on the selenium (Se)-isotope compositions of pyroxenites, as presented in **König *et al.* 2026**.

The provided version of *userModelComponents.m* and setups corresponds to this application to &delta;<sup>82</sup>Se:

- *default_MC.txt* corresponds to the main Monte Carlo simulation with 20,000 iterations, which is aimed as a joint inversion of the Se-isotope compositions and Se concentrations of the studied pyroxenites and global MORB database (*data.xlsx*).
- *default_forward.txt* corresponds to a single forward model used to generate detailed model inputs and outputs.

To reproduce **König *et al.* 2026** Figures 2 and 3 and all the tabulated outputs:

1. Set _setup_name_ to *default_MC.txt* in *userFolderPaths.m*.
2. Run *mainMixedSourceMelting.m*.
3. Run *processMonteCarlo.m* to generate the figures and tabulated outputs.

Alternatively, because no pre-loaded calculation outputs are provided due to storage limitations, running the full Monte Carlo simulation from the scratch will take a fairly long time. The figures can therefore be also obtained from pre-run results provided in the *output* folder:

1. Set _setup_name_ to *default_MC_prerun* in *userFolderPaths.m*.
2. Run *processMonteCarlo.m* to generate the figures and tabulated outputs.

To reproduce Electronic Appendix Figures EA1-5:

1. Set _setup_name_ to *default_forward.txt* in *userFolderPaths.m*
2. Run *mainMixedSourceMelting.m*.
3. The figures will be automatically generated.

## Credits

By [**Romain Tilhac**](https://romaintilhac.github.io), CNRS researcher at Géosciences Montpellier, France. <br />
Updated on October 21th, 2026 <br />
Contact: [romain.tilhac@umontpellier.fr](mailto:romain.tilhac@umontpellier.fr)

## References

König, S., Varas-Reus, M.I., Tilhac, R., Garrido, C. J. Selenium Isotopes Uncover the Mechanism of the Deep Plate Tectonic Cycle in the MORB Mantle. *Science Advances*, aec9625 (preliminary acceptance).

