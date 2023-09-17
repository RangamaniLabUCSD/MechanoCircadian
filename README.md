# MechanoCircadian repository
MATLAB-based implementation of the coupling between the cell Circadian clock and mechanotransduction, associated with the manuscript by Emmet Francis and Padmini Rangamani, "Computational modeling predicts mechanotransduction-mediated changes to Circadian oscillations in mammalian cells".

To run the code in this repository, you need to first install the following:
* DDE-BIFTOOL - https://sourceforge.net/projects/ddebiftool/ (download code and add to path in MATLAB)
* UQLab - https://www.uqlab.com/download (download code and add to path in MATLAB)
* Violinplot-MATLAB - https://github.com/bastibe/Violinplot-Matlab (download code and add to path)
* linspecer - https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap (download code and add to path)
* Required MATLAB toolboxes: Communications Toolbox, Curve Fitting Toolbox, Global Optimization Toolbox, Optimization Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox.

Scripts used to generate the figures in the paper are organized into 5 main folders, "CircadianAnalysisAndBifurcation" (Fig 3, Fig S4), "MechanoAnalysis" (Fig 1B, Fig S3), "SensitivityAndFitting" (Fig 2, Fig S1, Fig S2), and "MechanoCircadian" (Fig 1B-C, Fig 4, Fig 5, Fig S5). We break down the associated files below.

## CircadianAnalysisAndBifurcation
* CircadianAnalysis: Includes the 2 DDE system to model Circadian oscillations, used to plot phase diagrams of oscillation period in the YAP/TAZ-MRTF phase plane for different parameter combinations.
* ddeBifCircadian: File using DDE-BIFTOOL for numerical bifurcation analysis of the DDE system, with the goal of plotting out the Hopf bifurcation in the YAP/TAZ-MRTF phase plane. Developed using code from the demos for DDE-BIFTOOL v3 (https://ddebiftool.sourceforge.net/demos/neuron/).

## MechanoAnalysis
* MechanoOnlyModel: function defining the YAP/TAZ and MRTF mechanotransduction model without any inclusion of the Circadian parameters.
* MechanoOnly_main: Script used to generate Fig 1B and Fig S3 and plot different test cases for YAP/TAZ and MRTF nuclear localization following different treatments.
* MechanoSS: function that returns steady state values for all state variables associated with the mechanotransduction model. Note that this is called in the MechanoCircadianModel function to compute steady state nuclear concentrations of YAP/TAZ and MRTF. 

## SensitivityAndFitting
* MechanoCircadian_sensitivity: driver script for Sobol sensitivity analysis using UQLab. (Fig 2A-B)
* uq_MechanoCircadian: "Forward model" for use in UQLab applications. Takes in a set of parameters and outputs quantities of interest (in this case, oscillation period, amplitude, and decay rate)
* MechanoCircadian_fit: script used for Bayesian parameter estimation and then plotting resulting distributions in oscillation period and amplitude for different treatment conditions (Fig 2, Fig S1, Fig S2).

## MechanoCircadian
* MechanoCircadianModel: Function implementing DDE system for modeling mechanotransduction-Circadian coupling. Calls MechanoSS to find SS nuclear concentrations of YAP/TAZ and MRTF. Includes a stochastic version of the DDEs that is not implemented in current paper.
* MechanoCircadian_main: Driver script used to test various cases of the MechanoCircadian model in the paper (Fig 1B-C, Fig 4, Fig 5, Fig S5). Includes basic testing of the mechano-Circadian model and implementation of population-level testing.
* conditionToOutputs: function taking in conditions for a MechanoCircadian test and outputting the oscillation period, amplitude, and decay rate, along with the oscillation dynamics. (calls MechanoCircadianModel)

## Utility functions
In addition to these main folders, we have several utility functions sorted into subfolders:
* external: function for computing integrated autocorrelation time (downloaded from https://www.physik.hu-berlin.de/de/com/UWerr_fft.m)
* osc_analysis: contains the function, circOscAnalysis, which takes in a time vector and state variable vector and outputs the oscillation period, amplitude, and decay rate (rate of change in amplitude), along with a vector specifying the indices associated with peaks.
* plotting: contains two functions, prctilePlot and prettyGraph. prctilePlot plots shaded regions associated with interquartile range and range from 2.5 quantile to 97.5 quantile. prettyGraph contains specifications for generating publication-ready plots.

## Contributing guidelines
Please feel free to submit any issues to this Github repository or suggest modifications by submitting a pull request.