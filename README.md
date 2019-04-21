
21 April 2019

This repo supports the paper "Biological Mechanisms for Learning: A Computational Model of Olfactory Learning in the Manduca sexta Moth, with Applications to Neural Nets", by CB Delahunt, JA Riffell, and JN Kutz, Frontiers in Computational Neuroscience, December 2018.

1. In vivo datasets
The folder 'inVivoDataPlottingCode' contains: 
(a) Data files of experiments (electrode readings from Antennal Lobes of live moths during learning); 
(b) Matlab code to plot the timecourses of neural firing rates (spontaneous FRs and odor responses).
Please consult the various 'inVivoDataPlottingCode/processAndPlot*.m' script headers for details about each in vivo experiment dataset. 


2. Computational model 
A full version of the moth simulation codebase is available now at github/charlesDelahunt/PuttingABugInML. 
That version, "MothNet", simulates a moth olfactory network as it learns to read MNIST digits. It generates a model of the Manduca sexta moth olfactory network, then runs a time-stepped evolution of SDEs to train the network to read downsampled MNIST digits.
The mechanics of neural architecture generation and SDE evolution in MothNet are almost identical to those of the version used in "Biological Mechanisms for Learning".

If you have questions or comments, please email Charles Delahunt at delahunt@uw.edu.

Thank you!
