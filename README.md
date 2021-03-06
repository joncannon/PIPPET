# PIPPET/PATIPPET
 Phase (and tempo) Inference from Point Process Event Timing
Jonathan Cannon, MIT, 2020
 
 *See "Manuscript drafts > Revision_1_PIPPET.pdf" for essential context!*
 
 Brief guide to code:
-  Scripts to produce figures for the manuscript are in PIPPET.mlx and PATIPPET.mlx.
-  They call "run_PIPPET" and "run_PATIPPET" functions, respectively, to run simulations.
-  These functions return a history of posterior mean and variance/covariance over the simulation.
-  Default parameter sets for the simulations are generated by "PIPPET_params" and "PATIPPET_params", which take optional arguments to set different values of parameters.
-  These functions call subordinate functions "PIPPET_stream_params" and "PATIPPET_stream_params" once for each separate stream of point process events (only one stream for all simulations in the paper, but this code should generalize nicely to multiple streams).
-  To run a simulation with multiple streams, parameters can be entered for each stream as a cell array, or entered as scalars to apply to all streams.
 
 
 11/30/2020:
 Several important changes have been made to the code and manuscript since biorxiv posting:
 
 - PATIPPET math changed in order to scale the generative model's rate function with tempo. This is sensical because a faster tempo shouldn't mean we have overall weaker expectations at each expected time point, and it is practical because it allows us to estimate tempo based on event density alone.
 - An error was corrected. At events, the old math calculating the variance update was not correcting for the change of mean, and as a result it was overestimating variance after each jump. In the process of correcting this, it became useful to include old mean and new mean as separate inputs to the function calculating the updated variance.
 - Notation was totally revamped to make for a better presentation of the proof and less confusing conflation of PIPPET and PATIPPET. Thus, what was C in both models is now a matrix Sigma in PATIPPET and a scalar V in PIPPET. What was "mu" is now "x bar" in PATIPPET and "phi bar" in PIPPET. What was a bar is now a hat. Oy vei.

The images in the manuscript have not yet been updated to show the new notation or display the behavior of the new filters, but everything is basically the same.
