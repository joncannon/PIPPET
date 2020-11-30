# PIPPET/PATIPPET
 Phase (and tempo) Inference from Point Process Event Timing
 
 *See "PLOS pippet submission > PIPPET manuscript 11-6-20" for context!*
 
 Scripts to produce figures for this manuscript are in PIPPET.mlx and PATIPPET.mlx.
 
 Documentation of functions to follow...
 
 11/30/2020:
 Several important changes have been made to the code and manuscript since biorxiv posting:
 
 - PATIPPET math changed in order to scale the generative model's rate function with tempo. This is sensical because a faster tempo shouldn't mean we have overall weaker expectations at each expected time point, and it is practical because it allows us to estimate tempo based on event density alone.
 - An error was corrected. At events, the old math calculating the variance update was not correcting for the change of mean, and as a result it was overestimating variance after each jump. In the process of correcting this, it became useful to include old mean and new mean as separate inputs to the function calculating the updated variance.
 - Notation was totally revamped to make for a better presentation of the proof and less confusing conflation of PIPPET and PATIPPET. Thus, what was C in both models is now a matrix Sigma in PATIPPET and a scalar V in PIPPET. What was "mu" is now "x bar" in PATIPPET and "phi bar" in PIPPET. What was a bar is now a hat. Oy vei.

The images in the manuscript have not yet been updated to show the new notation or display the behavior of the new filters, but everything is basically the same.
