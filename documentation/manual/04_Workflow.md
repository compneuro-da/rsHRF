<a name="top">
<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

:pencil2: Workflow
----

# Preprocessing
##Input?
The input is voxelwise/vertexwise BOLD signal, already preprocessed according to your favorite recipe. something like:
* nuisance variable regression 
* bandpass filter in the 0.01-0.08 Hz interval
* despike
(These denoising steps are also provided in the SPM plugin.)
The input can be images (3D or 4D), mesh (2D), or directly matrices of [observation x voxels(vertices)].
It is possible to use a temporal mask to exclude some time points (for example after scrubbing).
The demos allow you to run the analyses on several formats of input data.
temporal_mask: generated from scrubbing.
## Resting-state retrieval and deconvolution
--> standalone: script demo 
--> SPM plugin:
 SPM plugin
-------------

The script spm_rsHRF.m is the main one, and it calls rsHRF.m. These two files are specific to the SPM plugin. 

See [rsHRF_toolbox.pptx](https://github.com/guorongwu/rsHRF/raw/master/rsHRF_toolbox.pptx) for more details (Installation/Usage/Outputs).
![rsHRF GUI](https://github.com/guorongwu/rsHRF_data/raw/master/rsHRF_GUI.png)
    --> two videos
    --> can be visualized: how? one video; plus image batch
    --> batch demo (in .zip)
    
## Connectivity analysis
The connectivity analysis (functional connectivity: Pearson/Spearman correlation, Pairwise/Conditional/Partially Conditioned Granger causality) is only provided in the rsHRF SPM plugin. 
