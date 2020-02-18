<a name="top">
<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

:pencil2: Workflow
----
__IMPORTANT__: Please use Google Chrome to browse the _Workflow_ Page! If you don't, you won't be able to expand the collapsibles. [<abbr title="Work In Progress"><i>WIP</i></abbr>]

# Preprocessing 
## Input:

<details><summary>What is the prefered data format?</summary> <!-- FAQ -->
<br> <!-- insert image: batch: scans -->
The <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox allows you to run the analyses on several formats of input data: <ul>
 <li>4D NIfTI</li>
 <li>3D NIfTI</li>
 <li>extracted signals – [observation x voxels/vertices] (.mat)</li>
 <li>2D surface-based (.gii)</li></ul>
 
<!-- The input can be images (3D or 4D), mesh (2D), or directly matrices of [observation x voxels/vertices]. The demos allow you to run the analyses on several formats of input data. -- as shown in the Flowchart [insert]-->

</details>
 
 Outlier removal is only legit when conducting a whole-brain analysis. Whole-brain vs. ROI? 


<!-- collapsibles -->
### File format

### z-scoring <!-- to check; already included [?] -->
### Denoising
The input is voxelwise/vertexwise BOLD signal, already preprocessed according to your favorite recipe. something like: <!-- cf. e-mail OHBM - what to add? -->
* nuisance variable regression 
* bandpass filter in the 0.01-0.08 Hz interval
* despike
(These denoising steps are also provided in the SPM plugin.)
It is possible to use a temporal mask to exclude some time points (for example after scrubbing).
temporal_mask: generated from scrubbing.

## rsHRF retrieval and deconvolution
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
