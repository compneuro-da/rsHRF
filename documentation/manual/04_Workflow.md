<a name="top">
<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

:pencil2: Workflow
----
__IMPORTANT__: Please use Google Chrome to browse the _Workflow_ Page! If you don't, you won't be able to expand the collapsibles. [<abbr title="Work In Progress"><i>WIP</i></abbr>]

# Preprocessing 
## Input:

<details><summary><b>What is the prefered input data format?</b></summary> <!-- FAQ -->
<br> <!-- insert image: batch: scans -->
<!-- <img align="right" src="https://github.com/compneuro-da/rsHRF/blob/update/img/input_01.png" alt="Input_Format" width="200"/> -->
<p align="justify">The <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox allows you to run the analyses on several formats of input data: <i>3D NIfTI</i>, <i>4D NIfTI</i>, <i>2D surface-based (.gii) files</i>, and <i>extracted signals (.mat) – [observation x voxels/vertices]</i>.</p>
<!-- The input can be images (3D or 4D), mesh (2D), or directly matrices of [observation x voxels/vertices]. The demos allow you to run the analyses on several formats of input data. As shown in the Flowchart [insert] -->
<!-- example for every kind of input; .mat is ok (1 van de drie?); not for other two 
examples are tested using MATLAB R2015b + which spm version
<!--
- nifti (3d & 4d)
- mat: stand (ok)
- 2g -->

</details>

<details><summary><b>Should the input data be standardized (i.e. z-scored) a priori?</b></summary><br> <!-- FAQ -->

<p align="justify">No, the standardization of the resting-state <abbr title="functional Magnetic Resonance Imaging">fMRI</abbr> <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> time series has already been included in the code for the <abbr title="hemodynamic response function">HRF</abbr> basis functions which you can find in the <code>rsHRF/code/</code> folder (i.e. <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code>, <code>wgr_rsHRF_FIR.m</code>, <code>rsHRF_estimation_FIR.m</code>, <code>rsHRF_estimation_temporal_basis.m</code></a>, and <code>rsHRF_estimation_impulseest.m</code>).</p>

</details>

<details><summary><b>Should the input data already be denoised?</b></summary><br> <!-- FAQ -->

<p align="justify">The input data consists of voxelwise/vertexwise BOLD signal, which you can already preprocesss according to your favorite recipe; however, the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox also provides the following denoising steps implemented in the <abbr title="statistical parametric mapping">SPM</abbr> plugin: <ul><li> 
 nuisance variable regression;</li>
<li> polynomial detrending;</li>
<li> band-pass filter (e.g. in the 0.01 - 0.1 Hz interval);</li>
<li> despiking. </li></ul></p>

<p align="justify">It is also possible to use a <code>temporal mask</code> to exclude some time points using the <code>Temporal mask for event detection</code> included in the <abbr title="statistical parametric mapping">SPM</abbr> plugin.</p>
<!-- temporal_mask: generated from scrubbing. -->

</details>

# Analysis:
<details><summary><b>Whole-brain or ROI analysis?</b></summary><br> <!-- FAQ -->

<p align="justify">The <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox consists of two main analysis options: 1) <i><abbr title="resting-state hemodynamic response function">rsHRF</abbr> retrieval and deconvolution</i> and 2) <i><abbr title="resting-state hemodynamic response function">rsHRF</abbr> connectivity analysis</i>. Both analyses are supported on either the whole-brain (i.e. <code>Voxels</code>/<code>Vertices</code> button in the <abbr title="statistical parametric mapping">SPM</abbr> <abbr title="graphical user interface">GUI</abbr>) or ROI (i.e. <code>ROIs-volume</code>/<code>ROIs-surface</code> button in the <abbr title="statistical parametric mapping">SPM</abbr> <abbr title="graphical user interface">GUI</abbr>) level. However, outlier removal is only legit when conducting a whole-brain analysis. 
 <!-- Both analyses can be performed  However, outlier removal denoted by OMrl (see output example) is only legit when conducting whole brain analysis. rshrf retrieval and deconvl is available both in the matlab standalone as well as in the spm plugin; however connectivity analyis is currently only avalibale in the spm plugin. Here below, you can find an outline for workflow examples for both the standalone and the spm plugin.</p>-->

</details>

## Examples:
__REMARK__: Examples for <i><abbr title="resting-state hemodynamic response function">rsHRF</abbr> connectivity analysis</i> (functional connectivity: Pearson/Spearman correlation, Pairwise/Conditional/Partially Conditioned Granger causality) are only provided for the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> <abbr title="statistical parametric mapping">SPM</abbr> plugin.

<!-- for both the matlab standalone as well as the spm plugin, you can find a few examples; 
matlab standalone; you can find the main scripts in; these will use the subfunctions (scripts) provide in the code folder
the data needed, can be found; three different data types are used (remove bilgin?); one voxel of the a human connectome participant is used (which one?); which again demonstrates one of the values, the voxel-wise/vertex-wise level of the script. Vertex-wise exaple?
The examples are demonstrated with the five different HRF basis functions: compare all five of them - same? document them; 
-- standalone: script demo: test all + input&output
-- SPM plugin: slide; demo batches
SPM plugin
The script spm_rsHRF.m is the main one, and it calls rsHRF.m. These two files are specific to the SPM plugin. 
See [rsHRF_toolbox.pptx](https://github.com/guorongwu/rsHRF/raw/master/rsHRF_toolbox.pptx) for more details (Installation/Usage/Outputs).
![rsHRF GUI](https://github.com/guorongwu/rsHRF_data/raw/master/rsHRF_GUI.png)
    -- two videos
    -- can be visualized: how? one video; plus image batch
    -- batch demo (in .zip) - slide X till X
--> 

<!-- F9 for matlab to run -->
