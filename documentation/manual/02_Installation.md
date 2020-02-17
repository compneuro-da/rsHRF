<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

⚙️ Installation
----
__IMPORTANT__: Please use Google Chrome to browse the _Installation_ Page! If you don't, you won't be able to expand the collapsibles. [<abbr title="Work In Progress"><i>WIP</i></abbr>]

# Dependencies
<p align="justify">To successfully run the MATLAB code for <i>rsHRF deconvolution and connectivity analysis</i>, <a href="https://nl.mathworks.com/help/install/"><b>MATLAB</b></a> and <a href="https://www.fil.ion.ucl.ac.uk/spm/software/download/"><abbr title="Statistical Parametric Mapping"><b>SPM</b></abbr></a> should be installed. <abbr title="Statistical Parametric Mapping">SPM</abbr> is still necessary because the MATLAB code uses some of its basis functions (e.g. <code>spm_vol.m</code>, <code>spm_read_vols.m</code>...). After completing the installation, open MATLAB and add the <abbr title="Statistical Parametric Mapping">SPM</abbr> directory including all its folders and subfolders to the MATLAB search path:</p>
	
``` matlab
spmdir = '/full/path/to/spm12/';
addpath(genpath(spmdir));
 ```

# Download the release version of your choice
<p align="justify">To download the release version of your choice (i.e. either v1.0, v2.0, or v2.2) as a <i>rsHRF.zip</i> folder in your <code>Downloads</code> folder, right-click on the corresponding 🏷 here below and open the link in a new tab. For each release version, the main <i>code</i> modifications are listed; open the collapsibles to have a closer look. For more information, head over to the HISTORY PAGE!</p> <!-- The <a href="https://github.com/compneuro-da/rsHRF"><abbr title="resting-state hemodynamic response function">rsHRF</abbr> GitHub repository</a> will always contain the latest version of the <abbr title="statistical parametric mapping">SPM</abbr> plugin (<i>Jan 9, 2019</i>: <b>v2.0</b>).--> 

<b>Release version</b>: 

<details><summary><i>rsHRF v2.2</i> <a href="">🏷</a> </summary>
<br>

```diff
!  Main modifications (M):  
``` 

* <p align="justify"><b>surface-based analysis</b>: a surface-based analysis module has been added to the processing pipeline which you can select by clicking on either the <code>vertices</code> (whole-brain analysis) or <code>ROI-surface</code> (ROI analyis) panel in the GUI.</p>
* <p align="justify"><b>visualization of <abbr title="resting-state hemodynamic response function">rsHRF</abbr> shapes</b>: the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> shapes can be visualized by clicking on the <code>Display</code> panel in the GUI. The underlying script (<a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_viewer.m"><code>rsHRF_viewer.m</code></a>) has been added to the <code>code</code> folder.</p>
* <b><abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation method</b>: the HRF basis functions have been updated, i.e. a Gamma/Fourier basis function (<a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_temporal_basis.m"><code>rsHRF_estimation_temporal_basis.m</code></a>) and a m-file (<a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_impulseest.m"><code>rsHRF_estimation_impulseest.m</code></a>) for non-parametric impulse response estimation (which is not included in the rsHRF GUI) have been added, along with an an update of the (s)FIR model (<a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_FIR.m"><code>rsHRF_estimation_FIR.m</code></a>).</p> 
* <b>connectivity analysis</b>: a m-file (<a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_mvgc.m"><code>rsHRF_mvgc.m</code></a>) for multivariate Granger causality analysis  has been added to the processing pipeline.</p>
<br>

</details>

<details><summary><i>rsHRF v2.0</i> <a href="https://github.com/compneuro-da/rsHRF/archive/v2.0.zip">🏷</a></summary>
<br>

```diff
!  Main modifications (M):  
``` 

* <p align="justify"><b>functional connectivity</b>: a functional connectivity analysis module has been added to the processing pipeline, including the Pearson and Spearman correlation.</p>
* <p align="justify"><b>effective connectivity</b>: an effective connectivity analysis module has been added to the processing pipeline, including the Pairwise/Conditional/Partially Conditioned Granger causality methods.</p>
* <p align="justify"><b>rsHRF_install_SPM.m</b>: <a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_install_SPM.m"><code>rsHRF_install_SPM.m</code></a> has been added to the <code>code</code> folder to facilitate the installation of the <abbr title="statistical parametric mapping">SPM</abbr> plugin.</p>
<br>

</details>

<details><summary><i>rsHRF v1.0</i> <a href="https://github.com/compneuro-da/rsHRF/archive/v1.0.zip">🏷</a></summary>
<br>

```diff
!  Main modifications (M):  
``` 

* <p align="justify"><b>outlier removal</b>: outliers based on the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> <abbr title="response height">RH</abbr> can be deleted and interpolated accordingly by respectively using <a href="https://github.com/compneuro-da/rsHRF/blob/master/deleteoutliers.m"><code>deleteoutliers.m</code></a> and <a href="https://github.com/compneuro-da/rsHRF/blob/master/inpaint_nans3.m"><code>inpaint_nans3.m</code></a>; the output files will then contain the <abbr title="OutLier ReMoval"><i>Olrm</i></abbr> abbreviation. Outlier removal is only legit when conducting a whole-brain analysis.</p>
* <p align="justify"><b>local peak detection</b>: the parameter used for local peak detection (<code>localK</code>) has been modified with its value depending on the <abbr title="repetition time">TR</abbr>.</p>
* <p align="justify"><b>global parameter modification</b>: some global parameters such as the interpolation method for outlier removal, can be adapted in <a href="https://github.com/compneuro-da/rsHRF/blob/master/wgr_rsHRF_global_para.m"><code>wgr_rsHRF_global_para.m</code></a>.</p>
* <p align="justify"><b><abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation method</b>: the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation method can be set to either <a href="https://github.com/compneuro-da/rsHRF/blob/master/wgr_rshrf_estimation_canonhrf2dd_par2.m"><abbr title="canonical HRF with its delay and dispersion derivatives"><i>canon2dd</i></abbr></a> or <a href="https://github.com/compneuro-da/rsHRF/blob/master/wgr_rsHRF_FIR.m"><abbr title="smoothed Finite Impulse Response basis functions"><i>(s)FIR</i></abbr></a>.</p>
<br>

</details>

# Setup
## Install the SPM toolbox
<p align="justify">After downloading the release version of your choice, you can either choose to use the scripts in the <code>rsHRF/code/</code> folder as a MATLAB Standalone or to install the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox as an <abbr title="statistical parametric mapping">SPM</abbr> plugin.</p>

<p align="justify"><b>MATLAB Standalone.</b> The MATLAB Standalone is currently only available for voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation and deconvolution as well as parameter retrieval (<abbr title="response height">RH</abbr>, <abbr title="time to peak">TTP</abbr>, and <abbr title="full width at half maximum">FWHM</abbr>). Demo scripts to use the functions in the <code>rsHRF/code/</code> folder can be found in the <code>rsHRF/documentation/demo/</code> folder. The demo scripts can easly be adapted into a function which you can use on your High Performance Computer. [<abbr title="Work In Progress"><i>WIP</i></abbr>] [EXAMPLE] </p>
	
<!-- The demo script can also be tranformed for usage on the cluster (see: XXX). However, for now, the connectivity scripts are not yet available as a standalone as they are incorporated in the SPM script.  use separet functions (for cluster usage); connectivity only available as part of SPM plugin  Follow the instructions. For more information, look at the ppt or watch the narrated video! -->

<!-- <img align="right" src="https://github.com/compneuro-da/rsHRF/blob/update/img/install_02.png" alt="Download" width="200"/> -->

<p align="justify"><b><abbr title="statistical parametric mapping">SPM</abbr> plugin.</b> If you want to install the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox as an <abbr title="statistical parametric mapping">SPM</abbr> plugin, extract the <i>rsHRF.zip</i> folder in your <code>Downloads</code> folder and move its contents into <code>./SPM/toolbox/rsHRF/</code>. To do so, open the downloaded <abbr title="resting-state hemodynamic response function"><code>rsHRF</code></abbr> folder and run the <a href="https://github.com/compneuro-da/rsHRF/blob/master/rsHRF_install_SPM.m" title="rsHRF_install_SPM.m"><code>rsHRF_install_SPM.m</code></a> script in the MATLAB Command Window. All scripts within the downloaded <abbr title="resting-state hemodynamic response function"><code>rsHRF</code></abbr> folder will be copied into a folder named <abbr title="resting-state hemodynamic response function"><code>rsHRF</code></abbr> located in the <code>spm/toolbox/</code> folder.</p>  

## Open the rsHRF GUI

2. <p align="justify"> <b>Open the rsHRF GUI.</b> <br> There are two separate GUIs: one for <i>rsHRF retrieval and deconvolution</i> and one for <i>rsHRF connectivity analysis</i>. In order to open them, type respectively <code>rsHRF</code> and <code>rsHRF conn</code> in the MATLAB Command Window. You can also open them using the SPM GUI. To do so, go to the SPM folder and type spm fmri in the MATLAB Command Window. Once the SPM GUI has appeared, you go to the drop-down toolbox option under the display button where you select the rsHRF toolbox. By doing so, the rsHRF GUI will open.
	
<img align="center" src="https://github.com/compneuro-da/rsHRF/blob/update/img/rsHRF_GUI_v2.2.png?raw=true" alt="Download" width="600"/>

<div style="width:100%; padding-bottom:56.25%; position:relative;">

  <iframe src="https://github.com/compneuro-da/rsHRF/blob/update/img/installation/tutorial.html" style="position:absolute; top:0px; left:0px; width:100%; height:100%; border: none; overflow: hidden;"></iframe>

</div>

<!--
1. Extract all
2. to spm tooloxes
3. open matlab
4. in command line
which spm version?? 
    * Run code in Command Window (within <abbr title="resting-state hemodynamic response function">rsHRF</abbr> folder)
      >> rsHRF_install_SPM
       -- SPM should be installed
       -- Codes will be copied to ./SPM/toolbox/rsHRF -- link to folder?? 
			             = spm('Dir')
    * Or add it into ./SPM/toolbox/
3. Start <abbr title="resting-state hemodynamic response function">rsHRF</abbr>
4. Start connectivity analysis
-->

<!--
```matlab
rsHRF_install_SPM %haha
```
```xml
<myxml>
   <someElement />  
</myxml>
```
-->

[![alt text](https://github.com/compneuro-da/rsHRF/blob/update/img/example_hrf.png)](https://github.com/compneuro-da/rsHRF/blob/update/img/installation/rsHRF01_install_SPM_plugin_200217.mp4 "title")
