<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

⚙️ Installation
----

# Dependencies
<p align="justify">In order to successfully run the MATLAB code for <i>rsHRF deconvolution and connectivity analysis</i>, <a href="https://nl.mathworks.com/help/install/"><b>MATLAB</b></a> and <a href="https://www.fil.ion.ucl.ac.uk/spm/software/download/"><abbr title="Statistical Parametric Mapping"><b>SPM</b></abbr></a> should be installed. <abbr title="Statistical Parametric Mapping">SPM</abbr> is still necessary because the MATLAB code uses some of its basis functions (e.g. <code>spm_vol.m</code>, <code>spm_read_vols.m</code>...). After completing the installation, open MATLAB and add the <abbr title="Statistical Parametric Mapping">SPM</abbr> directory to the MATLAB search path (<code>addpath.m</code>).</p>

# Download the release version of your choice
<p align="justify">Click on the corresponding 🏷 to download one of the release versions (i.e. <a href="#v1.0">v1.0</a>, <a href="#v2.0">v2.0</a>, and <a href="#v2.2">v2.2</a>) of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> <abbr title="statistical parametric mapping">SPM</abbr> plugin as a <i>.zip</i> folder in <code>Downloads</code>. For each release version, the main modifications are listed. For more information, head over to the HISTORY PAGE! <!-- The <a href="https://github.com/compneuro-da/rsHRF"><abbr title="resting-state hemodynamic response function">rsHRF</abbr> GitHub repository</a> will always contain the latest version of the <abbr title="statistical parametric mapping">SPM</abbr> plugin (<i>Jan 9, 2019</i>: <b>v2.0</b>).--></p> 

<b>Release version</b>: 

<details><summary><i>rsHRF v2.2</i> <a href="">🏷</a> <a href="">📅</a></summary>
<br>

```diff
!  Main modifications (M):  
``` 

<!--
* <p align="justify"><b>sFIR</b>: </p>
* <p align="justify"><b>regularization</b>: </p> -->
<br>

</details>

<details><summary><i>rsHRF v2.0</i> <a href="https://github.com/compneuro-da/rsHRF/archive/v2.0.zip">🏷</a> <a href="https://github.com/sofievdbos/rsHRF/wiki/01.-History-and-Development:-MATLAB-(standalone-and-SPM-plugin)#v2.0">📅</a></summary>
<br>

```diff
!  Main modifications (M):  
``` 

<!--
Two types of connectivity analyses have been added to the processing pipeline: 

* <p align="justify"><b>functional connectivity</b>: functional connectivity analyses have been added to the processing pipeline, including the Pearson and Spearman correlation.</p>
* <p align="justify"><b>effective connectivity</b>: effective connectivity analyses have been added to the processing pipeline; more specifically the Pairwise/Conditional/Partially Conditioned Granger causality methods.</p> -->
<br>

</details>

<details><summary><i>rsHRF v1.0</i> <a href="https://github.com/compneuro-da/rsHRF/archive/v1.0.zip">🏷</a> <a href="https://github.com/sofievdbos/rsHRF/wiki/01.-History-and-Development:-MATLAB-(standalone-and-SPM-plugin)#v1.0">📅</a></summary>
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

# Setup: install and open as an SPM toolbox

Ater downloading the release version of your choice, you can either use the scripts in ´code´ as a standalone, in combination with a demo script you will find in ´documentation/demo´ or as an SPM plugin with the GUI. The demo script can also be tranformed for usage on the cluster (see: XXX). However, for now, the connectivity scripts are not yet available as a standalone as they are incorporated in the SPM script. <!-- use separet functions (for cluster usage); connectivity only available as part of SPM plugin --> Follow the instructions. For more information, look at the ppt or watch the narrated video!

<img align="right" src="https://github.com/compneuro-da/rsHRF/blob/update/img/install_02.png" alt="Download" width="250"/>

1. <p align="justify"> <b>Extract the <i>.zip</i> folder and move it to <code>./SPM/toolbox/rsHRF</code>.</b> <br> After downloading a specific <abbr title="statistical parametric mapping">SPM</abbr> plugin <a title="release version">🏷</a> as a <i>.zip</i> folder in <code>Downloads</code>, you first have to extract the <i>.zip</i> folder and move the extracted folder, along with its content to <code>./SPM/toolbox/rsHRF</code>. You can do this by running <a href="https://github.com/compneuro-da/rsHRF/blob/master/rsHRF_install_SPM.m" title="rsHRF_install_SPM.m"><code>rsHRF_install_SPM.m</code></a> in the Command Window opened within the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> folder. 

2. <p align="justify"> <b>Open the rsHRF GUI.</b> <br> There are two separate GUIs: one for <i>rsHRF retrieval and deconvolution</i> and one for <i>rsHRF connectivity analysis</i>. In order to open them, type respectively <code>rsHRF</code> and <code>rsHRF conn</code> in the MATLAB Command Window. You can also open them using the SPM GUI. To do so, go to the SPM folder and type spm fmri in the Command Window. Once the SPM GUI has appeared, you go to the drop-down toolbox option under the display button where you select the rsHRF toolbox. By doing so, the rsHRF GUI will open.
	
<img align="center" src="https://github.com/compneuro-da/rsHRF/blob/update/img/rsHRF_GUI_v2.2.png?raw=true" alt="Download" width="600"/>

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
