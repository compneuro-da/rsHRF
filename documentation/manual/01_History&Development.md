<a name="top">
<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

📅 GitHub Commit History: MATLAB (standalone + <abbr title="Statistical Parametric Mapping">SPM</abbr> plugin) 
---- 
__IMPORTANT__: Please use Google Chrome to browse the _History and Development_ Page! If you don't, you won't be able to expand the collapsibles. [<abbr title="Work In Progress"><i>WIP</i></abbr>]

# Introduction
<img align="right" src="https://user-images.githubusercontent.com/23309041/73989536-66d51380-4946-11ea-8d2a-e1aa0d5f2eec.png" alt="SPM_plugin_history" width="400"/>

__REMARK__: We use <a href="https://semver.org/">SemVer</a> to define the different release versions. For the currently available release versions, head over to the <a href="https://github.com/compneuro-da/rsHRF/releases">tags on this GitHub repository</a>.

<p align="justify">Here below, you can find a chronological overview [newest --> oldest] of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox's main modifications through the last years [2019 --> 2016], resulting in the latest version of the <abbr title="Statistical Parametric Mapping"></abbr> plugin (<b>v2.2</b>). The three main releases of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> <abbr title="Statistical Parametric Mapping">SPM</abbr> plugin (i.e. <a href="#v1.0">v1.0</a>, <a href="#v2.0">v2.0</a>, and <a href="#v2.2">v2.2</a>) can be downloaded by right-clicking on the corresponding 🏷 next to the commit date and opening the link in a new tab; these collapsibles are opened by default. For more information, head over to the INSTALLATION PAGE!</p> 
<p align="justify">Modifications are either related to code or documentation/demo files which will be indicated by respectively a 💻 or 📖 emoji next to the commit date. Please note that a few smaller commits (e.g. minor debugging, lay-out tweaks...) are omitted. To consult the complete overview of a file's commit history, right-click on the corresponding superscript number and choose to open the link in a new tab. For example, to consult the GitHub commit history of <code>wgr_get_parameters.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_get_parameters.m" title="wgr_get_parameters.m">2</a></sup>, right-click on the <sup><i>2</i></sup> superscript and open the link in a new tab. You will now see a list of two commits associated with the <code>wgr_get_parameters.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_get_parameters.m" title="wgr_get_parameters.m">2</a></sup> script, i.e. respectively on <i>Oct 28, 2016</i> and <i>Jul 20, 2016</i>. After clicking on a commit (e.g. <abbr title="canonical HRF with its delay and dispersion derivatives"><i>canon2dd</i></abbr> on <i>Jul 20, 2016</i>), you will be able to see the number of changed files, along with the number of additions and deletions. In this particular example, there are 2 changed files with 58 additions and 4 deletions.</p> 

<!-- <img align="center" src="https://user-images.githubusercontent.com/23309041/74041807-a3485400-49c6-11ea-9038-10eedb805678.png" alt="commit_history_01" width="800"/> -->

# Legend

| Emoji        | Meaning            |  
| -------------|:-------------------| 
| 💻          | CODE               | 
| 📖          | DOCUMENTATION/DEMO |  
| 🔑          | LICENSE            |    
| 🏷           | RELEASE VERSION    |

# Overview

Files with:

```diff
+  Additions (A)
!  Modifications (M)
-  Deletions (D)
``` 

### 2019

<details open><summary><a href="#top">🔝</a> <a name="v2.2"><abbr title="SPM plugin (v2.2)"><i>Nov 15, 2019</i></abbr> 💻 📖 – <a href="">🏷</a> v2.2</summary>
<br>

```diff
+  💻 CODE + 📖 DOCUMENTATION/DEMO: SPM plugin (v2.2)
``` 

<p align="justify"><b>Version 2.2</b> of the <abbr title="Statistical Parametric Mapping">SPM</abbr> plugin is released on GitHub, adding a surface-based analysis module (<i>Aug, 2018</i>: <b>v2.1</b>), changing the <abbr title="graphical user interface">GUI</abbr> (with the addition of a surface analysis panel and <code>Display</code>), adding a <code>rsHRF_viewer.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_viewer.m" title="rsHRF_viewer.m">36</a></sup> for the visualization of <abbr title="hemodynamic response function">HRF</abbr> shapes, adding a m-file (<code>rsHRF_mvgc.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_mvgc.m" title="rsHRF_mvgc.m">35</a></sup>) for multivariate Granger causality connectivity analysis, updating the <abbr title="hemodynamic response function">HRF</abbr> basis functions with the addition of Gamma/Fourier basis functions<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_estimation_temporal_basis.m" title="rsHRF_estimation_temporal_basis.m">32</a></sup> (which are more flexible and support a finer temporal grid), updating the <abbr title="(smoothed) finite impulse response"><i>(s)FIR</i></abbr> model using <abbr title="kth-order autoregression">AR(k)</abbr> for auto-correlated noise modeling<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_estimation_FIR.m" title="rsHRF_estimation_FIR.m">33</a></sup>, and adding a m-file (<code>rsHRF_estimation_impulseest.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_estimation_impulseest.m" title="rsHRF_estimation_impulseest.m">34</a></sup>, see code for help) for non-parametric impulse response estimation (which is not included in the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> <abbr title="graphical user interface">GUI</abbr>). These updates are also available in: <code>update_log.txt</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/update_log.txt" title="update_log.txt">37</a></sup>. Corresponding modifications are inserted into <code>rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.m" title="rsHRF.m">20</a></sup>, <code>tbx_cfg_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/tbx_cfg_rsHRF.m" title="tbx_cfg_rsHRF.m">21</a></sup>, and <code>rsHRF_install_SPM.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_install_SPM.m" title="rsHRF_install_SPM.m">31</a></sup>. A few global parameters are added to <code>wgr_rsHRF_global_para.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_global_para.m" title="wgr_rsHRF_global_para.m">25</a></sup>. The pipeline for voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation and deconvolution as well as parameter retrieval using the updated <abbr title="hemodynamic response function">HRF</abbr> basis functions is illustrated by means of three separate hands-on demos in MATLAB: <code>demo_rsHRF_temporal_basis.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_code/demo_rsHRF_temporal_basis.m" title="demo_rsHRF_temporal_basis.m">38</a></sup>, <code>demo_rsHRF_FIR_sFIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_code/demo_rsHRF_FIR_sFIR.m" title="demo_rsHRF_FIR_sFIR.m">39</a></sup>, and <code>demo_rsHRF_impulseest.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_code/demo_rsHRF_impulseest.m" title="demo_rsHRF_impulseest.m">40</a></sup>.</p>

</details>

<details><summary><i>Apr 13, 2019</i> 📖</summary>
<br>

```diff
+  📖 DOCUMENTATION/DEMO 
``` 

<p align="justify">The pipeline for voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> deconvolution using <code>wgr_deconv_canonhrf_par.m</code><sup><a href="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m" title="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m">4</a></sup>, along with a brief theoretical framework, is available on <a href="https://guorongwu.github.io/HRF/rsHRF_deconv.html">https://guorongwu.github.io/HRF/rsHRF_deconv.html</a>.</p>

</details>

<details open><summary><a href="#top">🔝</a> <a name="v2.0"><abbr title="SPM plugin (v2.0)"><i>Jan 9, 2019</i></abbr> 💻 📖 – <a href="https://github.com/compneuro-da/rsHRF/archive/v2.0.zip">🏷</a> v2.0</summary>
<br> 

```diff
!  💻 CODE: SPM plugin (v2.0)
``` 

<p align="justify"><b>Version 2.0</b> of the <abbr title="Statistical Parametric Mapping">SPM</abbr> plugin is released on GitHub, adding both functional (Pearson/Spearman correlation) and effective (Pairwise/Conditional/Partially Conditioned Granger causality) connectivity analyses to the processing pipeline. Corresponding modifications are inserted into <code>rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.m" title="rsHRF.m">20</a></sup>, <code>rsHRF.man</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.man" title="rsHRF.man">22</a></sup>, and <code>tbx_cfg_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/tbx_cfg_rsHRF.m" title="tbx_cfg_rsHRF.m">21</a></sup> with all connectivity subfunctions being incorporated into <code>rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.m" title="rsHRF.m">20</a></sup>. A few global parameters are added to/omitted from <code>wgr_rsHRF_global_para.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_global_para.m" title="wgr_rsHRF_global_para.m">25</a></sup>.</p>

<details><summary></summary>

```diff
+  💻 CODE: rsHRF_install_SPM.m
``` 

<p align="justify">Running <code>rsHRF_install_SPM.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_install_SPM.m" title="rsHRF_install_SPM.m">31</a></sup> in the MATLAB Command Window within the downloaded <abbr title="resting-state hemodynamic response function">rsHRF</abbr> GitHub repository will copy all files included in the <code>rsHRF_file</code> cell array to <code>path/to/spm12/toolbox/rsHRF</code>.</p>

</details>

<details><summary></summary>

```diff
!  📖 DOCUMENTATION/DEMO
``` 

<p align="justify">The folder <code>demo_jobs</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_jobs.zip" title="demo_jobs.zip">26</a></sup> is expanded and compressed (<i>.zip</i>) by adding several new MATLAB batch jobs containing connectivity analyses to illustrate the pipeline for voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation and deconvolution as well as parameter retrieval using the <abbr title="Statistical Parametric Mapping">SPM</abbr> plugin (v2.0). Corresponding information concerning the new batch job examples can be found in an updated version of <code>rsHRF_toolbox.pptx</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_toolbox.pptx" title="rsHRF_toolbox.pptx">27</a></sup>.</p> 

</details>
</details>

### 2018

<details><summary><abbr title="<DEBUGGING>"><i>Aug 29 + Nov 4 + Dec 9, 2018</i></abbr> 💻 📖</summary>
<br>

```diff
!  💻 CODE + 📖 DOCUMENTATION/DEMO: <DEBUGGING>
``` 

* <code><b>knee_pt.m</b></code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/knee_pt.m" title="knee_pt.m">12</a></sup>: 
    * <i>Nov 4, 2018</i>: If <code>(length(lag) < 3)</code>, <code>min()</code> is used instead of the <code>knee_pt.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/knee_pt.m" title="knee_pt.m">12</a></sup> function for the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> lag estimation; the <code>knee_pt.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/knee_pt.m" title="knee_pt.m">12</a></sup> function is modified accordingly. 
    * <i>Dec 9, 2018</i>: To uniform <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> and <code>wgr_rsHRF_FIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_FIR.m" title="wgr_rsHRF_FIR.m">11</a></sup>, the <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> function is revised. Concretely, <code>beta(:, id+1)</code> and <code>lag(id+1)</code> are selected instead of <code>beta(:, id)</code> and <code>lag(id)</code> respectively. 
* <b><i>datatype</i></b>: 
    * <i>Nov 4, 2018</i>: <code>v1(i).dt = [16,0]</code> is included as the NIfTI file’s (.nii) datatype in <code>demo_4d_data.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_4D_data.m" title="demo_4d_rsHRF.m">7</a></sup> and <code>rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.m" title="rsHRF.m">20</a></sup>. 

</details>

<details><summary><abbr title="LICENSE"><i>Aug 14, 2018</i></abbr> 🔑</summary>
<br>

```diff
+  🔑 LICENSE
``` 

<p align="justify">A BSD 3-Clause <code><b>License</b></code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/LICENSE" title="LICENSE">30</a></sup> is added to the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> GitHub repository.</p>

</details>

<details><summary><abbr title="Tikhonov regularization"><i>Aug 4 – 5, 2018</i></abbr> 💻 📖</summary>
<br>

```diff
!  💻 CODE + 📖 DOCUMENTATION/DEMO: Tikhonov regularization
``` 

<p align="justify">The <code>knee_pt.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/knee_pt.m" title="knee_pt.m">12</a></sup> function is inserted into <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> as an alternative to <code>min()</code> for the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> lag estimation. To avoid division by zero when deconvolving the BOLD signal using the Inverse fast Fourier transform (<code>ifft.m</code>), the <abbr title="data_deconv = ifft(conj(H).*M./(H.*conj(H)+.1*mean(H.*conj(H))));"><b><i>Tikhonov regularization</i></b></abbr> was implemented and adapted accordingly in <code>demo_4d_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_4D_data.m" title="demo_4d_rsHRF.m">7</a></sup>, <code>demo_voxel.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_voxel.m" title="demo_voxel.m">15</a></sup>, <code>demo_voxel_calcium.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_voxel_calcium.m" title="demo_voxel_calcium.m">17</a></sup>, and <code>rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.m" title="rsHRF.m">20</a></sup>.</p> 

</details>

<details open><summary><a href="#top">🔝</a> <a name="v1.0"><abbr title="SPM plugin (v1.0) + batch jobs"><i>Jul 21 – Jul 31, 2018</i></abbr> 💻 📖 – <a href="https://github.com/compneuro-da/rsHRF/archive/v1.0.zip">🏷</a> v1.0</summary> 

<br>

```diff
+  💻 CODE: SPM plugin (v1.0)
``` 

<p align="justify"><b>Version <b><i>1.0</i></b></a></b> of the <abbr title="Statistical Parametric Mapping">SPM</abbr> plugin is committed to GitHub with the main <code>spm_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/spm_rsHRF.m" title="spm_rsHRF.m">19</a></sup> script calling <code>rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.m" title="rsHRF.m">20</a></sup> along with its configuration file (<code>tbx_cfg_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/tbx_cfg_rsHRF.m" title="tbx_cfg_rsHRF.m">21</a></sup>). The <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox’s version and description can be found in <code>rsHRF.man</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.man" title="rsHRF.man">22</a></sup>. Henceforth, <b><i>outliers</i></b> based on the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> <abbr title="response height">RH</abbr> can be deleted and interpolated accordingly by respectively using <code>deleteoutliers.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/deleteoutliers.m" title="deleteoutliers.m">23</a></sup> and <code>inpaint_nans3.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/inpaint_nans3.m" title="inpaint_nans3.m">24</a></sup>. The parameter used for local peak detection (<code>localK</code>) is modified with its value depending on the <abbr title="repetition time">TR</abbr>. Some <b><i>global parameters</i></b> such as the interpolation method, can be adapted in <code>wgr_rsHRF_global_para.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_global_para.m" title="wgr_rsHRF_global_para.m">25</a></sup>, while the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation method can be set to either <abbr title="canonical HRF with its delay and dispersion derivatives"><i>canon2dd</i></abbr><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> or <abbr title="(smoothed) finite impulse response"><i>(s)FIR</i></abbr><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_FIR.m" title="wgr_rsHRF_FIR.m">11</a></sup>.</p> 

<details><summary></summary>

```diff
+  📖 DOCUMENTATION/DEMO: batch jobs
``` 

<p align="justify">The folder <code>demo_jobs</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_jobs.zip" title="demo_jobs.zip">26</a></sup> is created containing a few different MATLAB <b><i>batch jobs</i></b> to illustrate the pipeline for voxel-wise resting-state <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation and deconvolution as well as parameter retrieval using the <abbr title="Statistical Parametric Mapping">SPM</abbr> plugin (v1.0). More detailed information concerning the installation and the batch job examples can be found in <code>rsHRF_toolbox.pptx</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_toolbox.pptx" title="rsHRF_toolbox.pptx">27</a></sup>. In addition, two more data examples are available: <code>eve.mat</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/eve.mat" title="eve.mat">28</a></sup>, and <code>voxelsample_bilgin.mat</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/voxelsample_bilgin.mat" title="voxelsample_bilgin.mat">29</a></sup>.</p>

</details>
</details>

<details><summary><abbr title="zscoring + filtering + examples"><i>Jun 18, 2018</i></abbr> 📖</summary>
<br>

```diff
!  📖 DOCUMENTATION/DEMO: zscoring and filtering
``` 

<p align="justify">The hands-on demo in MATLAB<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_4D_data.m" title="demo_4d_rsHRF.m">7</a></sup> is expanded by explicitly including <b><i>zscoring</i></b> and <b><i>filtering</i></b><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rest_IdealFilter.m" title="rest_IdealFilter.m">13</a></sup> into the preprocessing pipeline. You should not forget to include these steps. The <code>README.md</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/README.md" title="README.md">5</a></sup> file is adapted accordingly; but these and further updates are not transferred to <a href="https://guorongwu.github.io/HRF/">https://guorongwu.github.io/HRF/</a>.</p>

```diff
+  📖 DOCUMENTATION/DEMO: examples
``` 

<p align="justify">The expanded demo applied to two specific <b><i>examples</i></b>, i.e., a sample voxel from the Human Connectome Project (HCP; <code>voxelsample_hcp.mat</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/voxelsample_hcp.mat" title="voxelsample_hcp.mat">16</a></sup>; XXX) and <code>calcium.mat</code><sup><a href="https://github.com/compneuro-da/rsHRF_data/commits/master/calcium.mat" title="calcium.mat">18</a></sup>, is committed to GitHub.</p> 

</details>

<details><summary><abbr title="write back the number of pseudo point process events"><i>Jun 7 – 8, 2018</i></abbr> 💻 📖</summary>
<br>

```diff
!  💻 CODE + 📖 DOCUMENTATION/DEMO: number of pseudo point process events
``` 

<p align="justify">The event_bold parameter is added/updated in(to) the <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> and <code>wgr_rsHRF_FIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_FIR.m" title="wgr_rsHRF_FIR.m">11</a></sup> functions. As a result, each voxel’s <b><i>number of pseudo point process events</i></b> can be written back into a NIfTI (.nii) file (<code>demo_4d_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_4D_data.m" title="demo_4d_rsHRF.m">7</a></sup>). The <code>knee_pt.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/knee_pt.m" title="knee_pt.m">12</a></sup> function has not yet been inserted into <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> as an alternative to <code>min()</code> for the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> lag estimation. <sup>M: <i>Aug 4 – 5, 2018</i></sup></p> 

</details>

<details><summary><abbr title="<MERGE + REMOVE REDUNDANT CODE> + filtering functions"><i>May 31 – Jun 3, 2018</i></abbr> 💻 📖</summary>
<br>

```diff
+  💻 CODE: <MERGE + REMOVE REDUNDANT CODE> + filtering functions
``` 

<p align="justify">The <code>hrf_retrieval_and_deconvolution_para.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commit/5a72732b6ed9fd34987fff68ceacbcd7a0196264#diff-b182bf5beb61ef4c399c0ff0ac61c982" title="hrf_retrieval_and_deconvolution_para.m">9</a></sup> script is removed, along with the <i>rbeta</i> <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation method. The function is replaced by a standalone MATLAB function to estimate the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> using the <abbr title="(smoothed) finite impulse response"><i>(s)FIR</i></abbr> basis functions<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_FIR.m" title="wgr_rsHRF_FIR.m">11</a></sup> without the option to upsample the time resolution. The temporal mask used to exclude pseudo point process events induced by motion artifacts is included as a parameter. In addition, the <code>knee_pt.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/knee_pt.m" title="knee_pt.m">12</a></sup> function is committed to GitHub and inserted into <code>wgr_rsHRF_FIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_FIR.m" title="wgr_rsHRF_FIR.m">11</a></sup> as an alternative to <code>min()</code> for the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> lag estimation. Two <b><i>filtering functions</i></b><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rest_IdealFilter.m" title="rest_IdealFilter.m">13</a></sup><sup>,</sup><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rest_nextpow2_one35.m" title="rest_nextpow2_one35.m">14</a></sup> from the REST toolbox (<a href="#ref3">Song et al., 2011</a>) are committed to GitHub as well.</p>

```diff
!  📖 DOCUMENTATION/DEMO: <MERGE + REMOVE REDUNDANT CODE>
``` 

<p align="justify">The <code>demo_main_deconvolution_FIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commit/a422be9d1e5bc406fbe8ad37eac25c5529c5c159" title="demo_main_deconvolution_FIR.m">10</a></sup> script is removed. Instead, the <abbr title="finite impulse response"><i>FIR</i></abbr> and <abbr title="(smoothed) finite impulse response"><i>(s)FIR</i></abbr> <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation methods (<code>wgr_rsHRF_FIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_FIR.m" title="wgr_rsHRF_FIR.m">11</a></sup>) are included into <code>demo_4d_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_4D_data.m" title="demo_4d_rsHRF.m">7</a></sup>.</p>

</details>

### 2017

<details><summary><abbr title="sFIR and rbeta rsHRF estimation method"><i>Dec 18, 2017</i></abbr> 💻 📖</summary>
<br>

```diff
+  💻 CODE: sFIR and rbeta
``` 

<p align="justify">An updated and expanded version<sup><a href="https://github.com/compneuro-da/rsHRF/commit/5a72732b6ed9fd34987fff68ceacbcd7a0196264#diff-b182bf5beb61ef4c399c0ff0ac61c982" title="hrf_retrieval_and_deconvolution_para.m">9</a></sup> of <code>wgr_deconv_canonhrf_par.m</code><sup><a href="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m" title="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m">4</a></sup> with the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation method set by the <i>flag</i> parameter, being either <abbr title="canonical HRF with its delay and dispersion derivatives"><i>canon2dd</i></abbr>, the smoothed Finite Impulse Response basis functions (denoted as <abbr title="smoothed finite impulse response"><b><i>sFIR</i></b></abbr>), or the <b><i>rbeta</i></b> function (<a href="#ref2">Tagliazucchi et al., 2012</a>). The option to upsample the time resolution and the <i>temporal_mask</i> argument used to exclude pseudo point process events induced by motion artifacts have not yet been included. Again, some subfunctions of the MATLAB code are modified from the <abbr title="hemodynamic response function">HRF</abbr> Estimation Toolbox<sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> (<a href="#ref1">Lindquist, Waugh, & Wager, 2007</a>).</p>

```diff
+  📖 DOCUMENTATION/DEMO
``` 

<p align="justify">A separate hands-on demo in MATLAB<sup><a href="https://github.com/compneuro-da/rsHRF/commit/a422be9d1e5bc406fbe8ad37eac25c5529c5c159" title="demo_main_deconvolution_FIR.m">10</a></sup> to illustrate the pipeline for voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> deconvolution using <code>hrf_retrieval_and_deconvolution_para.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commit/5a72732b6ed9fd34987fff68ceacbcd7a0196264#diff-b182bf5beb61ef4c399c0ff0ac61c982" title="hrf_retrieval_and_deconvolution_para.m">9</a></sup>.</p> 

</details>

<details><summary><abbr title="reading and writing 4D NIfTI files"><i>Apr 28, 2017</i></abbr> 📖</summary>
<br>

```diff
+  📖 DOCUMENTATION/DEMO: 4D NIfTI
``` 

<p align="justify">An updated and expanded version<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_4D_data.m" title="demo_4d_rsHRF.m">7</a></sup> of the more detailed hands-on demo in MATLAB including the option to read and write <b><i>4D NIfTI</i></b> (.nii) files and implementing <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> and <code>wgr_get_parameters.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_get_parameters.m" title="wgr_get_parameters.m">2</a></sup> instead of <code>wgr_deconv_canonhrf_par.m</code><sup><a href="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m" title="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m">4</a></sup>, is committed to GitHub. A data structure example<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/data_structure_example" title="data_structure_example.txt">8</a></sup> is committed as well.</p>

</details>

<details><summary><i>Dec 16, 2016</i> 📖 <sub><sup>(reference in <code>README.md</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/README.md" title="README.md">5</a></sup> on <i>Jan 13, 2017</i>)</sup></sub></summary>
<br>

```diff
!  📖 DOCUMENTATION/DEMO  
``` 

<p align="justify">The pipeline for voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation as well as parameter retrieval as illustrated by means of a concise overview in Markdown<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/README.md" title="README.md">5</a></sup>, along with a brief theoretical framework, is available on <a href="https://guorongwu.github.io/HRF/">https://guorongwu.github.io/HRF/</a>.</p>

</details>

### 2016

<details><summary><abbr title="canon2dd rsHRF estimation method"><i>Jul 20, 2016</i></abbr> 💻 📖</summary>
<br>

```diff
+  💻 CODE: canon2dd
```    

<p align="justify">MATLAB code for voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation and deconvolution<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> as well as parameter retrieval<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_get_parameters.m" title="wgr_get_parameters.m">2</a></sup> (<abbr title="response height">RH</abbr>, <abbr title="time to peak">TTP</abbr>, and <abbr title="full width at half maximum">FWHM</abbr>) using a canonical <abbr title="hemodynamic response function">HRF</abbr> with its delay and dispersion derivatives (denoted as <abbr title="canonical HRF with its delay and dispersion derivatives"><i><b>canon2dd</b></i></abbr>); also including an option to upsample the time resolution and the temporal mask used to exclude pseudo point process events induced by motion artifacts. Some subfunctions of the MATLAB code are modified from the <abbr title="hemodynamic response function">HRF</abbr> Estimation Toolbox<sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> (<a href="#ref1">Lindquist, Waugh, & Wager, 2007</a>).</p>

```diff
+  📖 DOCUMENTATION/DEMO
```
 
<p align="justify">Illustrating the pipeline for voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> estimation as well as parameter retrieval by means of a concise overview in Markdown<sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/README.md" title="README.md">5</a></sup> using <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> and <code>wgr_get_parameters.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_get_parameters.m" title="wgr_get_parameters.m">2</a></sup> and a reference to a more detailed hands-on demo in MATLAB<sup><a href="https://users.ugent.be/~dmarinaz/wgr_voxelwise_HRF_deconvolution_demo.m" title="https://users.ugent.be/~dmarinaz/wgr_voxelwise_HRF_deconvolution_demo.m">6</a></sup> using <code>wgr_deconv_canonhrf_par.m</code><sup><a href="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m" title="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m">4</a></sup>; the latter also illustrates the voxel-wise <abbr title="resting-state hemodynamic response function">rsHRF</abbr> deconvolution.</p> 

</details>

# Files

<table ><tbody ><tr></tr><tr><td><details ><summary><sub><b>Click to see more:</b></sub><br><br>

|       | File Name    | Additional information |
|-------|:-------------|:-----------------------| 
| 1.    | <code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rshrf_estimation_canonhrf2dd_par2.m" title="wgr_rshrf_estimation_canonhrf2dd_par2.m">1</a></sup> | modified from an older version (<code>wgr_deconv_canonhrf_par.m</code><sup><a href="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m" title="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m">4</a></sup>) with the modified <code>Fit_Canonical_HRF()</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> and <code>CanonicalBasisSet()</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> subfunctions; the older version<sup><a href="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m" title="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m">4</a></sup> is not committed to the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> GitHub repository |

</summary><hr>

|       | File Name    | Additional information |
|-------|:-------------|:-----------------------|
| 2.    | <code>wgr_get_parameters.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_get_parameters.m" title="wgr_get_parameters.m">2</a></sup> | modified from <code>get_parameters2.m</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> |
| 3.    | Estimation Toolbox<sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> <br> (<a href="#ref1">Lindquist, Waugh, & Wager, 2007</a>) | |   
| 4.    | <code>wgr_deconv_canonhrf_par.m</code><sup><a href="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m" title="https://users.ugent.be/~dmarinaz/wgr_deconv_canonhrf_par.m">4</a></sup> | | 
| 5.    | <code>README.md</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/README.md" title="README.md">5</a></sup> | |
| 6.    | <code>wgr_voxelwise_HRF_deconvolution_demo.m</code><sup><a href="https://users.ugent.be/~dmarinaz/wgr_voxelwise_HRF_deconvolution_demo.m" title="https://users.ugent.be/~dmarinaz/wgr_voxelwise_HRF_deconvolution_demo.m">6</a></sup> | |  
| 7.    | <code>demo_4d_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_4D_data.m" title="demo_4d_rsHRF.m">7</a></sup> | M: <i>Aug 6, 2018</i> (renamed to <code>demo_4D_data.m</code>)
| 8.    | <code>data_structure_example.txt</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/data_structure_example" title="data_structure_example.txt">8</a></sup> | M: <i>Jun 19, 2018</i> |
| 9.    | <code>hrf_retrieval_and_deconvolution_para.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commit/5a72732b6ed9fd34987fff68ceacbcd7a0196264#diff-b182bf5beb61ef4c399c0ff0ac61c982" title="hrf_retrieval_and_deconvolution_para.m">9</a></sup> | D: <i>May 31, 2018</i> (with the modified  <code>Fit_Canonical_HRF2()</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup>, <code>CanonicalBasisSet()</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup>, <code>Fit_sFIR()</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup>, and <code>tor_make_deconv_mtx3()</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> subfunctions) |
| 10.    | <code>demo_main_deconvolution_FIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commit/a422be9d1e5bc406fbe8ad37eac25c5529c5c159" title="demo_main_deconvolution_FIR.m">10</a></sup> | D: <i>May 31, 2018</i> |
| 11.    | <code>wgr_rsHRF_FIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_FIR.m" title="wgr_rsHRF_FIR.m">11</a></sup> | with the modified <code>Fit_sFIR2()</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> and <code>tor_make_deconv_mtx3()</code><sup><a href="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2" title="https://github.com/canlab/CanlabCore/tree/master/CanlabCore/HRF_Est_Toolbox2">3</a></sup> subfunctions |
| 12.    | <code>knee_pt.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/knee_pt.m" title="knee_pt.m">12</a></sup> | |  
| 13.    | <code>rest_IdealFilter.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rest_IdealFilter.m" title="rest_IdealFilter.m">13</a></sup> | |
| 14.    | <code>rest_nextpow2_one35.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rest_nextpow2_one35.m" title="rest_nextpow2_one35.m">14</a></sup> | |
| 15.    | <code>demo_voxel.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_voxel.m" title="demo_voxel.m">15</a></sup> | |
| 16.    | <code>voxelsample_hcp.mat</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/voxelsample_hcp.mat" title="voxelsample_hcp.mat">16</a></sup> | |
| 17.    | <code>demo_voxel_calcium.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_voxel_calcium.m" title="demo_voxel_calcium.m">17</a></sup> | |
| 18.    | <code>calcium.mat</code><sup><a href="https://github.com/compneuro-da/rsHRF_data/commits/master/calcium.mat" title="calcium.mat">18</a></sup> | A: <i>Jul 22, 2018</i>; M: <i>Dec 25, 2018</i> (moved to <code>rsHRF_data<code/>) |
| 19.    | <code>spm_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/spm_rsHRF.m" title="spm_rsHRF.m">19</a></sup> | |
| 20.    | <code>rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.m" title="rsHRF.m">20</a></sup> | |
| 21.    | <code>tbx_cfg_rsHRF.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/tbx_cfg_rsHRF.m" title="tbx_cfg_rsHRF.m">21</a></sup> | |
| 22.    | <code>rsHRF.man</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF.man" title="rsHRF.man">22</a></sup> | |
| 23.    | <code>deleteoutliers.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/deleteoutliers.m" title="deleteoutliers.m">23</a></sup> | |
| 24.    | <code>inpaint_nans3.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/inpaint_nans3.m" title="inpaint_nans3.m">24</a></sup> | |
| 25.    | <code>wgr_rsHRF_global_para.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/wgr_rsHRF_global_para.m" title="wgr_rsHRF_global_para.m">25</a></sup> | |
| 26.    | <code>demo_jobs</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_jobs.zip" title="demo_jobs.zip">26</a></sup> | contains <code>ROIsig_job01.mat</code>, <code>ROIwise_normalized_space_job01.mat</code>, <code>voxelwise_native_space_job01.mat</code>, <code>voxelwise_normalized_space_job01.mat</code>|
| 27.    | <code>rsHRF_toolbox.pptx</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_toolbox.pptx" title="rsHRF_toolbox.pptx">27</a></sup> | |
| 28.    | <code>eve.mat</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/eve.mat" title="eve.mat">28</a></sup> | |
| 29.    | <code>voxelsample_bilgin.mat</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/voxelsample_bilgin.mat" title="voxelsample_bilgin.mat">29</a></sup> | |
| 30.    | <code>License</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/LICENSE" title="LICENSE">30</a></sup> | |
| 31.    | <code>rsHRF_install_SPM.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_install_SPM.m" title="rsHRF_install_SPM.m">31</a></sup> | |
| 32.    | <code>rsHRF_estimation_temporal_basis.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_estimation_temporal_basis.m" title="rsHRF_estimation_temporal_basis.m">32</a></sup> | |
| 33.    | <code>rsHRF_estimation_FIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_estimation_FIR.m" title="rsHRF_estimation_FIR.m">33</a></sup> | |
| 34.    | <code>rsHRF_estimation_impulseest.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_estimation_impulseest.m" title="rsHRF_estimation_impulseest.m">34</a></sup> | |
| 35.    | <code>rsHRF_mvgc.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_mvgc.m" title="rsHRF_mvgc.m">35</a></sup> | |
| 36.    | <code>rsHRF_viewer.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/rsHRF_viewer.m" title="rsHRF_viewer.m">36</a></sup> | |
| 37.    | <code>update_log.txt</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/update_log.txt" title="update_log.txt">37</a></sup> | |
| 38.    | <code>demo_rsHRF_temporal_basis.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_code/demo_rsHRF_temporal_basis.m" title="demo_rsHRF_temporal_basis.m">38</a></sup> | |
| 39.    | <code>demo_rsHRF_FIR_sFIR.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_code/demo_rsHRF_FIR_sFIR.m" title="demo_rsHRF_FIR_sFIR.m">39</a></sup> | |
| 40.    | <code>demo_rsHRF_impulseest.m</code><sup><a href="https://github.com/compneuro-da/rsHRF/commits/master/demo_code/demo_rsHRF_impulseest.m" title="demo_rsHRF_impulseest.m">40</a></sup> | |

</details></td></tr></tbody>
</table>

# References: 

* <a name="ref1">Lindquist, M. A., Waugh, C., & Wager, T. D. (2007). Modeling state-related fMRI activity using change-point theory. <i>NeuroImage</i>, <i>34</i>(3), 1125-1141. <https://doi.org/10.1016/j.neuroimage.2007.01.004> 
* <a name="ref3">Song, X. W., Dong, Z. Y., Long, X. Y., Li, S. F., Zuo, X. N., Zhu, C. Z., ... & Zang, Y. F. (2011). REST: A toolkit for resting-state functional magnetic resonance imaging data processing. <i>PloS one</i>, <i>6</i>(9), e25031. <https://doi.org/10.1371/journal.pone.0025031> 
* <a name="ref2">Tagliazucchi, E., Balenzuela, P., Fraiman, D., & Chialvo, D. R. (2012). Criticality in large-scale brain fMRI dynamics unveiled by a novel point process analysis. <i>Frontiers in Physiology</i>, <i>3</i>, 15. <https://doi.org/10.3389/fphys.2012.00015> 

<a href="#top">🔝</a>
