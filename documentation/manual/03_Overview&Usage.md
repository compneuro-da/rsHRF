<a name="top">
<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

📖  Overview and Usage
----
__IMPORTANT__: Please use Google Chrome to browse the _Overview and Usage_ Page! If you don't, you won't be able to expand the collapsibles. [<abbr title="Work In Progress"><i>WIP</i></abbr>]

# Overview 
 <img align="right" src="https://github.com/guorongwu/rsHRF/raw/master/docs/example_hrf.png" alt="BOLD_HRF" width="350"/> <!-- find other image to illustrate pseudo-point process + code to produce it -->
<p align="justify"><b>The basic idea.</b> According to the point process theory discrete <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> events (i.e. pseudo-events in the absence of an external stimulus) govern the brain dynamics at rest (e.g. Tagliazucchi et al. 2012). The <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox is aimed to retrieve the neuronal  onsets of these pseudo-events with no explicit stimulus and timing, along with the hemodynamic response (<abbr title="resting-state hemodynamic response function">rsHRF</abbr>) it set off. To this end, the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox first identifies the pseudo-events, i.e. when the <i>standardized</i> resting-state <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> signal crosses a given threshold (1 SD; see Figure). Thereafter, a model is fitted to retrieve: <ol>
<li>the <i>optimal lag</i> between the pseudo-events and the neuronal (<abbr title="resting-state hemodynamic response function">rsHRF</abbr>) onset; </li>
 <li>the <i>shape of the estimated <abbr title="resting-state hemodynamic response function">rsHRF</abbr></i> which will depend on the by-the-toolbox predefined <abbr title="hemodynamic response function">HRF</abbr> basis functions. Users of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox (<a href=""><b>v2.2</b></a>) can choose one of five options:</li><ul>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/master/wgr_rshrf_estimation_canonhrf2dd_par2.m"><code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code></a>: a canonical HRF with its delay and dispersion derivatives – <i>canon2dd</i>;</li>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/master/wgr_rsHRF_FIR.m"><code>wgr_rsHRF_FIR.m</code></a>: (smoothed) Finite Impulse Response basis functions – <i>(s)FIR</i>;</li>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_FIR.m"><code>rsHRF_estimation_FIR.m</code></a>: an update of the (smoothed) Finite Impulse Response basis functions – <i>updated (s)FIR</i>;</li>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_temporal_basis.m"><code>rsHRF_estimation_temporal_basis.m</code></a>: a Gamma/Fourier basis function;</li>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_impulseest.m"><code>rsHRF_estimation_impulseest.m</code></a>: non-parametric impulse response estimation (which is not included in the rsHRF <abbr title="graphical user interface">GUI</abbr>).</li>
</ul></ol></p>

# Usage 
<p align="justify">Once that the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> has been retrieved for each voxel/vertex in the brain, you can:
 
<details><summary><i>use the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> as a pathophysiological indicator</i> (by mapping the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> shape onto the brain surface and looking at the inter-subject variability);<!--[[4](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf)]--></summary><br>

<img align="right" src="https://github.com/guorongwu/rsHRF/raw/master/docs/FIR_Height_full_layout.png" alt="HRF_map" width="350"/>
<p align="justify">The shape of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> can be characterized by <b>three parameters</b>, namely response height (<abbr title="response height"><i>RH</i></abbr>), time to peak (<abbr title="time to peak"><i>TTP</i></abbr>), and Full Width at Half Maximum (<abbr title="Full Width at Half Maximum"><i>FWHM</i></abbr>). Each of these parameters can be mapped onto the brain surface (see Figure for an example: full brain map of the response height estimated using the Finite Impulse Response basis functions). Note that the full brain map covers the full brain surface, including white matter and <abbr title="cerebrospinal fluid">CSF</abbr>. For more information, head over to the WORKFLOW page!</p> <!-- With the <a href="https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dMVM.html">3dMVM function</a> embedded in AFNI, a multivariate analysis can be run in which these three parameters are modeled as multiple, simultaneous response variables (Chen, Adleman, Saad, Leibenluft, & Cox, 2014). -->
 
<!-- , including white matter,-->

</details>

<details><summary><i>deconvolve the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> from the resting-state <abbr title="functional Magnetic Resonance Imaging">fMRI</abbr> <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> signal</i> (for example to improve lag-based connectivity estimates).</summary><br>
 
<!-- <p align="justify">The shape of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr>, and thus the time to peak (<abbr title="time to peak"><i>TTP</i></abbr>), differs for each voxel/vertex in the brain. As a result, two different voxels/vertices with pseudo-events at the same time might have disparate neuronal <abbr title="resting-state hemodynamic response function">rsHRF</abbr> onsets (for a schematic example, see Figure). Using non-deconvolved resting-state <abbr title="functional Magnetic Resonance Imaging">fMRI</abbr> <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> signal might therefore interpose confounders on temporal precedence, which can deteriorate lag-based connectivity estimates (e.g. functional connectivity). --> <!--As functional connectivity analysis is built on associating BOLD events on two different spatial locations but at the same time; eliminating these confounders is of essence [REF]. -->
 
</details>

<!--
<b>rsHRF deconvolution.</b> Improve lag-based connectivity estimates.
<b>rsHRF retrieval.</b> The rsHRF shape as a pathophysiological indicator. -->
<!-- <details><summary><b></b></summary>
</details> -->
<!--rsHRF different for each voxel in the brain [REF].-->

<!--
<p align="justify"><b>use shape as a biomarker</b>to retrieve the shape of the HRF concretized by three HRF parameters (see Figure 4); thee then use them in multivariate model; or look at some of them sepatetely
, or one can map the shape parameters everywhere in the brain (including white matter), and use the shape as a pathophysiological indicator [[4](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf)]. </p>
-->
<!-- collapsible -->
<!--
<p align="justify"><b>eliminate time-lag confounds.</b> Once that the HRF has been retrieved for each voxel/vertex, it can be deconvolved from the time series (for example to improve lag-based connectivity estimates)
  to deconvolve the resting-state signal and ilimante timing confounds (as the HRF shape is different for each voxel in the brain, the time to peak is different as well; therefore even as two voxels would have a pseudo-event at the same time, the timing of the corresponding neuronal events might not coincide (see FIgure 3 as an example). As functional connectivity analysis is built on associating BOLD events on two different spatial locations but at the same time; elliminating such time confounds is of essence [REF].</p>
-->

<!--
  the canonical shape with two derivatives, Gamma functions, Fourier set (Hanning), or a (smooth) Finite Impulse Response. -->
  
  <!--
  [[4](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf)] -->

# References

<!-- Chen, G., Adleman, N. E., Saad, Z. S., Leibenluft, E., & Cox, R. W. (2014). Applications of multivariate modeling to neuroimaging group analysis: A comprehensive alternative to univariate general linear model. NeuroImage, 99, 571-588. https://doi.org/10.1016/j.neuroimage.2014.06.027  -->

* <a name="ref1">Wu, G. R., Liao, W., Stramaglia, S., Ding, J. R., Chen, H., & Marinazzo, D. (2013). A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data. <i>Medical Image Analysis</i>, <i>17</i>(3), 365-374. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2013_MIA.pdf) 
* <a name="ref2">Wu, G. R., & Marinazzo, D. (2016). Sensitivity of the resting-state haemodynamic response function estimation to autonomic nervous system fluctuations. <i>Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences</i>, <i>374</i>(2067), 20150190. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2016_PTA.pdf)
* <a name="ref3">Wu, G. R., & Marinazzo, D. (2015). <i>Retrieving the Hemodynamic Response Function in resting state fMRI: Methodology and applications</i> (No. e1621). PeerJ PrePrints. [Poster2016](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/poster_OHBM2016_HRF.pdf), [Poster2018](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/rs_HRF_OHBM2018_Daniele.pdf)
* <a name="ref4">Wu, G. R., Di Perri, C., Charland-Verville, V., Martial, C., Carrière, M., Vanhaudenhuyse, A., ... & Marinazzo, D. (2019). Modulation of the spontaneous hemodynamic response function across levels of consciousness. <i>NeuroImage</i>, <i>200</i>, 450-459. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf) 

<a href="#top">🔝</a>
