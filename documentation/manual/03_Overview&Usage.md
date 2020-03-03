<a name="top">
<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

📖  Overview and Usage
----
__IMPORTANT__: Please use Google Chrome to browse the _Overview and Usage_ Page! If you don't, you won't be able to expand the collapsibles. [<abbr title="Work In Progress"><i>WIP</i></abbr>]

# Overview 
<a href="https://github.com/compneuro-da/rsHRF/blob/update/dissemination/paper_2013_MIA.pdf"><img align="right" src="https://github.com/guorongwu/rsHRF/raw/master/docs/example_hrf.png" alt="BOLD_HRF" width="350"/></a> <!-- find other image to illustrate pseudo-point process + code to produce it -->
<p align="justify"><b>The basic idea.</b> According to the point process theory discrete <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> events (i.e. pseudo-events in the absence of an external stimulus) govern the brain dynamics at rest (e.g. <a href="#ref4">Tagliazucchi et al. 2012</a>). The <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox is aimed to retrieve the neuronal  onsets of these pseudo-events with no explicit stimulus and timing together with the hemodynamic response (<abbr title="resting-state hemodynamic response function">rsHRF</abbr>) it set off (<a href="#ref6">Wu et al., 2013</a>; <a href="#ref9">Wu & Marinazzo, 2015</a>; <a href="#ref8">Wu & Marinazzo, 2016</a>). To this end, the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox first identifies the pseudo-events, i.e. when the <i>standardized</i> resting-state <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> signal crosses a given threshold (1 SD; see Figure). Thereafter, a model is fitted to retrieve: <ol type="A">
<li>the <i>optimal lag</i> between the pseudo-events and the neuronal (<abbr title="resting-state hemodynamic response function">rsHRF</abbr>) onset; </li>
<li>the <i>shape of the estimated <abbr title="resting-state hemodynamic response function">rsHRF</abbr></i> which will depend on the by-the-toolbox predefined <abbr title="hemodynamic response function">HRF</abbr> basis functions. Users of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox (<a href=""><b>v2.2</b></a>) can choose one of eight options:<br><br>

| Script         | HRF basis functions         | 
| :------------- |:-------------| 
| <a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_temporal_basis.m"><code>rsHRF_estimation_temporal_basis.m</code> | <i>canontd</i>: a canonical HRF with its time derivative |
| | <i>canontdd</i>: a canonical HRF with its time and dispersion derivatives |
| | <i>Gamma Functions</i> with a variable number of basis functions (<i>k</i>), e.g. 3-5 |
| | <i>Fourier Set</i> with a default number of basis functions (<i>k</i>) equal to 3 |
| | <i>Fourier Set (Hanning)</i> with a default number of basis functions (<i>k</i>) equal to 3 |
| <a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_FIR.m"><code>rsHRF_estimation_FIR.m</code> | <i>FIR</i>: Finite Impulse Response |
| | <i>sFIR</i>: smoothed Finite Impulse Response |
| (<a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_impulseest.m"><code>rsHRF_estimation_impulseest.m</code>) | <i>non-parametric impulse response estimation</i>: not included in the rsHRF <abbr title="graphical user interface">GUI</abbr> |
</li></ol></p>

# Usage 
<p align="justify">Once that the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> has been retrieved for each voxel/vertex in the brain, you can:
 
<details open><summary><i>use the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> as a pathophysiological indicator</i> (by mapping the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> shape onto the brain surface and looking at the inter-subject variability);<!--[[4](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf)]--></summary>

<a href="https://figshare.com/articles/voxelwise_resting_state_HRF_shape_WM_and_GM_/7139702"><img align="right" src="https://github.com/guorongwu/rsHRF/raw/master/docs/FIR_Height_full_layout.png" alt="HRF_map" width="350"/></a>
<p align="justify">The shape of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> can be characterized by <b>three parameters</b>, namely response height (<abbr title="response height"><i>RH</i></abbr>), time to peak (<abbr title="time to peak"><i>TTP</i></abbr>), and Full Width at Half Maximum (<abbr title="Full Width at Half Maximum"><i>FWHM</i></abbr>). Each of these parameters can be mapped onto the brain surface (see Figure for an example: full brain map of the response height estimated using the Finite Impulse Response basis functions). Note that the full brain map covers the full brain surface, including <a href= "https://ndownloader.figshare.com/files/13141076">white matter</a> and <abbr title="cerebrospinal fluid">CSF</abbr>. To consult some example data, head over to <a href="https://neurovault.org/collections/3584/">NeuroVault</a>. The <b>number of pseudo-events</b> per voxel/vertex can also be mapped onto the brain surface. After mapping all parameters (i.e. <abbr title="response height"><i>RH</i></abbr>, <abbr title="time to peak"><i>TTP</i></abbr>, <abbr title="Full Width at Half Maximum"><i>FWHM</i></abbr>, number of pseudo-events) onto the brain surface for each voxel/vertex and subject, the subject-specific brain maps can be used to examine whether/how the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> is modulated by psycho-physiological factors (i.e. inter-subject hemodynamic variability; e.g. <a href="#ref2">post-traumatic stress disorder</a>, <a href="#ref10">autism spectrum disorder</a>, <a href="#ref7">chronic pain</a>, <a href="#ref5">consciousness</a>). With the <a href="https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dMVM.html">3dMVM function</a> embedded in AFNI, one can even run a multivariate analysis in which the three <abbr title="resting-state hemodynamic response function">rsHRF</abbr> parameters are modeled as multiple, simultaneous response variables (<a href="#ref1">Chen, Adleman, Saad, Leibenluft, & Cox, 2014</a>).</p> <!-- insert: nr. of events, examples for inter-subject hemodynamic variability, 3dMVM; With the <a href="https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dMVM.html">3dMVM function</a> embedded in AFNI, a multivariate analysis can be run in which these three parameters are modeled as multiple, simultaneous response variables (Chen, Adleman, Saad, Leibenluft, & Cox, 2014). -->
 
<!-- , including white matter,-->

</details>

<details open><summary><i>deconvolve the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> from the resting-state <abbr title="functional Magnetic Resonance Imaging">fMRI</abbr> <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> signal</i> (for example to improve lag-based connectivity estimates).</summary><br>
 
<!-- <p align="justify">The shape of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr>, and thus the time to peak (<abbr title="time to peak"><i>TTP</i></abbr>), differs for each voxel/vertex in the brain. As a result, two different voxels/vertices with pseudo-events at the same time might have disparate neuronal <abbr title="resting-state hemodynamic response function">rsHRF</abbr> onsets (for a schematic example, see Figure). Using non-deconvolved resting-state <abbr title="functional Magnetic Resonance Imaging">fMRI</abbr> <abbr title="Blood Oxygenation Level Dependent">BOLD</abbr> signal might therefore interpose confounders on temporal precedence, which can deteriorate lag-based connectivity estimates (e.g. functional connectivity). --> <!--As functional connectivity analysis is built on associating BOLD events on two different spatial locations but at the same time; eliminating these confounders is of essence [REF]. --> <!-- insert non-deconvolved signal in connectivity analysis toolbox -->
 
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

* <a name="ref1">Chen, G., Adleman, N. E., Saad, Z. S., Leibenluft, E., & Cox, R. W. (2014). Applications of multivariate modeling to neuroimaging group analysis: A comprehensive alternative to univariate general linear model. NeuroImage, 99, 571-588. https://doi.org/10.1016/j.neuroimage.2014.06.027
* <a name="ref2">Rangaprakash, D., Dretsch, M. N., Yan, W., Katz, J. S., Denney Jr, T. S., & Deshpande, G. (2017). Hemodynamic variability in soldiers with trauma: Implications for functional MRI connectivity studies. <i>NeuroImage: Clinical</i>, <i>16</i>, 409-417. <https://doi.org/10.1016/j.nicl.2017.07.016>
* <a name="ref3">Rangaprakash, D., Wu, G. R., Marinazzo, D., Hu, X., & Deshpande, G. (2018). Hemodynamic response function (HRF) variability confounds resting‐state fMRI functional connectivity. <i>Magnetic Resonance in Medicine</i>, <i>80</i>(4), 1697-1713. <https://doi.org/10.1002/mrm.27146>
* <a name="ref4">Tagliazucchi, E., Balenzuela, P., Fraiman, D., & Chialvo, D. R. (2012). Criticality in large-scale brain fMRI dynamics unveiled by a novel point process analysis. <i>Frontiers in Physiology</i>, <i>3</i>, 15. <https://doi.org/10.3389/fphys.2012.00015> 
* [:paperclip:](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf) <a name="ref5">Wu, G. R., Di Perri, C., Charland-Verville, V., Martial, C., Carrière, M., Vanhaudenhuyse, A., ... & Marinazzo, D. (2019). Modulation of the spontaneous hemodynamic response function across levels of consciousness. <i>NeuroImage</i>, <i>200</i>, 450-459. <https://doi.org/10.1016/j.neuroimage.2019.07.011> 
* [:paperclip:](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2013_MIA.pdf) <a name="ref6">Wu, G. R., Liao, W., Stramaglia, S., Ding, J. R., Chen, H., & Marinazzo, D. (2013). A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data. <i>Medical Image Analysis</i>, <i>17</i>(3), 365-374. <https://doi.org/10.1016/j.media.2013.01.003>  
* <a name="ref7">Wu, G. R., & Marinazzo, D. (2015). Point-process deconvolution of fMRI BOLD signal reveals effective connectivity alterations in chronic pain patients. <i>Brain Topography</i>, <i>28</i>(4), 541-547.
* [:paperclip:](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2016_PTA.pdf) <a name="ref8">Wu, G. R., & Marinazzo, D. (2016). Sensitivity of the resting-state haemodynamic response function estimation to autonomic nervous system fluctuations. <i>Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences</i>, <i>374</i>(2067), 20150190. <https://doi.org/10.1098/rsta.2015.0190> 
* [:paperclip:](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/poster_OHBM2016_HRF.pdf)[:paperclip:](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/rs_HRF_OHBM2018_Daniele.pdf) <a name="ref9">Wu, G. R., & Marinazzo, D. (2015). <i>Retrieving the Hemodynamic Response Function in resting state fMRI: Methodology and applications</i> (No. e1621). PeerJ PrePrints. 
 * <a name="ref10">Yan, W., Rangaprakash, D., & Deshpande, G. (2018). Aberrant hemodynamic responses in autism: Implications for resting state fMRI functional connectivity studies. <i>NeuroImage: Clinical</i>, <i>19</i>, 320-330. <https://doi.org/10.1016/j.nicl.2018.04.013>

<a href="#top">🔝</a>
<a href="https://github.com/compneuro-da/rsHRF/blob/update/documentation/manual/02_Installation.md">:arrow_backward:</a>
<a href="https://github.com/compneuro-da/rsHRF/blob/update/documentation/manual/04_Workflow.md">:arrow_forward:</a>
