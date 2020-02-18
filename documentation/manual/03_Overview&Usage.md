<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

📖  Overview and Usage
----
__IMPORTANT__: Please use Google Chrome to browse the _Overview and Usage_ Page! If you don't, you won't be able to expand the collapsibles. [<abbr title="Work In Progress"><i>WIP</i></abbr>]

# Overview 
<figure>
 <img align="right" src="https://github.com/guorongwu/rsHRF/raw/master/docs/example_hrf.png" alt="BOLD_HRF" width="350"/> <!-- find other image to illustrate pseudo-point process + code to produce it -->
 <figcaption align="right">Fig.1 - Trulli, Puglia, Italy.</figcaption>
</figure>
<!-- <img align="center" src="https://github.com/guorongwu/rsHRF/raw/master/docs/FIR_Height_full_layout.png" alt="HRF_map" width="800"/> -->
<!-- add text here -->
<p align="justify"><b>The basic idea.</b> According to the point process theory discrete BOLD events (i.e. pseudo-events in the absence of an external stimulus) govern the brain dynamics at rest (e.g. Tagliazucchi et al. 2012). The <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox is aimed to retrieve the neuronal onsets of these pseudo-events with no explicit stimulus and timing for the <abbr title="hemodynamic response function">HRF</abbr> onset. To this end, the toolbox first identifies the pseudo-events, i.e. when the <i>standardized</i> resting-state BOLD signal crosses a given threshold (1 SD; see Figure 1). Thereafter, a model is fitted to retrieve: <ol>
<li>the <i>optimal lag</i> between the pseudo-events and the neuronal (<abbr title="hemodynamic response function">HRF</abbr>) onset; </li>
 <li>the <i>shape of the estimated <abbr title="hemodynamic response function">HRF</abbr></i> which will depend on the by-the-toolbox predefined <abbr title="hemodynamic response function">HRF</abbr> basis function. Users of the <abbr title="resting-state hemodynamic response function">rsHRF</abbr> toolbox (<a href=""><b>v2.2</b></a>) can choose one of five options:</li><br> 

<ul>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/master/wgr_rshrf_estimation_canonhrf2dd_par2.m"><code>wgr_rshrf_estimation_canonhrf2dd_par2.m</code></a>: a canonical HRF with its delay and dispersion derivatives – <i>canon2dd</i>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/master/wgr_rsHRF_FIR.m"><code>wgr_rsHRF_FIR.m</code></a>: (smoothed) Finite Impulse Response basis functions – <i>(s)FIR</i></li>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_FIR.m"><code>rsHRF_estimation_FIR.m</code></a>: an update of the (smoothed) Finite Impulse Response basis functions – <i>updated (s)FIR</i></li>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_temporal_basis.m"><code>rsHRF_estimation_temporal_basis.m</code></a>: a Gamma/Fourier basis function</li>
<li><a href="https://github.com/compneuro-da/rsHRF/blob/update/code/rsHRF_estimation_impulseest.m"><code>rsHRF_estimation_impulseest.m</code></a>: non-parametric impulse response estimation (which is not included in the rsHRF <abbr title="graphical user interface">GUI</abbr>)</li>
</ul>
 
 <!-- OHBM uitleg-->
<p align="justify">In the Figure (for code see ...), NORMALIZED resting-state BOLD signal is illustrated; whenever the normalized signal exceeds a predefined treshold (i.e. SD = 1), that point in time is being marked as pseudo-event, indicated by SYMBOL. Therefafeter, the question remains how far the pseudo-event lags after the neuronal event.  A model WHICH ONE is used to retrieve the optimal lag between the events and the <abbr title="hemodynamic response function">HRF</abbr> onset, as well as the <abbr title="hemodynamic response function">HRF</abbr> shape. </p>

# Usage 
## Objectives
<p align="justify">HRF different for each voxel in the brain [REF] Then we can use it for two purposes:</p>

<p align="justify"><b>use shape as a biomarker</b>to retrieve the shape of the HRF concretized by three HRF parameters (see Figure 4); thee then use them in multivariate model; or look at some of them sepatetely
, or one can map the shape parameters everywhere in the brain (including white matter), and use the shape as a pathophysiological indicator [[4](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf)]. </p>

<!-- collapsible -->
<p align="justify"><b>eliminate time-lag confounds.</b> Once that the HRF has been retrieved for each voxel/vertex, it can be deconvolved from the time series (for example to improve lag-based connectivity estimates)
  to deconvolve the resting-state signal and ilimante timing confounds (as the HRF shape is different for each voxel in the brain, the time to peak is different as well; therefore even as two voxels would have a pseudo-event at the same time, the timing of the corresponding neuronal events might not coincide (see FIgure 3 as an example). As functional connectivity analysis is built on associating BOLD events on two different spatial locations but at the same time; elliminating such time confounds is of essence [REF].</p>


  the canonical shape with two derivatives, Gamma functions, Fourier set (Hanning), or a (smooth) Finite Impulse Response.

# References

1. Guo-Rong Wu, Wei Liao, Sebastiano Stramaglia, Ju-Rong Ding, Huafu Chen, Daniele Marinazzo. "A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data." Medical Image Analysis, 2013, 17:365-374. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2013_MIA.pdf)

2. Guo-Rong Wu, Daniele Marinazzo. "Sensitivity of the resting state hemodynamic response function estimation to autonomic nervous system fluctuations." Philosophical Transactions of the Royal Society A, 2016, 374: 20150190. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2016_PTA.pdf)

3. Guo-Rong Wu, Daniele Marinazzo. "Retrieving the Hemodynamic Response Function in resting state fMRI: methodology and applications." PeerJ PrePrints, 2015. [Poster2016](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/poster_OHBM2016_HRF.pdf),[Poster2018](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/rs_HRF_OHBM2018_Daniele.pdf)

4. Guo-Rong Wu, Carol Di Perri, Vanessa Charland-Verville, Charlotte Martial, Manon Carriere, Audrey Vanhaudenhuyse, Steven Laureys, Daniele Marinazzo. “Modulation of the spontaneous hemodynamic response function across levels of consciousness.” Neuroimage, 2019(200), 450–459. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf)
