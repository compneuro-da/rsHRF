<a href="https://github.com/compneuro-da/rsHRF/blob/update/README.md#table-of-contents">:leftwards_arrow_with_hook:</a> <br>

📖  Overview and Usage
-------------

The basic idea 
-------------

This toolbox is aimed to retrieve the onsets of pseudo-events triggering an hemodynamic response from resting state fMRI BOLD voxel-wise (vertex-wise) signal.
It is based on point process theory, and fits a model to retrieve the optimal lag between the events and the HRF onset, as well as the HRF shape, using either the canonical shape with two derivatives, Gamma functions, Fourier set (Hanning), or a (smooth) Finite Impulse Response.

![BOLD HRF](https://github.com/guorongwu/rsHRF/raw/master/docs/example_hrf.png)

Once that the HRF has been retrieved for each voxel/vertex, it can be deconvolved from the time series (for example to improve lag-based connectivity estimates), or one can map the shape parameters everywhere in the brain (including white matter), and use the shape as a pathophysiological indicator [[4](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf)].

![HRF map](https://github.com/guorongwu/rsHRF/raw/master/docs/FIR_Height_full_layout.png)

How to use the toolbox - Matlab standalone
-------------
N.B. it is still necessary to have SPM in the path, since some of the functions there are used.

The input is voxelwise/vertexwise BOLD signal, already preprocessed according to your favorite recipe. something like:

* nuisance variable regression 
* bandpass filter in the 0.01-0.08 Hz interval
* despike

(These denoising steps are also provided in the SPM plugin.)

The input can be images (3D or 4D), mesh (2D), or directly matrices of [observation x voxels(vertices)].

It is possible to use a temporal mask to exclude some time points (for example after scrubbing).

The demos allow you to run the analyses on several formats of input data.

How to use the toolbox - SPM plugin
-------------

The script spm_rsHRF.m is the main one, and it calls rsHRF.m. These two files are specific to the SPM plugin. 

See [rsHRF_toolbox.pptx](https://github.com/guorongwu/rsHRF/raw/master/rsHRF_toolbox.pptx) for more details (Installation/Usage/Outputs).
![rsHRF GUI](https://github.com/guorongwu/rsHRF_data/raw/master/rsHRF_GUI.png)

The connectivity analysis (functional connectivity: Pearson/Spearman correlation, Pairwise/Conditional/Partially Conditioned Granger causality) is only provided in the rsHRF SPM plugin. 

Flowcharts & videos!

(Version 1.0 (2018) can be downloaded [here](https://github.com/compneuro-da/rsHRF_data/raw/master/rsHRF_v1_2018.zip))

### add DEMOs (for each option)

#### Start resting state hrf (in set-up) linked there!!

which data to use? 

1. input format: voxels, rois, signals
2. 

**References**
--------

1. Guo-Rong Wu, Wei Liao, Sebastiano Stramaglia, Ju-Rong Ding, Huafu Chen, Daniele Marinazzo. "A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data." Medical Image Analysis, 2013, 17:365-374. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2013_MIA.pdf)

2. Guo-Rong Wu, Daniele Marinazzo. "Sensitivity of the resting state hemodynamic response function estimation to autonomic nervous system fluctuations." Philosophical Transactions of the Royal Society A, 2016, 374: 20150190. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2016_PTA.pdf)

3. Guo-Rong Wu, Daniele Marinazzo. "Retrieving the Hemodynamic Response Function in resting state fMRI: methodology and applications." PeerJ PrePrints, 2015. [Poster2016](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/poster_OHBM2016_HRF.pdf),[Poster2018](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/rs_HRF_OHBM2018_Daniele.pdf)

4. Guo-Rong Wu, Carol Di Perri, Vanessa Charland-Verville, Charlotte Martial, Manon Carriere, Audrey Vanhaudenhuyse, Steven Laureys, Daniele Marinazzo. “Modulation of the spontaneous hemodynamic response function across levels of consciousness.” Neuroimage, 2019(200), 450–459. [PDF](https://github.com/compneuro-da/rsHRF_data/raw/master/docs/2019_NI.pdf)
