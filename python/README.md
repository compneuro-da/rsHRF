Resting state HRF estimation and deconvolution.
========
Please refer to https://github.com/guorongwu/rsHRF/ for MATLAB version

![BOLD HRF](https://github.com/guorongwu/rsHRF/raw/master/docs/BOLD_HRF.png)


The basic idea 
-------------

This toolbox is aimed to retrieve the onsets of pseudo-events triggering an hemodynamic response from resting state fMRI BOLD voxel-wise signal.
It is based on point process theory, and fits a model to retrieve the optimal lag between the events and the HRF onset, as well as the HRF shape, using either the canonical shape with two derivatives, or a (smoothed) Finite Impulse Response.

![BOLD HRF](http://users.ugent.be/~dmarinaz/example_hrf.png)

Once that the HRF has been retrieved for each voxel, it can be deconvolved from the time series (for example to improve lag-based connectivity estimates), or one can map the shape parameters everywhere in the brain (including white matter), and use the shape as a pathophysiological indicator.

![HRF map](http://users.ugent.be/~dmarinaz/FIR_Height_full_layout.png)

How to use the toolbox 
-------------

The input is voxelwise BOLD signal, already preprocessed according to your favorite recipe. Important thing are:

* bandpass filter in the 0.01-0.08 Hz interval (or something like that)
* z-score the voxel BOLD time series

To be on the safe side, these steps are performed again in the code.

The input can be images (3D or 4D), or directly matrices of [observation x voxels].

It is possible to use a temporal mask to exclude some time points (for example after scrubbing).

The demos allow you to run the analyses on several formats of input data.

Towards Python and BIDS-app 
-------------
Currently we are working to translate the toolbox to Python, and to build a [BIDS-app](https://bids-apps.neuroimaging.io/) out of it.

**The Command Line Interface rsHRF**

A command line interface has been built in python to facilitate easy use of the toolbox.

The CLI can be installed with ``pip install rsHRF`` and requires ``Python>=3.5``.

This standalone command takes care of all the necessary dependencies so that the tool is
usable straight out of the box.

Once installed, run ``rsHRF --help`` to see the required positional and optional arguments.

In essence, the whole usage of the application can be broken down to 7 major steps:

1. **The input:**

There are 2 ways one can input data to this application.

* A standalone ``.nii / .nii.gz`` file. This option can be accessed using the
``--input_file`` optional argument followed by the name of the file.

* A ``bids_dir`` positional argument which is the location of the BIDS formatted 
data-set directory.

Out of the above 2 options, one of them is always required and both cannot be supplied
at once.

2. **The mask / atlas files:**

There are 2 ways one can supply the corresponding mask / atlas files.

* A standalone ``.nii / .nii.gz`` file. This option can be accessed using the
``--atlas`` optional argument followed by the name of the file.

* The ``--brainmask`` argument which tells the application that the mask files
are present within the BIDS formatted data-set directory itself (which was supplied
as ``bids_dir`` as a positional argument).
 
Out of the above 2 options, one of them is always required and both cannot be supplied
at once. Also, ``--input_file`` and ``--brainmask`` together are an invalid combination.

The other 3 use-cases are explained below:

* ``--input_file`` and ``--atlas`` : The standalone ``atlas`` and the standalone
 ``input_file`` are passed to the application for the analysis and the outputs are
 determined accordingly.
 
* ``bids_dir`` and ``--atlas`` : The standalone ``atlas`` is used with ALL the input
files present in the ``bids_dir`` directory. Thus, the ``atlas`` serves as a common mask
for the whole BIDS formatted data-set.

* ``bids_dir`` and ``--brainmask`` : This should be used when for each input file present
in the BIDS formatted data-set, the corresponding mask file exists within the same data-set.
The application then pairs the input_files with their corresponding masks provided that
the 2 files share a common prefix.

3. **The output directory:** 

The output directory ``output_dir`` is the folder under which all the resulting
``.nii`` files will be stored. The application further makes folders for each of 
the participants / subjects if the supplied arguments are ``bids_dir`` and ``--brainmask``.

4. **The Analysis Level:**

There are 2 kinds of analysis that can be performed.

* ``participant`` : participant level analysis performs the analysis for each (or some) subject(s) 
present in the BIDS formatted data-set. This positional argument is mandatory when the type of 
input used is ``bids_dir`` as group level analysis is not supported yet. This should not be 
supplied when input is supplied with ``--input_file`` argument. Doing so shall result in an error.

* ``group`` : Coming Soon! - Use ``participant`` for now.

5. **The Analysis Method:** 

The analysis can be carried out using 3 estimation methods.

These are ``canon2dd``, ``sFIR`` and ``FIR``.

One of them needs to be supplied using the ``--estimation`` argument followed by
one of the above 3 choices.

6. **The input parameters:** 

There are many input parameters that can be supplied to customize the analysis.
Please see all of them under the ``Parameters`` heading under the documentation
by running ``rsHRF --help``.

7. **The participant labels:**

Specifying which subjects to perform the analysis on can be given as a space separated
list after specifying the ``--participant_label`` argument. Only the valid subjects
in the list (which are actually present in the BIDS directory) will be considered.
The ``sub-`` prefix should not be supplied.

**Few Examples:**

* Running the analysis with a single input file and a single mask file.

``rsHRF --input_file input.nii results --atlas mask.nii --estimation canon2dd``.

In the above example, the input file is ``input.nii``. The ``output_dir`` is ``results``
directory. The corresponding mask file supplied is ``mask.nii``.
The estimation method passed is ``canon2dd``. The analysis level is not to be supplied here.

* Running the analysis with a BIDS formatted data-set and a common mask file
to be used for all the input files present in the data-set.

Note: All input files in the BIDs directory need to be of the type ``*_preproc.nii`` or 
``*_preproc.nii.gz``. Also, they must be present in the ``func`` directory under their
respective subject / session folder.

``rsHRF input_dir results participant --atlas mask.nii --estimation sFIR``

In the above example, the ``output_dir`` is ``results`` directory. The 
corresponding mask file supplied is ``mask.nii``. The BIDS formatted data-set
lies in the ``input_dir`` directory. The analysis level is ``participant``.

* Running the analysis with a BIDS formatted data-set that also includes a
unique mask file for each of the input file present. 

Note: All input files in the BIDs directory need to be of the type ``*_preproc.nii`` or 
``*_preproc.nii.gz``. The corresponding mask files in the BIDs directory need to
be of the type ``*_brainmask.nii`` or ``*_brainmask.nii.gz``. Also, they must be 
present in the ``func`` directory under their respective subject / session folder.
Furthermore, two corresponding input and mask files need to have the same prefix.

For example, 2 corresponding input and mask files can be ``input_preproc.nii`` and
``input_brainmask.nii``. These 2 will then be paired up for analysis.

``rsHRF input_dir results participant --brainmask --estimation canon2dd --participant_label 001 002``.

In the above example, the ``output_dir`` is ``results`` directory. The BIDS formatted data-set
lies in the ``input_dir`` directory. The associated mask files also lie within the BIDS dataset.
The analysis level is ``participant``. The analysis is performed only for ``sub-001`` & ``sub-002``
out of all the available subjects in the BIDS dataset.

Collaborators 
-------------
* Guorong Wu
* Nigel Colenbier
* Sofie Van Den Bossche
* Daniele Marinazzo

* Madhur Tandon (Python - BIDS)
* Asier Erramuzpe (Python - BIDS)


**References**
--------

1. Guo-Rong Wu, Wei Liao, Sebastiano Stramaglia, Ju-Rong Ding, Huafu Chen, Daniele Marinazzo*. "A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data." Medical Image Analysis, 2013, 17:365-374. [PDF](https://github.com/guorongwu/rsHRF/raw/master/docs/2013_MIA.pdf)

2. Guo-Rong Wu, Daniele Marinazzo. "Sensitivity of the resting state hemodynamic response function estimation to autonomic nervous system fluctuations." Philosophical Transactions of the Royal Society A, 2016, 374: 20150190. [PDF](https://github.com/guorongwu/rsHRF/raw/master/docs/2016_PTA.pdf)

3. Guo-Rong Wu, Daniele Marinazzo. "Retrieving the Hemodynamic Response Function in resting state fMRI: methodology and applications." PeerJ PrePrints, 2015. [PDF](https://github.com/guorongwu/rsHRF/raw/master/docs/poster_OHBM2016_HRF.pdf)
