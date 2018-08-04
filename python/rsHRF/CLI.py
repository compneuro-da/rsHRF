import os.path as op
from argparse import ArgumentParser
from bids.grabbids import BIDSLayout
import numpy as np
import warnings
from rsHRF import spm_dep, fourD_rsHRF

__version__ = "0.1.7"

def get_parser():
    parser = ArgumentParser(description='retrieves the onsets of pseudo-events triggering a '
                                        'haemodynamic response from resting state fMRI BOLD '
                                        'voxel-wise signal')

    group_input = parser.add_mutually_exclusive_group(required=True)

    group_input.add_argument('--input_file', action='store', type=op.abspath,
                             help='the absolute path to a single data file')

    group_input.add_argument('bids_dir', nargs='?', action='store', type=op.abspath,
                             help='the root folder of a BIDS valid dataset '
                                  '(sub-XXXXX folders should be found at the '
                                  'top level in this folder).')

    parser.add_argument('output_dir', action='store', type=op.abspath,
                        help='the output path for the outcomes of processing')

    parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                        'Multiple participant level analyses can be run independently '
                        '(in parallel) using the same output_dir.', choices=['participant'], nargs='?')

    parser.add_argument('--participant_label',
                        help='The label(s) of the participant(s) that should be analyzed. The label '
                             'corresponds to sub-<participant_label> from the BIDS spec '
                             '(so it does not include "sub-"). If this parameter is not '
                             'provided all subjects should be analyzed. Multiple '
                             'participants can be specified with a space separated list.',
                        nargs="+")

    group_mask = parser.add_mutually_exclusive_group(required=True)

    group_mask.add_argument('--atlas', action='store', type=op.abspath,
                            help='the absolute path to a single atlas file')

    group_mask.add_argument('--brainmask', action='store_true',
                            help='to enable the use of mask files present in the BIDS '
                                 'directory itself')

    group_para = parser.add_argument_group('Parameters')

    group_para.add_argument('--estimation', action='store',
                            choices=['canon2dd', 'sFIR', 'FIR'], required=True,
                            help='Choose the estimation procedure from '
                                 'canon2dd (canonical shape with 2 derivatives), '
                                 'sFIR (smoothed Finite Impulse Response) , '
                                 'FIR (Finite Impulse Response)')

    group_para.add_argument('--passband', action='store', type=float, nargs=2, metavar=('LOW_FREQ','HIGH_FREQ'),
                            default=[0.01, 0.08],
                            help='set intervals for bandpass filter, default is 0.01 - 0.08')

    group_para.add_argument('-T', action='store', type=float, nargs=1, default=3,
                            help='set T parameter')

    group_para.add_argument('-T0', action='store', type=float, nargs=1, default=3,
                            help='set T0 parameter')

    group_para.add_argument('-TD_DD', action='store', type=int, nargs=1, default=2,
                            help='set TD_DD parameter')

    group_para.add_argument('-AR_lag', action='store', type=float, nargs=1, default=1,
                            help='set AR_lag parameter')

    group_para.add_argument('--thr', action='store', type=float, nargs=1, default=1,
                            help='set thr parameter')

    group_para.add_argument('--len', action='store', type=int, nargs=1, default=24,
                            help='set len parameter')

    group_para.add_argument('--min_onset_search', action='store', type=int, nargs=1, default=4,
                            help='set min_onset_search parameter')

    group_para.add_argument('--max_onset_search', action='store', type=int, nargs=1, default=8,
                            help='set max_onset_search parameter')

    return parser


def run_rsHRF():
    parser = get_parser()

    args = parser.parse_args()

    arg_groups = {}

    for group in parser._action_groups:
        group_dict = {a.dest: getattr(args, a.dest, None) for a in group._group_actions }
        arg_groups[group.title] = group_dict

    para = arg_groups['Parameters']

    if args.input_file is not None and args.analysis_level:
        parser.error('analysis_level cannot be used with --input_file, do not supply it')

    if args.input_file is not None and args.participant_label:
        parser.error('participant_labels are not to be used with --input_file, do not supply it')

    if args.input_file is not None and args.brainmask:
        parser.error('--brainmask cannot be used with --input_file, use --atlas instead')

    if args.bids_dir is not None and not args.analysis_level:
        parser.error('analysis_level needs to be supplied with bids_dir, choices=[participant]')

    if args.input_file is not None and (not args.input_file.endswith(('.nii', '.nii.gz'))):
        parser.error('--input_file should end with .nii or .nii.gz')

    if args.atlas is not None and (not args.atlas.endswith(('.nii', '.nii.gz'))):
        parser.error('--atlas should end with .nii or .nii.gz')

    if args.input_file is not None and args.atlas is not None:
        # carry analysis with input_file and atlas
        TR = spm_dep.spm.spm_vol(args.input_file).header.get_zooms()[-1]
        para['TR'] = TR
        para['dt'] = para['TR'] / para['T']
        para['lag'] = np.arange(np.fix(para['min_onset_search'] / para['dt']),
                                np.fix(para['max_onset_search'] / para['dt']) + 1,
                                dtype='int')
        fourD_rsHRF.demo_4d_rsHRF(args.input_file, args.atlas, args.output_dir, para, mode='input w/ atlas')

    if args.bids_dir is not None and args.atlas is not None:
        # carry analysis with bids_dir and 1 atlas
        layout = BIDSLayout(args.bids_dir)

        if args.participant_label:
            input_subjects = args.participant_label
            subjects_to_analyze = layout.get_subjects(subject=input_subjects)
        else:
            subjects_to_analyze = layout.get_subjects()

        if not subjects_to_analyze:
            parser.error('Could not find participants. Please make sure the BIDS data '
                         'structure is present and correct. Datasets can be validated online '
                         'using the BIDS Validator (http://incf.github.io/bids-validator/).')

        all_inputs = layout.get(modality='func', subject=subjects_to_analyze, task='rest', type='preproc', extensions=['nii', 'nii.gz'])
        if not all_inputs != []:
            parser.error('There are no files of type *preproc.nii / *preproc.nii.gz '
                         'Please make sure to have at least one file of the above type '
                         'in the BIDS specification')
        else:
            for file_count in range(len(all_inputs)):
                try:
                    TR = layout.get_metadata(all_inputs[file_count].filename)['RepetitionTime']
                except KeyError as e:
                    TR = spm_dep.spm.spm_vol(all_inputs[file_count].filename).header.get_zooms()[-1]
                para['TR'] = TR
                para['dt'] = para['TR'] / para['T']
                para['lag'] = np.arange(np.fix(para['min_onset_search'] / para['dt']),
                                        np.fix(para['max_onset_search'] / para['dt']) + 1,
                                        dtype='int')
                fourD_rsHRF.demo_4d_rsHRF(all_inputs[file_count], args.atlas, args.output_dir, para, mode='bids w/ atlas')

    if args.bids_dir is not None and args.brainmask:
        # carry analysis with bids_dir and brainmask
        layout = BIDSLayout(args.bids_dir)

        if args.participant_label:
            input_subjects = args.participant_label
            subjects_to_analyze = layout.get_subjects(subject=input_subjects)
        else:
            subjects_to_analyze = layout.get_subjects()

        if not subjects_to_analyze:
            parser.error('Could not find participants. Please make sure the BIDS data '
                         'structure is present and correct. Datasets can be validated online '
                         'using the BIDS Validator (http://incf.github.io/bids-validator/).')

        all_inputs = layout.get(modality='func', subject=subjects_to_analyze, task='rest', type='preproc', extensions=['nii', 'nii.gz'])
        all_masks = layout.get(modality='func', subject=subjects_to_analyze, task='rest', type='brainmask', extensions=['nii', 'nii.gz'])

        if not all_inputs != []:
            parser.error('There are no files of type *preproc.nii / *preproc.nii.gz '
                         'Please make sure to have at least one file of the above type '
                         'in the BIDS specification')
        if not all_masks != []:
            parser.error('There are no files of type *brainmask.nii / *brainmask.nii.gz '
                         'Please make sure to have at least one file of the above type '
                         'in the BIDS specification')
        if len(all_inputs) != len(all_masks):
            parser.error('The number of *preproc.nii / .nii.gz and the number of '
                         '*brainmask.nii / .nii.gz are different. Please make sure that '
                         'there is one mask for each input_file present')

        all_inputs.sort()
        all_masks.sort()

        all_prefix_match = False
        prefix_match_count = 0
        for i in range(len(all_inputs)):
            input_prefix = all_inputs[i].filename.split('/')[-1].split('_preproc')[0]
            mask_prefix = all_masks[i].filename.split('/')[-1].split('_brainmask')[0]
            if input_prefix == mask_prefix:
                prefix_match_count += 1
            else:
                all_prefix_match = False
                break
        if prefix_match_count == len(all_inputs):
            all_prefix_match = True

        if not all_prefix_match:
            parser.error('The mask and input files should have the same prefix for correspondence. '
                         'Please consider renaming your files')
        else:
            for file_count in range(len(all_inputs)):
                try:
                    TR = layout.get_metadata(all_inputs[file_count].filename)['RepetitionTime']
                except KeyError as e:
                    TR = spm_dep.spm.spm_vol(all_inputs[file_count].filename).header.get_zooms()[-1]
                para['TR'] = TR
                para['dt'] = para['TR'] / para['T']
                para['lag'] = np.arange(np.fix(para['min_onset_search'] / para['dt']),
                                        np.fix(para['max_onset_search'] / para['dt']) + 1,
                                        dtype='int')
                fourD_rsHRF.demo_4d_rsHRF(all_inputs[file_count], all_masks[file_count], args.output_dir, para, mode='bids')


def main():
    warnings.filterwarnings("ignore")
    run_rsHRF()


if __name__ == '__main__':
    raise RuntimeError("CLI.py should not be run directly;\n"
                       "Please `pip install` rsHRF and use the `rsHRF` command")
