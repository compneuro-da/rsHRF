import matplotlib
matplotlib.use('agg')
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import stats, signal
from scipy.sparse import lil_matrix
import scipy.io as sio
import warnings
from rsHRF import spm_dep, processing, canon, sFIR, parameters

warnings.filterwarnings("ignore")


def demo_4d_rsHRF(input_file, mask_file, output_dir, para, mode='bids'):
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    if mode == 'bids':
        name = input_file.filename.split('/')[-1].split('.')[0]
        v = spm_dep.spm.spm_vol(mask_file.filename)
    elif mode == 'bids w/ atlas':
        name = input_file.filename.split('/')[-1].split('.')[0]
        v = spm_dep.spm.spm_vol(mask_file)
    else:
        name = input_file.split('/')[-1].split('.')[0]
        v = spm_dep.spm.spm_vol(mask_file)
    brain = spm_dep.spm.spm_read_vols(v)

    voxel_ind = np.where(brain > 0)[0]

    temporal_mask = []

    if mode == 'bids' or mode == 'bids w/ atlas':
        v1 = spm_dep.spm.spm_vol(input_file.filename)
    else:
        v1 = spm_dep.spm.spm_vol(input_file)

    if v1.header.get_data_shape()[:-1] != v.header.get_data_shape():
        print('The dimension of your mask is different than '
              'the one of your fMRI data!')
        return
    else:
        data = v1.get_data()
        nobs = data.shape[3]
        data1 = np.reshape(data, (-1, nobs), order='F').T
        bold_sig = stats.zscore(data1[:, voxel_ind], ddof=1)
        bold_sig = np.nan_to_num(bold_sig)
        bold_sig = processing. \
            rest_filter. \
            rest_IdealFilter(bold_sig, para['TR'], para['passband'])
        data_deconv = np.zeros(bold_sig.shape)
        event_number = np.zeros((1, bold_sig.shape[1]))

        print('Retrieving HRF ...')

        if 'canon' in para['estimation']:
            beta_hrf, bf, event_bold = \
                canon.canon_hrf2dd.wgr_rshrf_estimation_canonhrf2dd_par2(
                    bold_sig, para, temporal_mask
                )
            hrfa = np.dot(bf, beta_hrf[np.arange(0, bf.shape[1]), :])
        elif 'FIR' in para['estimation']:
            para['T'] = 1
            hrfa, event_bold = sFIR. \
                smooth_fir. \
                wgr_rsHRF_FIR(bold_sig, para, temporal_mask)

        nvar = hrfa.shape[1]
        PARA = np.zeros((3, nvar))

        for voxel_id in range(nvar):
            hrf1 = hrfa[:, voxel_id]
            PARA[:, voxel_id] = \
                parameters.wgr_get_parameters(hrf1, para['TR'] / para['T'])

        print('Done')

        print('Deconvolving HRF ...')

        T = np.around(para['len'] / para['TR'])

        if para['T'] > 1:
            hrfa_TR = signal.resample_poly(hrfa, 1, para['T'])
        else:
            hrfa_TR = hrfa

        for voxel_id in range(nvar):
            hrf = hrfa_TR[:, voxel_id]
            H = np.fft.fft(
                np.append(hrf,
                          np.zeros((nobs - max(hrf.shape), 1))), axis=0)
            M = np.fft.fft(bold_sig[:, voxel_id])
            data_deconv[:, voxel_id] = \
                np.fft.ifft(H.conj() * M / (H * H.conj() + .1*np.mean((H * H.conj()))))
            event_number[:, voxel_id] = np.amax(event_bold[voxel_id].shape)

        if mode == 'bids' or mode == 'bids w/ atlas':
            try:
                sub_save_dir = os.path.join(
                    output_dir, 'sub-' + input_file.subject,
                    'session-' + input_file.session,
                    input_file.modality
                )
            except AttributeError as e:
                sub_save_dir = os.path.join(
                    output_dir, 'sub-' + input_file.subject,
                    input_file.modality
                )
        else:
            sub_save_dir = output_dir

        if not os.path.isdir(sub_save_dir):
            os.makedirs(sub_save_dir, exist_ok=True)

        sio.savemat(os.path.join(sub_save_dir, name + '_hrf.mat'),
                    {'para': para, 'hrfa': hrfa,
                     'event_bold': event_bold, 'PARA': PARA})
        HRF_para_str = ['Height.nii', 'Time2peak.nii', 'FWHM.nii']
        data = np.zeros(v.get_data().shape).flatten(order='F')

        for i in range(3):
            fname = os.path.join(sub_save_dir,
                                 name + '_' + HRF_para_str[i])
            data[voxel_ind] = PARA[i, :]
            data = data.reshape(v.get_data().shape, order='F')
            spm_dep.spm.spm_write_vol(v, data, fname)
            data = data.flatten(order='F')

        fname = os.path.join(sub_save_dir, name + '_event_number.nii')
        data[voxel_ind] = event_number
        data = data.reshape(v.get_data().shape, order='F')
        spm_dep.spm.spm_write_vol(v, data, fname)

        data = np.zeros(v1.get_data().shape)
        dat3 = np.zeros(v1.header.get_data_shape()[:-1]).flatten(order='F')
        for i in range(nobs):
            fname = os.path.join(sub_save_dir, name + '_deconv')
            dat3[voxel_ind] = data_deconv[i, :]
            dat3 = dat3.reshape(v1.header.get_data_shape()[:-1], order='F')
            data[:, :, :, i] = dat3
            dat3 = dat3.flatten(order='F')
        spm_dep.spm.spm_write_vol(v1, data, fname)

        event_plot = lil_matrix((1, nobs))
        event_plot[:, event_bold[0]] = 1
        event_plot = np.ravel(event_plot.toarray())

        plt.figure()
        plt.plot(para['TR'] * np.arange(1, np.amax(hrfa[:, 0].shape) + 1),
                 hrfa[:, 0], linewidth=1)
        plt.xlabel('time (s)')
        plt.savefig(os.path.join(sub_save_dir, name + '_plot_1.png'))

        plt.figure()
        plt.plot(para['TR'] * np.arange(1, nobs + 1),
                 np.nan_to_num(stats.zscore(bold_sig[:, 0], ddof=1)),
                 linewidth=1)
        plt.plot(para['TR'] * np.arange(1, nobs + 1),
                 np.nan_to_num(stats.zscore(data_deconv[:, 0], ddof=1)),
                 color='r', linewidth=1)
        markerline, stemlines, baseline = \
            plt.stem(para['TR'] * np.arange(1, nobs + 1), event_plot)
        plt.setp(baseline, 'color', 'k', 'markersize', 1)
        plt.setp(stemlines, 'color', 'k', 'markersize', 1)
        plt.setp(markerline, 'color', 'k', 'markersize', 3, 'marker', 'd')
        plt.legend(['BOLD', 'deconvolved', 'events'])
        plt.xlabel('time (s)')
        plt.savefig(os.path.join(sub_save_dir, name + '_plot_2.png'))
