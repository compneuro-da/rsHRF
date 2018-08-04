import numpy as np
import warnings

warnings.filterwarnings("ignore")

def nextpow2(n):
    return np.ceil(np.log2(np.abs(n))).astype('long')


def rest_nextpow2_one35(n):
    if type(n) is np.ndarray:
        if(np.amax(n.shape)) > 1:
            n = np.amax(n.shape)
    if n < 16:
        Result = 2 ** nextpow2(n)
        return Result
    limit = nextpow2(n)
    tbl = np.arange((2 ** (limit - 1)), ((2 ** limit) + 1))
    tbl = np.extract(tbl >= n, tbl)
    for x in range(0, np.amax(tbl.shape)):
        Result = tbl[x]
        f, p = np.frexp(Result)
        if f.size != 0 and f == 0.5:
            return Result
        if np.remainder(Result, 3 * 5) == 0:
            y = Result / (3 * 5)
            f, p = np.frexp(y)
            if f.size != 0 and f == 0.5:
                return Result
        if np.remainder(Result, 3) == 0:
            y = Result / 3
            f, p = np.frexp(y)
            if f.size != 0 and f == 0.5:
                return Result
        if np.remainder(Result, 5) == 0:
            y = Result / 5
            f, p = np.frexp(y)
            if f.size != 0 and f == 0.5:
                return Result
    Result = np.nan
    return Result


def rest_IdealFilter(Data, SamplePeriod, Band):
    sampleFreq = 1 / SamplePeriod
    sampleLength = Data.shape[0]
    paddedLength = rest_nextpow2_one35(sampleLength)
    LowCutoff_HighPass = Band[0]
    HighCutoff_LowPass = Band[1]
    if (LowCutoff_HighPass >= sampleFreq / 2):
        idxLowCutoff_HighPass = int(paddedLength / 2 + 1)
    else:
        idxLowCutoff_HighPass = int(
            np.ceil(LowCutoff_HighPass * paddedLength * SamplePeriod + 1)
        )

    if (HighCutoff_LowPass >= sampleFreq / 2) or (HighCutoff_LowPass == 0):
        idxHighCutoff_LowPass = int(paddedLength / 2 + 1)
    else:
        idxHighCutoff_LowPass = int(
            np.fix(HighCutoff_LowPass * paddedLength * SamplePeriod + 1)
        )

    FrequencyMask = np.zeros((paddedLength, 1))
    FrequencyMask[idxLowCutoff_HighPass - 1:idxHighCutoff_LowPass, 0] = 1
    FrequencyMask[paddedLength - idxLowCutoff_HighPass + 1:
                  paddedLength - idxHighCutoff_LowPass:-1, 0] = 1
    FrequencySetZero_Index = np.nonzero(FrequencyMask == 0)

    Data = Data - np.tile(np.mean(Data, axis=0), (Data.shape[0], 1))
    Data = np.concatenate(
        (Data, np.zeros((paddedLength - sampleLength, Data.shape[1]))), axis=0
    )
    Data = np.fft.fft(Data, axis=0)
    Data[FrequencySetZero_Index, :] = 0
    Data = np.fft.ifft(Data, axis=0)
    Data_Filtered = Data[0:sampleLength, :]
    return Data_Filtered
