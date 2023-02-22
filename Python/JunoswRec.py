from math import sqrt
from SaveFiles import read4npy, write2npy

class Deconvolution():
    name = 'Deconvolution'
    tot_LPMT = 17612
    n_pmt_types = 2

    def __init__(self, L):
        import os
        import numpy as np
        from ROOT import TFile, TVirtualFFT
        import ctypes
        from tqdm import trange

        ptr_L = ctypes.byref(ctypes.c_int(L))
        print('Generating FFT plan: forward')
        self.fft_forward = TVirtualFFT.FFT(1, ptr_L, "R2C EX K")
        print('Done\nGenerating FFT plan: back')
        self.fft_back    = TVirtualFFT.FFT(1, ptr_L, "C2R EX K")
        print('Done')
        self.wf   = np.zeros(L)
        self.wfre = np.zeros(L)
        self.wfim = np.zeros(L)

        envs = os.environ
        if 'JUNOTOP' in envs.keys():
            junotop=envs['JUNOTOP']
        else:
            print(f'Warn: No JUNO env found, now to use J22.2.0-rc2')
            junotop='/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J22.2.0-rc2'
        junover=junotop.split('/')[-1]

        tag = 'filter'
        data_tag = f'{junover}_{self.name}_{tag}'
        self.Filter = read4npy(data_tag)
        if self.Filter.size == 0:
            print(f'Reading {data_tag}...')
            filterF = TFile(f'{junotop}/data/Reconstruction/Deconvolution/share/filter3_m.root')
            tmp_filter = []
            for i_type in range(self.n_pmt_types):
                fname = 'fh0' if i_type else 'fn0'
                filter1d = filterF.Get(fname)
                tmp_filter.append([])
                for i in trange(filter1d.GetNbinsX()):
                    tmp_filter[i_type].append(filter1d.GetBinContent(i+1))
            filterF.Close()
            self.Filter = np.array(tmp_filter)
            write2npy(self.Filter, data_tag)
        self.n_filter = self.Filter.shape[-1]
        print(f'Read filter: {self.Filter.shape}')

        tag = 'SPEfreq'
        data_tag = f'{junover}_{self.name}_{tag}'
        self.Freq = read4npy(data_tag)
        if self.Freq.size == 0:
            print(f'Reading {data_tag}...')
            freqF = TFile(f'{junotop}/data/Reconstruction/Deconvolution/share/SPE_v20.root')
            freq_names = ('RE', 'IM')
            tmp_freq = []
            for j in range(len(freq_names)):
                tmp_freq.append([])
                for i in trange(self.tot_LPMT):
                    tmp_freq[j].append([])
                    tmp1d = freqF.Get(f"PMTID_{i}_SPE{freq_names[j]}")
                    for n in range(tmp1d.GetNbinsX()):
                        tmp_freq[j][i].append(tmp1d.GetBinContent(n+1))
            freqF.Close()
            self.Freq = np.array(tmp_freq)
            write2npy(self.Freq, data_tag)
        self.n_freq = self.Freq.shape[-1]
        print(f'Read filter: {self.Freq.shape}')


def subBSL_NTW(raw_wf, N=3, L_bsl=50):
    if isinstance(raw_wf, (list, tuple)):
        bsls     = tuple( raw_wf[i*L_bsl:(i+1)*L_bsl] for i in range(N) )
        bsl_sums = tuple(                sum(bsls[i]) for i in range(N) )

        index = bsl_sums.index(min(bsl_sums))
        final_bsl = bsls[index]
        bsl = bsl_sums[index]/L_bsl
        bsl_sigma = sqrt(sum(tuple( (final_bsl[i]-bsl)*(final_bsl[i]-bsl) for i in range(L_bsl) )) / (L_bsl-1))
        return raw_wf-bsl, bsl_sigma
    else:
        print('Warning subBSL_3TW: unknown input type!')

def getNPE_AdcSum(wf, SPEadcSum):
    return { 'TTQ':sum(wf)/SPEadcSum }

def getNPE_OverZeroSum(wf, SPEOverZeroSum):
    return { 'TTQ':sum(filter(lambda x:x>0, wf))/SPEOverZeroSum }

def getNPE_AB(wf, amp_threshold, itg_threshold, SPEadcSum):
    N, passT, startT, endT, i = len(wf), 0, 0, 0, 1
    Q, T = [], []
    while i<N:
        if wf[i]>=amp_threshold and wf[i-1]<=amp_threshold:
            passT = i
            for j in range(passT,0,-1):
                if wf[j]<=0 or j==1:
                    startT = j
                    break
            for j in range(passT, N):
                if wf[j]<=0 or j==N-1:
                    endT = j
                    break
            tmpAdcSum = sum(wf[startT:endT])
            if tmpAdcSum>itg_threshold:
                Q.append(tmpAdcSum/SPEadcSum)
                T.append(startT)
            i = endT+1
        else:
            i += 1
    if not T:
        Q.append(getNPE_AdcSum(wf, SPEadcSum)['TTQ'])
        T.append(0)
    return {'TTQ':sum(Q), 'FHT':T[0], 'Q':Q, 'T':T}

if __name__ == "__main__":
    a = Deconvolution(1000)                
    print('test finished!')
