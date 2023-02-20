from math import sqrt

class Deconv():
    name = 'Deconvolution'
    tot_LPMT = 17612

    def __init__(self, L):
        import os
        from tqdm import trange
        from array import array
        from ROOT import TFile
        from ROOT.TVirtualFFT import FFT
        self.L  = array('d', [L])
        self.fft_forward = FFT(1, self.L, "R2C EX K")
        self.fft_back    = FFT(1, self.L, "C2R EX K")
        self.wf  = array('d',[0 for _ in range(L)])
        self.wfre= array('d',[0 for _ in range(L)])
        self.wfim= array('d',[0 for _ in range(L)])

        envs = os.environ
        if 'JUNOTOP' in envs.keys():
            junotop=envs['JUNOTOP']
        else:
            print(f'Warn: No JUNO env found, now to use J22.2.0-rc1')
            junotop='/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J22.2.0-rc1'
        junover=junotop.split('/')[-1]

        tag = 'filter'
        data_tag = f'{junover}_{self.name}_{tag}'
        data = read4txt(data_tag)
        if not data:
            print(f'Reading {data_tag}...')
            filterF = TFile(f'{junotop}/data/Reconstruction/Deconvolution/share/filter3_m.root')
            tmp_filter = {}
            for fname in ['fn0', 'fh0']:
                filter1d = filterF.Get(fname)
                pmttype = 'hmmt' if fname=='fh0' else 'nnvt'
                tmp_filter[pmttype] = []
                for i in range(filter1d.GetNbinsX()):
                    tmp_filter[pmttype].append(filter1d.GetBinContent(i+1))
            filterF.Close()
            write2txt(tmp_filter, data_tag)
            self.Filter = (array('d', tmp_filter['nnvt']), array('d', tmp_filter['hmmt']))
        else:
            self.Filter = (array('d',       data['nnvt']), array('d',       data['hmmt']))
        self.N_filter = self.Filter[0].itemsize

        tag = 'SPEfreq'
        data_tag = f'{junover}_{self.name}_{tag}'
        data = read4txt(data_tag)
        if not data:
            print(f'Reading {data_tag}...')
            freqF = TFile(f'{junotop}/data/Reconstruction/Deconvolution/share/SPE_v20.root')
            freq_names = ('RE', 'IM')
            tmp_freq = {}
            for i in trange(self.tot_LPMT):
                for j in freq_names:
                    tmp1d = freqF.Get(f"PMTID_{i}_SPE{j}")
                    for n in range(tmp1d.GetNbinsX()):
                        tmp_freq[j].append(tmp1d.GetBinContent(n+1))
            freqF.Close()
            write2txt(tmp_freq, data_tag)
            self.SPERE = array('d', tmp_freq['RE'])
            self.SPEIM = array('d', tmp_freq['IM'])


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

                
    
