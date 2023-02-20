from SaveFiles import *
from ROOT import TFile
import os
from tqdm import trange

class LPMTINFO():
    maxpmtid = 17612
    type_names = ('nnvt','hmmt')
    calib_names = []

    pmt_types = []
    calib_paras = {}

    # def __init__(self, enableCalib=True, enableCalib_coti=False, enableSPEIntegral=True, enableSPEOverZeroSum=True, enableUpdate=False):
    def __init__(self, enableCalib=True, enableCalib_coti=False, enableSPEIntegral=True, enableSPEOverZeroSum=True):
        envs = os.environ
        if 'JUNOTOP' in envs.keys():
            junotop=envs['JUNOTOP']
        else:
            print(f'Warn: No JUNO env found, now to use J22.1.0-rc2')
            junotop='/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J22.1.0-rc2'
        junover=junotop.split('/')[-1]

        with open(f'{junotop}/data/Detector/Geometry/PMTType_CD_LPMT.csv', 'r') as infile:
            for line in infile:
                self.pmt_types.append(line.split()[1][-1]=='u')
        self.maxpmtid = len(self.pmt_types)
        print(f'total LPMT: {self.maxpmtid}, number of hmmt: {sum(self.pmt_types)}')

        if enableCalib:
            tags = ['relativeDE', 'SPEratio', 'timeOffset', 'darkNoiseRate', 'meanRatio']
            for tag in tags:
                self.calib_paras[tag] = []
            self.calib_names += tags
            print(f'Reading calibration parameters of Deconvolution...')
            with open(f'{junotop}/data/Calibration/PMTCalibSvc/data/PmtPrtData_deconv.txt') as infile:
                for line in infile:
                    calib_info = line.split()
                    for i in range(len(tags)):
                        self.calib_paras[tags[i]].append(float(calib_info[i+1]))
        
        if enableCalib_coti:
            tags = ['relativeDE_coti', 'SPEratio_coti', 'timeOffset_coti', 'darkNoiseRate_coti', 'meanRatio_coti']
            for tag in tags:
                self.calib_paras[tag] = []
            self.calib_names += tags
            print(f'Reading calibration parameters of COTI...')
            with open(f'{junotop}/data/Calibration/PMTCalibSvc/data/PmtPrtData_inte.txt') as infile:
                for line in infile:
                    calib_info = line.split()
                    for i in range(len(tags)):
                        self.calib_paras[tags[i]].append(float(calib_info[i+1]))

        if enableSPEIntegral:
            tag = 'SPEIntegral'
            self.calib_names.append(tag)
            self.calib_paras[tag] = []

            data_tag = f'{junover}_{tag}'
            data = read4txt(data_tag)
            if not data:
                print(f'Reading {data_tag}...')
                SPEf = TFile.Open(f'{junotop}/data/Reconstruction/Deconvolution/share/SPE_v20.root')
                for pmtid in trange(self.maxpmtid):
                    self.calib_paras[tag].append( SPEf.Get(f'PMTID_{pmtid}_MeanSpec').Integral() )
                SPEf.Close()
                write2txt([tag, self.calib_paras[tag]], data_tag)
            else:
                self.calib_paras[tag] = data[tag]
                print(f'Info: Done to read {data_tag} from exist file')
        
        if enableSPEOverZeroSum:
            tag = 'SPEOverZeroSum'
            self.calib_names.append(tag)
            self.calib_paras[tag] = []

            data_tag = f'{junover}_{tag}'
            data = read4txt(data_tag)
            if not data:
                print(f'Reading {data_tag}...')
                SPEf = TFile.Open(f'{junotop}/data/Reconstruction/Deconvolution/share/SPE_v20.root')
                for pmtid in trange(self.maxpmtid):
                    adcSum = 0.
                    spe1d = SPEf.Get(f'PMTID_{pmtid}_MeanSpec')
                    for i in range(spe1d.GetNbinsX()):
                        if spe1d.GetBinContent(i+1) > 0:
                            adcSum += spe1d.GetBinContent(i+1)
                    self.calib_paras[tag].append( adcSum )
                SPEf.Close()
                write2txt([tag, self.calib_paras[tag]], data_tag)
            else:
                self.calib_paras[tag] = data[tag]
                print(f'Info: Done to read {data_tag} from exist file')
        
        self.calib_names = trans2tuple(self.calib_names)
        self.pmt_types     = trans2tuple(self.pmt_types)
        self.calib_paras = trans2tuple(self.calib_paras)

    def getNPMT(self):
        return self.maxpmtid
    
    def isHmmt(self, pmtid):
        return self.pmt_types[pmtid]

    def getType(self, pmtid):
        return self.type_names[self.isHmmt(pmtid)]

    def getCalib(self, pmtid, type):
        if type in self.calib_names:
            return self.calib_paras[type][pmtid]
        else:
            print(f'{type} not found!')
            return 0

def CdID2pmtId(cdid):
    return (cdid-(0x10<<24))>>8

class Hep_Jobs():
    hep_jobs = {}
    def __init__(self, j_user='') -> None:
        hepjobs = os.popen(f'hep_q -u {j_user}').read().split('\n')
        jobLines= tuple(filter(lambda x: x[-3:]=='.sh', hepjobs))
        for jl in jobLines:
            print(jl)
            jinfo = jl.strip('.sh').split()
            jtid, jst, jname = jinfo[0], jinfo[5], jinfo[8]
            tmpindex = max(jname.rfind('_'), jname.rfind('-'))+1
            if tmpindex<1:
                continue
            if jname[tmpindex:].isdigit():
                jid = int(jname[tmpindex:])
            else:
                continue
            jnb = jname[:tmpindex]
            if jnb not in self.hep_jobs.keys():
                self.hep_jobs[jnb] = {}
                for st in ('R', 'H', 'I'):
                    self.hep_jobs[jnb][st] = {}
                    for ids in ('tid', 'jid'):
                        self.hep_jobs[jnb][st][ids] = []
            self.hep_jobs[jnb][jst]['tid'].append(jtid)
            self.hep_jobs[jnb][jst]['jid'].append( jid)
        
        if not len(self.hep_jobs.keys()):
            print(f"No jobs got!")
    
    def getIDs(self, jnb, st):
        if jnb not in self.hep_jobs.keys():
            print(f'No {jnb} found!')
            return []
        return self.hep_jobs[jnb][st]['jid']
    
