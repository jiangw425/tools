from JunoswAna import *
from os.path import exists
from ROOT import TObject, TH1D, TH2D, TGraph

def trans2tuple(indata):
    if not indata or isinstance(indata, tuple):
        return indata
    elif isinstance(indata, dict):
        return {tag:trans2tuple(indata[tag]) for tag in indata}
    elif isinstance(indata, list):
        return tuple(indata)
    else:
        print(f'Warning: unknown input type!')

def write2root(outf, rtobj):
    if not rtobj:
        pass
    elif isinstance(rtobj, dict):
        for tmpobj in rtobj.values():
            write2root(outf, tmpobj)
    elif isinstance(rtobj, (list, tuple)):
        for tmpobj in rtobj:
            write2root(outf, tmpobj)
    elif isinstance(rtobj, TObject):
        outf.cd()
        rtobj.Write()
    else:
        print(f'Warning: unknown input type!')

def getCachePath(fname):
    def addtxt(s):
        # print(s if s[-4:]=='.txt' else s+'.txt')
        return s if s[-4:]=='.txt' else s+'.txt'

    if fname[0] == '/' or fname[:2]=='./':
        return addtxt(fname)
    else:
        cache_path = '/junofs/users/jiangw/include/cache'
        specified_names = ( 'elecTruth', )
        for sn in specified_names:
            if sn in fname:
                new_fname = fname.replace(sn, '')
                return addtxt(f'{cache_path}/{sn}/{new_fname}')
        return addtxt(f'{cache_path}/{fname}')

def fWriteNum(f, num, width=10):
    if isinstance(num, int):
        f.write(f'{num:>{width}} ')
    else:
        f.write(f'{num:>{width}.2f} ')

def write2txt(val, oname, width=10):
    outf = getCachePath(oname)

    if isinstance(val, (list, tuple)):
        if isinstance(val[0], (list, tuple)):
            # val : [tag, value]
            ntag = len(val)
            n = len(val[0][1])
            with open(outf,'w') as f:
                for i in range(ntag):
                    f.write(f'{val[i][0]} ')
                f.write('\n')
                for j in range(n):
                    for i in range(ntag):
                        fWriteNum(f, val[i][j], width)
                    f.write('\n')
        elif isinstance(val[0], str):
            n = len(val[1])
            with open(outf,'w') as f:
                f.write(f'{val[0]}\n')
                for j in range(n):
                    fWriteNum(f, val[1][j], width)
    elif isinstance(val, dict):
        tags = tuple(val.keys())
        ntag = len(tags)
        n = len(val[tags[0]])
        with open(outf,'w') as f:
            for i in range(ntag):
                f.write(f'{tags[i]} ')
            f.write('\n')
            for j in range(n):
                for i in range(ntag):
                    fWriteNum(f, val[tags[i]][j], width)
                f.write('\n')
    else:
        print(f'Warning: unknown input type!')

def read4txt(iname):
    val = {}
    tags = []
    inf = getCachePath(iname)
    # print('read from ',inf)
    isFirst = True
    if exists(inf):
        with open(inf,'r') as f:
            for line in f:
                infos = line.split()
                if isFirst:
                    for tag in infos:
                        tags.append(tag)
                        val[tag] = []
                    isFirst = False
                    ntag = len(tags)
                    continue
                for i in range(ntag):
                    val[tags[i]].append(float(infos[i]))
    else:
        print(f'Warn! {inf} not exist!')
    return trans2tuple(val)
    # val: { tag0:data0, tag1:data1, ... }

def generateTH1Ds(args, setData=False):
    th1ds = {}
    if isinstance(args, dict):
        # { name:[title, bins, low, up] }
        for name in args:
            cf = args[name]
            th1ds[name] = TH1D(name, cf[0], cf[1], cf[2], cf[3])
            if setData and len(cf[4])==cf[1]:
                for i in range(cf[1]):
                    th1ds[name].SetBinContent(i+1, cf[4][i])
    elif isinstance(args, (list, tuple)):
        # [ [name, title, bins, low, up] ]
        for cf in args:
            th1ds[cf[0]] = TH1D(cf[0], cf[1], cf[2], cf[3], cf[4])
            if setData and len(cf[5])==cf[2]:
                for i in range(cf[2]):
                    th1ds[cf[0]].SetBinContent(i+1, cf[5][i])
    else:
        print(f'Warning: unknown input type!')
    return th1ds

def generateTGraphs(args):
    tgraphs = {}
    if isinstance(args, dict):
        for name in args:
            cf = args[name]
            tgraphs[name] = TGraph()
            tgraphs[name].SetNameTitle(name,name)
            if not cf[0]:
                continue
            tgraphs[name].SetMinimum(min_more(cf[1]))
            tgraphs[name].SetMaximum(max_more(cf[1]))
            for i in range(len(cf[0])):
                tgraphs[name].SetPoint(i, cf[0][i], cf[1][i])
    elif isinstance(args, (list, tuple)):
        for cf in args:
            name = cf[0]
            tgraphs[name] = TGraph()
            tgraphs[name].SetNameTitle(name,name)
            if not cf[1]:
                continue
            tgraphs[name].SetMinimum(min_more(cf[2]))
            tgraphs[name].SetMaximum(max_more(cf[2]))
            for i in range(len(cf[1])):
                tgraphs[name].SetPoint(i, cf[1][i], cf[2][i])
    else:
        print(f'Warning: unknown input type!')
    return tgraphs

def generateTH2Ds(args):
    th2ds = {}
    if isinstance(args, dict):
        for name in args:
            cf = args[name]
            th2ds[name] = TH2D(name, cf[0], cf[1], cf[2], cf[3], cf[4], cf[5], cf[6])
    elif isinstance(args, (list, tuple)):
        for cf in args:
            th2ds[cf[0]] = TH2D(cf[0], cf[1], cf[2], cf[3], cf[4], cf[5], cf[6], cf[7])
    else:
        print(f'Warning: unknown input type!')
    return th2ds

def generateTProfiles(th2ds, axis='X'):
    tprofiles = {}
    if isinstance(th2ds, dict):
        for name in th2ds:
            tprofiles[name] = th2ds[name].ProfileX() if axis=='X' else th2ds[name].ProfileY()
    elif isinstance(th2ds, (list, tuple)):
        for th2d in th2ds:
            if isinstance(th2d, TH2D):
                tprofiles[th2d.GetName()] = th2d.ProfileX() if axis=='X' else th2d.ProfileY()
            else:
                print(f'Warning: not TH2D object!')
                return {}
    else:
        print(f'Warning: unknown input type!')
    return tprofiles

def getWF1D(adc, name, title):
    length = adc.size()
    wf1d = TH1D(name, title, length, 0, length)
    wf1d.SetXTitle('time/ns')
    wf1d.SetYTitle('adc')
    for i in range(length):
        wf1d.SetBinContent(i+1, adc[i])
    return wf1d