def min_more(num, scale=0.1):
    if isinstance(num, (int, float)):
        return (1-scale)*num if num>0 else (1+scale)*num
    elif isinstance(num, (tuple, list)):
        return min_more(min(num))
    else:
        print(f'Warning min_more: unknown input type!')

def max_more(num, scale=0.1):
    if isinstance(num, (int, float)):
        return (1+scale)*num if num>0 else (1-scale)*num
    elif isinstance(num, (tuple, list)):
        return max_more(max(num))
    else:
        print(f'Warning max_more: unknown input type!')

def minmax_more(nmin, nmax, scale=0.1):
    if (isinstance(nmin, (int, float)) and isinstance(nmax, (int, float))) or (isinstance(nmin, (tuple, list)) and isinstance(nmax, (tuple, list))):
        tmpmin = min_more(nmin)
        tmpmax = max_more(nmax)
        if tmpmin>0 and tmpmin<scale*tmpmax:
            return (0, tmpmax)
        elif tmpmax<0 and tmpmax>scale*tmpmin:
            return (tmpmin, 0)
        else:
            return (tmpmin, tmpmax)
    else:
        print(f'Warning minmax_more: unknown input type!')

def MC_merge_maxQT(raw_HT, raw_PE, merge_cut = 20):
    if isinstance(raw_HT, (list, tuple)):
        N = len(raw_HT)
        if not N:
            return 0

        index = sorted(range(len(raw_HT)), key=lambda k: raw_HT[k])
        in_ht = tuple(raw_HT[i] for i in index)
        in_pe = tuple(raw_PE[i] for i in index)
        # index = np.argsort(np.array(raw_HT))
        # in_ht = np.array(raw_HT)[index]
        # in_pe = np.array(raw_PE)[index]
        out_ht= [in_ht[0]]
        out_pe= [in_pe[0]]
        for i in range(1,len(in_ht)):
            if in_ht[i] - in_ht[i-1] < merge_cut:
                out_pe[-1] += in_pe[i]
            else:
                out_ht.append(in_ht[i])
                out_pe.append(in_pe[i])
        # return out_ht[np.argmax(out_pe)]
        return out_ht[out_pe.index(max(out_pe))]

def fitGausL(h, fitRange=4, xrange=None):
    if not xrange:
        h.Fit('gaus','Q')
        f = h.GetFunction('gaus')
        pars = tuple( f.GetParameter(i) for i in (1,2) )
        del f

        xmin = pars[0] - fitRange*pars[1]
        xmax = pars[0] + fitRange*pars[1]
        xrange = (xmin, xmax)

    h.Fit('gaus','LQ','',xrange[0],xrange[1])
    f = h.GetFunction('gaus')
    pars = tuple( f.GetParameter(i) for i in (1,2) )
    return pars

