import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict
from os import listdir
from matplotlib.colors import LinearSegmentedColormap, ColorConverter

supplpath = '/casa/bsc/projects/1_DCS/2004_paper_prep/supplementary'

## encode miRBase structure
def load_mirbase_str():
    strfile = '%s/mirbase-v21.str' % supplpath
    strinfos = defaultdict(list)
    mir = ''
    for l in open(strfile, 'rt'):
        if l.startswith('>'):
            mir = l.split()[0][1:]
        else:
            strinfos[mir].append(l.replace('\n',''))
    return strinfos


def concat_seq(infos, nucs):
    ss5, ds5, match, ds3, ss3 = infos
    seq = ''
    if ds5[0] in nucs:
        ss5 = '- '+ss5
    if ds3[0] in nucs:
        ss3 = '- '+ss3
        
    for s5, d5 in zip(ss5.split(), ds5.split()):
        seq += s5+d5
    seq += match[-1]
    for d3, s3 in zip(ds3.split()[::-1], ss3.split()[::-1]):
        seq += d3[::-1]+s3[::-1]    
    seq = seq.replace('-','').replace('|','').replace(' ','')
    return seq


def get_new_str(mir):
    nucs = 'ACGUacgu'
    ss5, ds5, match, ds3, ss3 = strinfos[mir][1:6]
    newinfo = [ '', '', '', '', '' ]
    
    for ss, ssnew, oppo, index in [ (ss5,'',ss3,0), (ss3,'',ss5,4) ]:
        for i,s in enumerate(ss):
            if s in nucs and i<match.find('|'): ssnew += 'F'
            elif s in nucs and i>match.rfind('|'): ssnew += 'L'
            elif s in nucs and oppo[i] in nucs: ssnew += 'S'
            elif s in nucs: ssnew += 'A'
            else: ssnew += s
        newinfo[index] = ssnew
    
    for ds, dsnew, index in [ (ds5,'',1), (ds3,'',3) ]:
        for i,s in enumerate(ds):
            if s in nucs and i>match.rfind('|'): dsnew += 'L'
            elif s in nucs: dsnew += 'M'
            else: dsnew += s
        newinfo[index] = dsnew
    
    if match[-1] in nucs: newinfo[2] = match[:-1] + 'L'
    else: newinfo[2] = match   
    return concat_seq(newinfo, 'FLSAM')


## Find pre-miRNA annotation
def parse_line(row):
    for s in row['attr'].split(';'):
        tag, value = s.split('=')
        row[tag] = value
    return row


def load_mirbase_annot():
    gff = '%s/human_mirbase-v21.gff3' % supplpath
    anntbl = pd.read_table(gff, header=12, sep='\t', usecols=[0,2,3,4,6,8], 
                           names=['chr','kind','start','end','strand','attr'])
    anntbl = anntbl.apply(parse_line, axis=1).drop(['attr'], axis=1)
    return anntbl

strinfos = load_mirbase_str()
anntbl = load_mirbase_annot()
annpri = anntbl[anntbl['kind']=='miRNA_primary_transcript'].set_index('Name')
annmat = anntbl[anntbl['kind']=='miRNA'].set_index('ID')
prif = '%s/hairpin_mirbase-v21.fa' % supplpath
hpnseqs = { s.id:str(s.seq) for s in SeqIO.parse(prif, 'fasta') }


def get_pri_mature():
    primat = {}
    for pri, row in annpri.iterrows():
        matureDerives = anntbl[anntbl['Derives_from']==row['ID']]
        primat[pri] = matureDerives['Name'].tolist()
    return primat


def count_len(strt):
    return strt.count('M')+strt.count('S')


def get_pre_seq(pri, relativepos, arm):
    overhang3 = 2
    priseq = hpnseqs[pri]
    pristr = get_new_str(pri)
    if arm=='5p':
        preend = [ i for i in range(len(pristr)+1) 
            if count_len(pristr[i:])<=max(0,count_len(pristr[:relativepos-1])-overhang3) ][0]
        return priseq[relativepos-1:preend]
    prestart = [ i+1 for i in range(len(pristr))
                 if count_len(pristr[:i])>=count_len(pristr[relativepos:])+overhang3 ][0]
    return priseq[prestart-1:relativepos]
    

def get_pre_annot(pri):
    chrom, strand, priid = annpri.loc[pri, ['chr','strand','ID']]
    matids = list(annmat[annmat['Derives_from']==priid].index)
    if len(matids)==2:
        start = min(annmat.loc[matids, 'start'])
        end = max(annmat.loc[matids, 'end'])
    elif strand=='+':
        pristart, priend = annpri.loc[pri, ['start','end']]
        matstart, matend = annmat.loc[matids[0], ['start','end']]
        if (matstart-pristart)<(priend-matend): # 5p
            start = matstart
            end = matstart + len(get_pre_seq(pri, matstart-pristart+1, '5p')) - 1
        else: # 3p
            end = matend
            start = matend - len(get_pre_seq(pri, matend-pristart+1, '3p')) + 1
    else: # strand=='-'
        pristart, priend = annpri.loc[pri, ['start','end']]
        matstart, matend = annmat.loc[matids[0], ['start','end']]
        if (matstart-pristart)<(priend-matend): # 3p
            start = matstart
            end = matstart + len(get_pre_seq(pri, priend-matstart+1, '3p')) - 1
        else: # 5p
            end = matend
            start = matend - len(get_pre_seq(pri, priend-matend+1, '5p')) + 1
    return chrom, start, end, strand


def get_phylop_scores(chrom, start, end):
    phystarts = [ int(f.split('.')[0]) for f in listdir('%s/phyloP100way/%s'%(supplpath, chrom)) 
                  if f.endswith('.gz') ]
    init = sorted([ s for s in phystarts if s<=start ])[-1]
    nextst = sorted([ s for s in phystarts if s>start])[0]
    rels, rele = start-init, end-init
    phys = open('%s/phyloP100way/%s/%s.phylop.gz'%(supplpath,chrom,init),'rt').read().strip().split('\n')
    if len(phys) <= rels:
        if nextst > end:
            return ['-']*(end-start+1)
        else:
            return ['-']*(nextst-start) + get_phylop_scores(chrom, nextst, end)
    elif len(phys) > rele:
        return map(float, phys[rels:rele+1])
    else: # rels < len(phys) <= rele
        return map(float, phys[rels:]) + get_phylop_scores(chrom, init+len(phys), end)


def calculate_dist(pair1, pair2):
    return min(abs(pair1[0]-pair2[0]), abs(pair1[1]-pair2[1]))


def list_apical_jcs(infos, minss, init): # [ (1, 124), (2, 123), ... ]
    jcpos, dist = [], init
    for uppair, lowpair in zip(infos[:-1], infos[1:]):
        jcsize = abs(uppair[0]-lowpair[0])+abs(lowpair[1]-uppair[1])-2
        if jcsize>=minss:
            jcpos.append((dist,uppair,jcsize))
        dist += calculate_dist(uppair,lowpair)
    jcpos.append((dist,infos[-1],infos[-1][1]-infos[-1][0]))
    return jcpos


def find_apical_junction(mir, minss, ps, pe, constseq, strinfo, optipos):
    hpnseq50 = hpnseqs[mir].replace('U','T')[:50]
    hpnstr = get_new_str(mir)
    if constseq.find(hpnseq50)>=0:
        loopmid = constseq.find(hpnseq50)+hpnstr.find('L')+hpnstr.count('L')//2
    else:
        loopmid = hpnstr.find('L')-hpnseq50.find(constseq[:10])+hpnstr.count('L')//2
    
    ustems = [ (s,e) for s,e in strinfo if ps<=s<loopmid<e<=pe ]
    if not ustems:
        return (0, (0,0))
    
    us1st5, us1st3 = ustems[0]
    lstems = [ (s,e) for s,e in strinfo if s<us1st5 and e>us1st3 ]
    if not lstems:
        init = 1-2*int(ps==0)
    else:
        ls1st5, ls1st3 = lstems[-1]
        if ps:
            init = min(us1st5-ps, ls1st3-us1st3-1)+1
        else:
            init = min(pe-us1st3, us1st5-ls1st5-1)-1
    jcposes = list_apical_jcs(ustems, minss, init)
    return min(jcposes, key=lambda x: (-x[2],abs(x[0]-optipos)))[:2]


def list_basal_jcs(infos, minss, lastunp, init): # [ (1, 124), (2, 123), ... ]
    jcpos, dist = [], init
    for uppair, lowpair in zip(infos[:-1], infos[1:]):
        jcsize = abs(uppair[0]-lowpair[0])+abs(lowpair[1]-uppair[1])-2
        if jcsize>=minss:
            jcpos.append((dist,uppair,jcsize))
        dist += calculate_dist(uppair, lowpair)
    jcpos.append((dist,infos[-1],lastunp))
    return jcpos


def find_basal_junction(mir, strinfo, minss, ajst, lastunp, optipos, init):
    stems = [ (s,e) for s,e in strinfo if s<=ajst<e ][::-1]
    jcposes = list_basal_jcs(stems, minss, lastunp, init)
    if jcposes:
        return min(jcposes, key=lambda x:(-x[2],abs(x[0]-optipos)))[:2]
    return (0, (0,0))


to_rgb = ColorConverter().to_rgb
def custom_cmap(colors, cmap_name="newmap", nspace=3, linear=True):
    if (type(colors) is str) or (len(colors) is 1):
        colors = [colors, "white"]
    ncolors = len(colors)
    sidx = map(int, map(np.around, np.linspace(0, nspace-1, num=ncolors)))
    intervals = np.linspace(0, 1.0, num=nspace)
    rgb = ["red", "green", "blue"]
    cdict = {e:None for e in rgb}
    for element, components in zip(rgb, zip(*[to_rgb(c) for c in colors])):
        intensities = [components[0]]
        for i, value in enumerate(components):
            if i + 1 == len(components): break
            v1, v2 = components[i:i+2]
            intensities += list(np.linspace(v1, v2, num=sidx[i+1] - sidx[i] + 1))[1:]
        cdict[element] =  zip(intervals, intensities, intensities)
    return LinearSegmentedColormap(cmap_name, cdict)


def adjust_ct(strinfo, ajs, aje):
    newinfo = [ 0 if max(i,pair)<ajs or min(i,pair)>aje or (ajs<i<aje and ajs<pair<aje)
                else pair for i,pair in strinfo ]
    return zip(range(1,126), newinfo)



