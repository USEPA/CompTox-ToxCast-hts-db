import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec
from textwrap import wrap
import scipy.optimize as optz
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.text as text
import matplotlib.font_manager as fm 
import pylab as pl
from matplotlib.patches import Ellipse, Circle
from collections import *
import matplotlib.cm as cm
import matplotlib as mpl
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from copy import copy
#from IPython import parallel as PP
import os
import subprocess
import mongoengine as mng
import numpy.linalg as LA

from util.misc import *
from bio.data.htsdb import *



FT1=['p53Act','StressKinase','OxidativeStress','MicrotubuleCSK','MitoMass','MitoMembPot',
    'MitoticArrest','CellCycleArrest','NuclearSize','CellLoss']
LB1=['p53 Act','Stress Kinase','Oxidative Stress','Microtubule','Mito Mass','Mito Memb Pot',
    'Mitotic Arrest','Cell Cycle Arrest','Nuclear Size','Cell Number']
FTLB1 = dict(zip(FT1,LB1))

def getAnalysisMethods():
    return [(i.name,i.desc) for i in HtsAnalysisMethod.objects]

def getChem(x):
    return [(i.name,i.eid) for i in HtsChem.objects(name__icontains=x)]
    
def getExpList():
    return [dict(exp_id=i.eid,name=i.name,cell=i.cell,org=i.org) for i in HtsExp.objects]

def getExpPlates(exp_id):
    Ex = HtsExp.objects(eid=exp_id)
    if Ex: Ex=Ex.first()
    else: return

    R = [dict(plate_id=i.eid,time=i.timeh) for i in HtsPlate.objects(exp=Ex)]

    return R

def getExpChems(exp_id,info=False,get_cid=False):
    Ex = HtsExp.objects(eid=exp_id)
    if Ex: Ex=Ex.first()
    else: return
    
    if info:
        Chems = []
        for HC in Ex.chems:
            if HC.cid:
                Chems.append(dict(chem_id=HC.cid,hchem_id=HC.eid,name=HC.name,casrn=HC.casrn))
        return Chems
        
    elif get_cid:
        CID = []
        for HC in Ex.chems:
            if HC.cid:
                CID.append(HC.cid)
        return CID
    else:
        return [hc.eid for hc in Ex.chems]

def getChemConcResp(hchem_id,assay_id=None,exp_id=None,plate_id=None,ares=False,timeh=None,
                    prows=['chem_id','hchem_id','chem_name','chem_casrn','lconc'],
                    pcols=['timeh','assay_id'],
                    add_t0=False,
                    ret=None):
    """
    ret == None -- returns a list of dicts
    ret == df   -- returns a dataframe that can be accessed as follows:-
    add_t0 just adds t=0 with slfc=0, pzsc=0, ppct=0,sm=0,raw=0 (only slfc and pzsc are useful)

    Usage:
    Get the conc-response for Captan/TX001434 in experiment APR-HepG2-PhI:
    
    DF = getChemConcResp(u'TX001434',exp_id='APR-HepG2-PhI',ret='df')
    
    Get the smoothed log2(foldchange) for 'CellLoss' at 72h:-
    DF.ix[:,'slfc'].ix[:,(72.0,'CellLoss')]
    
    Get the Z-score for 'StressKinase' at 1h:-
    DF.ix[:,'pzsc'].ix[:,(1.0,'StressKinase')]

    Get the pct change for all endpoints at 24h:-
    DF.ix[:,'ppct'].ix[:,(24.0)]

    Get the pzsc for 'MitoMembPot' across all times:-
    DF.ix[:,'pzsc'].ix[:,[(1.0,'MitoMembPot'),(24.0,'MitoMembPot'),(72.0,'MitoMembPot')]]    
    
    """
    HC = HtsChem.objects(eid=hchem_id)
    if not HC: return
    HC = HC.first()
    chem_id = HC.cid
    chem_casrn = HC.casrn
    chem_name = HC.name
    args = dict(chem=HC)
    if assay_id:
        A = HtsAssay.objects(eid=assay_id)
        if A:
            A = A.first()
            args['assay']=A

    if exp_id and plate_id:
        P = HtsPlate.objects(eid=plate_id,
                             exp__in=HtsExp.objects(eid=exp_id)).first()
        args['plate']=P
    elif plate_id:
        P = HtsPlate.objects(eid=plate_id).first()
        args['plate']=P
    elif exp_id:
        args['plate__in']=HtsPlate.objects(exp__in=HtsExp.objects(eid=exp_id))
            
    elif timeh:
        args['plate__in']=HtsPlate.objects(timeh=timeh)

    CRC = HtsConcRespCurveNrm.objects(**args)
    Res = []
    for Crc in CRC:
        R = dict(chem_id=chem_id,hchem_id=hchem_id,assay_id=Crc.assay.eid,plate_id=Crc.plate.eid,timeh=Crc.plate.timeh,chem_casrn=chem_casrn,chem_name=chem_name,
                 resp=[dict(i.to_mongo()) for i in Crc.resp])
        if ares:
            AR = HtsAssayResult.objects(crc=Crc)
            if AR:
                AR = AR.first()
                R['assay_result']=dict(hit_call=AR.hit_call,lec_mol=AR.lec_mol,eff_lfc=AR.eff_lfc)
        Res.append(R)
        
    if ret!='df': return Res

    # Construct dataframe
    RESP = []
    K1 = ['assay_id','chem_id','hchem_id','plate_id','timeh','chem_casrn','chem_name']
    for X in Res:
        x = dict(zip(K1,[X[k] for k in K1]))
        for r in X['resp']:
            r.update(x)
            RESP.append(r)
    RESP_df = pd.DataFrame(RESP)

    if add_t0:
        TM = RESP_df.timeh.unique()
        FT = RESP_df.assay_id.unique()
        LC = RESP_df.lconc.unique()
        RESP0 = []
        K1 = ['chem_id','hchem_id','plate_id','chem_casrn','chem_name']
        x = dict(zip(K1,[RESP[0][k] for k in K1]))
        for tm in TM:
            for lc in LC:
                for ft in FT:
                    x0 = copy(x)
                    x0.update(dict(assay_id=ft,lconc=lc,timeh=0.0,ppct=0.0,pzsc=0.0,slfc=0.0,sm=0.0,raw=0.0))
                    RESP0.append(x0)
        RESP_df = pd.DataFrame(RESP+RESP0)
    #return RESP,RESP0                      
    return pd.pivot_table(RESP_df,
                          index=prows,
                          columns=pcols,
                          values=['ppct','pzsc','raw','slfc','sm'],
                          aggfunc=np.average,fill_value=0.0)

def getChemAssayResults(hchem_id,exp_id=None,plate_id=None,meth_id=None,ret='df'):
    """
    ret == None -- returns a list of dicts
    ret == df   -- returns a dataframe that can be accessed as follows:-

    Usage:
    Get the assay results for Captan/TX001434 in experiment APR-HepG2-PhI:
    
    DF = getChemAssayResults(u'TX001434',exp_id='APR-HepG2-PhI',ret='df')
    
    Get the hit calls for 72h:-
    DF.ix[:,'hit_call'].ix[:,72.0]

    Get the log2(fc) efficacy for CellLoss:-
    DF['lfc_eff'].ix[:,[(1.0,'CellLoss'),(24.0,'CellLoss'),(72.0,'CellLoss')]]
    
    """
    HITS=[]
    Meth=None
    # Get the chemical
    HC = HtsChem.objects(eid=hchem_id)
    if not HC: return
    HC = HC.first()
    cid = HC.cid
    if not cid:
        cid = getCidFromHcid(HC)
    H = []

    args=dict(chem=HC)

    if meth_id:
        Meth=HtsAnalysisMethod.objects(name=meth_id)
        if Meth: args['meth__in']=Meth

    if exp_id and plate_id:
        P = HtsPlate.objects(eid=plate_id,
                             exp__in=HtsExp.objects(eid=exp_id))
        args['crc__in']=HtsConcRespCurveNrm.objects(plate__in=P)
    elif plate_id:
        P = HtsPlate.objects(eid=plate_id)
        args['crc__in']=HtsConcRespCurveNrm.objects(plate__in=P)
    elif exp_id:
        args['exp__in']=HtsExp.objects(eid=exp_id)


    for AR in HtsAssayResult.objects(**args):
        HITS.append(dict(chem_id=cid,hchem_id=HC.eid,chem_name=HC.name,chem_casrn=HC.casrn,
                         meth_id=meth_id,
                         assay_id=AR.assay.eid,timeh=AR.crc.plate.timeh,
                         exp_id=AR.crc.plate.exp.eid,plate_id=AR.crc.plate.eid,
                         hit_call=AR.hit_call,lec=AR.lec,llec=AR.llec,
                         lfc_eff=AR.lfc_eff,pc_eff=AR.pc_eff,llc_eff=AR.llc_eff))

    if ret=='df':
        return pd.pivot_table(pd.DataFrame(HITS),
                              index=['chem_id','hchem_id','chem_name','chem_casrn'],
                              columns=['timeh','assay_id'],
                              values=['hit_call','lec','llec','lfc_eff','pc_eff','llc_eff'],
                              aggfunc=np.average,
                              fill_value=0)
    else:
        return HITS


    
def plotHtsConcRespHM(hchem_id,exp_id,add_t0=False,use_resp='fc',
                      fgsz=(20,4),vmin=-3,vmax=3,fs=0.9,xyoff=[2,-6],ctrl=False,iconc=3,
                      FT=FT1,
                      FTLB=FTLB1,
                      lec=False,
                      CM=cm.RdYlBu_r,draw_chem=False,
                      save=None,loc=None,cb=False,sub_title=None):
    # Get the conc-resp data
    X_CR = getChemConcResp(hchem_id,exp_id=exp_id,ret='df',add_t0=add_t0)
    if X_CR.shape[0]==0: return
    # Get the fold change
    resp_col=None
    leg_txt =None
    if use_resp=='zs':
        resp_col='pzsc'
        leg_txt ='Z'
    elif use_resp=='fc':
        resp_col='slfc',
        leg_txt =r'$log_2(FC)$'
    elif use_resp=='pct':
        resp_col='ppct'
        leg_txt ='Prcnt'        
    else:
        resp_col='slfc'
        leg_txt =r'$log_2(FC)$'

    X_FC = X_CR[resp_col]

    # Get LEC
    LEC = None
    EFF = None
    X_AR = getChemAssayResults(hchem_id,exp_id=exp_id,ret='df')
    if X_AR.shape[0]==0: return
    try:
        LEC = X_AR['llec']
        EFF = X_AR['lfc_eff']
    except:
        print "No LEC data"

    if not lec: LEC=None
    fig = pl.figure(figsize=fgsz)
    FT0 = set([i[1] for i in X_FC.columns])
    TM0 = sorted(set([i[0] for i in X_FC.columns]))
    #gs1 = gridspec.GridSpec(2,6,width_ratios=[1,1,1,1,1,1])
    #if ctrl:
    #    gs1 = gridspec.GridSpec(1,8,width_ratios=[1,1,1,1,1,1,1,1])
    #else:
    gs1 = gridspec.GridSpec(2,8,width_ratios=[1,1,1,1,1,1,1,0.3])
    
    # Chemical name casrn
    # ['chem_id', 'hchem_id', 'chem_name', 'chem_casrn', 'lconc']
    ax = pl.subplot(gs1[0:2,0:2])
    chem = X_FC.index[0][2]
    casrn= X_FC.index[0][3]
    sid  = X_FC.index[0][1]
    cid =  X_FC.index[0][0]
    nm = '\n'.join(wrap(chem,30)) + '\n' + casrn + '\n' + sid
    if sub_title: nm += '\n' + sub_title
        
    txt = text.Text(60,100,nm,ha='center',va='center',color='#111111',
                    fontproperties=fm.FontProperties(size=10))
    ax.add_artist(txt)    
    ax.set_axis_off()
    ax.set_xlim(right=120)
    ax.set_ylim(top=120)
        
    k=0
    iconc=-1
    for ft in FT:
        if ft not in FT0: continue
        X = X_FC[[(tm,ft) for tm in TM0]]
        TM = ['%dh' % tm for tm in TM0]
        i = k % 5 + 2
        if k<5:
            j=0
        else:
            j=1
        ax = pl.subplot(gs1[j,i])
        x=ax.pcolor(X.T,cmap=CM,edgecolors='white',lw=1,vmin=vmin,vmax=vmax)
        ax.set_title(FTLB[ft])
        # Calc conc
        LC = [ii[iconc] for ii in X.index]
        uM_log = np.array(LC)
        uM     = np.round(10**np.array(uM_log + 6,dtype=np.float),decimals=2)
        
        try:
            LLC = [v[iconc] for v in X_FC.index]
            LLCx= np.arange(len(LLC))+0.5
            
            for it,tm in enumerate(TM0):
                if tm==0: continue
                lec_ft = LEC[(tm,ft)]
                if lec_ft==0 or np.isnan(lec_ft): continue
                lec_fc = EFF[(tm,ft)]
                col=ifthen(lec_fc<0,'red','blue')
                ilc=np.interp(lec_ft,LLC,LLCx)
                ax.text(ilc,it+0.5,r'$\bullet$',ha='center',va='center',color=col,
                        fontproperties=fm.FontProperties(size=20))
        except:
            pass
        #print i
        #if i == 0:
        #print i,j,ft
        if j==1 or ctrl:
            ax.xaxis.tick_bottom()
            ax.set_xticks(np.arange(len(uM))+0.5, minor=False)
            ax.set_xticklabels(uM,rotation=90)
            for tick in ax.get_xticklines(): tick.set_visible(True)
            for tick in ax.get_xticklabels(): tick.set_fontsize(9)
        else:
            for tick in ax.get_xticklines(): tick.set_visible(True)
            for lab in ax.get_xticklabels(): lab.set_visible(False)
    
        if i == 2:
            ax.yaxis.tick_left()
            ax.set_ylim(top=len(TM))
            ax.set_yticks(np.arange(len(TM))+0.5, minor=False)
            ax.set_yticklabels(TM)
            for tick in ax.get_yticklines(): tick.set_visible(True)
            for tick in ax.get_yticklabels(): 
                tick.set_fontsize(10)
        else:
            for tick in ax.get_yticklines(): tick.set_visible(False)
            for lab in ax.get_yticklabels(): lab.set_visible(False)
        k+=1
        
    ax = pl.subplot(gs1[:,7])
    ax.set_axis_off()
    if cb: 
        cbr = pl.colorbar(x,ax=ax,shrink=0.9,ticks=range(-vmax,vmax+1))
        cbr.ax.set_xlabel(leg_txt,fontsize=10)        
    if save and loc:
        fig.savefig(loc+"cr-%s-%s-%s-%s.svg" % (resp_col,exp_id.lower(),hchem_id.lower(),chem.lower().replace(' ','-')))
                
    pl.subplots_adjust(wspace=0.1,hspace=0.3,top=0.9,bottom=0.1)


def plotHtsTrajHM(hchem_id,exp_id,add_t0=False,use_resp='fc',
                  fgsz=(18,4),vmin=-3,vmax=3,fs=0.9,xyoff=[2,-6],ctrl=False,iconc=3,
                  FT=FT1,
                  FTLB=FTLB1,
                  LCi=None,draw_chem=False,
                  CM=cm.RdYlBu_r,suffix=1,
                  save=None,loc=None,cb=False,sub_title=None):
    # Get the conc-resp data
    X_CR = getChemConcResp(hchem_id,exp_id=exp_id,ret='df',add_t0=add_t0,
                           prows=['chem_id','hchem_id','chem_name','chem_casrn','timeh','lconc'],
                           pcols=['assay_id'])
    if X_CR.shape[0]==0: return

    resp_col=None
    leg_txt =None
    if use_resp=='zs':
        resp_col='pzsc'
        leg_txt ='Z'
    elif use_resp=='fc':
        resp_col='slfc',
        leg_txt =r'$log_2(FC)$'
    elif use_resp=='pct':
        resp_col='ppct'
        leg_txt ='Prcnt'        
    else:
        resp_col='slfc'
        leg_txt =r'$log_2(FC)$'

    X_FC = X_CR.ix[:,resp_col]

    # Get LEC
    LEC = None
    EFF = None
    #X_AR = getChemAssayResults(hchem_id,exp_id=exp_id,ret='df')
    #if X_AR.shape[0]==0: return
    #try:
    #    LEC = X_AR['llec']
    #    EFF = X_AR['lfc_eff']
    #except:
    #    print "No LEC data"
    
    fig = pl.figure(figsize=fgsz)
    FT0 = set([i for i in X_FC.columns])
    TM0 = sorted(set([i[-2] for i in X_FC.index]))
    LC0 = list(sorted(set([i[-1] for i in X_FC.index]))) 
    if LCi:
        LC0 = [LC0[i] for i in LCi]
        
    #gs1 = gridspec.GridSpec(2,6,width_ratios=[1,1,1,1,1,1])
    #if ctrl:
    #    gs1 = gridspec.GridSpec(1,8,width_ratios=[1,1,1,1,1,1,1,1])
    #else:

    gs1 = None
    if len(LC0)>5:
        gs1 = gridspec.GridSpec(2,8,width_ratios=[1,1,1,1,1,1,1,0.3])
        ax = pl.subplot(gs1[0:1,0:1])
    else:
        gs1 = gridspec.GridSpec(1,7,width_ratios=[1.2,1,1,1,1,1,0.3])
        ax = pl.subplot(gs1[0,0])
    
    # Chemical name casrn
    # ['chem_id', 'hchem_id', 'chem_name', 'chem_casrn', 'lconc']
    chem = X_FC.index[0][2]
    casrn= X_FC.index[0][3]
    sid  = X_FC.index[0][1]
    cid =  X_FC.index[0][0]
    nm = None
    if len(LC0)>5:
        nm = '\n'.join(wrap(chem,30)) + '\n' + casrn + '\n' + sid
        if sub_title: nm += '\n' + sub_title
    else:
        nm = '\n'.join(wrap(chem,25))
        
    txt = text.Text(50,50,nm,ha='center',va='center',color='#111111',
                    fontproperties=fm.FontProperties(size=10))
    ax.add_artist(txt)    
            
    ax.set_axis_off()
    ax.set_xlim(0,100)
    ax.set_ylim(0,100)
    
        
    k=0
    iconc=-1

    for lc in LC0:
        uM  = np.round(10**(lc + 6),decimals=2)
        X = X_FC[FT].select(lambda xx: xx[iconc]==lc,axis=0)
        #print lc,X_FC.shape,X.shape
        TM = ['%dh' % tm for tm in TM0]
        i = k % 5 + 1
        if k<5:
            j=0
        else:
            j=1
        ax = pl.subplot(gs1[j,i])
        
        x=ax.pcolor(X.as_matrix(),cmap=CM,edgecolors='white',lw=1,vmin=vmin,vmax=vmax)

        ax.set_title(r'%5.2f$\mu$M' % uM)
        
        if len(LC0)<=5 or j==1:
            ax.xaxis.tick_bottom()
            ax.set_xticks(np.arange(len(FT0))+0.5, minor=False)
            ax.set_xticklabels([FTLB[ft] for ft in X.columns] ,rotation=90)
            for tick in ax.get_xticklines(): tick.set_visible(True)
            for tick in ax.get_xticklabels(): tick.set_fontsize(9)
        else:
            for tick in ax.get_xticklines(): tick.set_visible(True)
            for lab in ax.get_xticklabels(): lab.set_visible(False)
    
        if (len(LC0)>5 and i == 2) or i==1:
            ax.yaxis.tick_left()
            ax.set_ylim(top=len(TM))
            ax.set_yticks(np.arange(len(TM))+0.5, minor=False)
            ax.set_yticklabels(TM)
            for tick in ax.get_yticklines(): tick.set_visible(True)
            for tick in ax.get_yticklabels(): 
                tick.set_fontsize(10)
        else:
            for tick in ax.get_yticklines(): tick.set_visible(False)
            for lab in ax.get_yticklabels(): lab.set_visible(False)
        k+=1
        
    if len(LC0)>5:
        ax = pl.subplot(gs1[:,7])
    else:
        ax = pl.subplot(gs1[:,6])
        
    ax.set_axis_off()
    if cb: 
        cbr = pl.colorbar(x,ax=ax,shrink=0.9,ticks=range(-vmax,vmax+1))
        cbr.ax.set_xlabel(leg_txt,fontsize=10)        

    pl.subplots_adjust(wspace=0.3,hspace=0.3,top=1,bottom=0)
    if save and loc:
        fig.savefig(loc+"traj-hm-%s-%s-%s-%s-%d.svg" % (resp_col,exp_id.lower(),hchem_id.lower(),
                                                        chem.lower().replace(' ','-'),suffix))
        fig.savefig(loc+"traj-hm-%s-%s-%s-%s-%d.png" % (resp_col,exp_id.lower(),hchem_id.lower(),
                                                        chem.lower().replace(' ','-'),suffix))
