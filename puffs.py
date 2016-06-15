# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 12:13:09 2016

@author: Kyle Ellefsen
"""
from __future__ import (absolute_import, division,print_function) #, unicode_literals) http://stackoverflow.com/a/28162261/4441139
import numpy as np
from PyQt4.QtCore import Qt
from PyQt4.QtGui import qApp
import pyqtgraph as pg
from .gaussianFitting import fitGaussian, fitRotGaussian


def scatterRemovePoints(scatterplot,idxs):
    i2=[i for i in np.arange(len(scatterplot.data)) if i not in idxs]
    points=scatterplot.points()
    points=points[i2]
    spots=[{'pos':points[i].pos(),'data':points[i].data(),'brush':points[i].brush()} for i in np.arange(len(points))]
    scatterplot.clear()
    scatterplot.addPoints(spots)
    
def scatterAddPoints(scatterplot,pos,data):
    points=scatterplot.points()
    spots=[{'pos':points[i].pos(),'data':points[i].data()} for i in np.arange(len(points))]
    spots.extend([{'pos':pos[i],'data':data[i]} for i in np.arange(len(pos))])
    scatterplot.clear()
    scatterplot.addPoints(spots)
    
class Puffs:
    def __init__(self,clusters,cluster_im,puffAnalyzer,persistentInfo=None):#weakfilt,strongfilt,paddingXY,paddingT_pre,paddingT_post,maxSigmaForGaussianFit,rotatedfit):
        self.puffAnalyzer=puffAnalyzer
        self.udc=puffAnalyzer.udc        
        self.puffs=[]
        self.index=0
        self.clusters=clusters
        self.normalized_window=puffAnalyzer.normalized_window
        self.data_window=puffAnalyzer.data_window
        self.cluster_im=cluster_im
        
        self.puffs=[]
        nClusters=len(self.clusters.clusters)
        for i in np.arange(nClusters):
            percent=i/nClusters
            self.puffAnalyzer.algorithm_gui.gaussianProgress.setValue(percent*100); qApp.processEvents();
            self.puffs.append(Puff(i,self.clusters,self,persistentInfo))
        self.puffAnalyzer.algorithm_gui.gaussianProgress.setValue(100); qApp.processEvents();

    def __getitem__(self, item):
        if len(self.puffs)>0:
            return self.puffs[item]
        else:
            return None
    def removeCurrentPuff(self):
        del self.puffs[self.index]
        if self.index==0:
            return self.index
        else:
            self.index-=1
        return self.index
    def getPuff(self):
        if len(self.puffs)>0:
            return self.puffs[self.index]
        else:
            return None
    def increment(self):
        self.index+=1
        if len(self.puffs)<self.index+1:
            self.index=0
    def decrement(self):
        self.index-=1
        if self.index<0:
            self.index=len(self.puffs)-1
    def setIndex(self,index):
        self.index=index
        if len(self.puffs)<self.index+1:
            self.index=0
        elif self.index<0:
            self.index=len(self.puffs)-1
    def removePuffs(self,puffs):
        idxs=[]
        for puff in puffs:
            idxs.append([point['data'] for point in self.puffAnalyzer.s1.data].index(puff))
            self.puffs.remove(puff)
        scatterRemovePoints(self.puffAnalyzer.s1,idxs)
        if self.index>=len(self.puffs):
            self.index=len(self.puffs)-1
            
    def addPuffs(self,puffs):
        s=self.puffAnalyzer.s1
        self.puffs.extend(puffs)
        pos=[[puff.kinetics['x'],puff.kinetics['y']] for puff in puffs]
        scatterAddPoints(s,pos,puffs)
        self.puffAnalyzer.updateScatter()
        #s.addPoints(pos=pos,data=puffs)

  
   
class Puff:
    def __init__(self,starting_idx,clusters,puffs,persistentInfo=None):
        print('Creating event {}/{}'.format(starting_idx, len(clusters.clusters)-1))
        self.starting_idx=starting_idx
        self.clusters=clusters
        self.puffs=puffs
        self.udc=puffs.udc
        self.color=(255,0,0,255)
        self.originalbounds=self.clusters.bounds[starting_idx] # 2x3 array: [[t_min,x_min,y_min],[t_max,x_max,y_max]]
        t0=self.originalbounds[0][0]
        t1=self.originalbounds[1][0]
        x0=self.originalbounds[0][1]-self.udc['paddingXY']
        x1=self.originalbounds[1][1]+self.udc['paddingXY']
        y0=self.originalbounds[0][2]-self.udc['paddingXY']
        y1=self.originalbounds[1][2]+self.udc['paddingXY']
        mt,mx,my=self.puffs.data_window.image.shape
        if t0<0: t0=0
        if y0<0: y0=0
        if x0<0: x0=0
        if t1>=mt: t1=mt-1
        if y1>=my: y1=my-1
        if x1>=mx: x1=mx-1
        self.bounds=[(t0,t1),(x0,x1),(y0,y1)]
        
        if persistentInfo is not None:
            puff=persistentInfo.puffs[starting_idx]
            self.trace=puff['trace']
            self.kinetics=puff['kinetics']
            self.gaussianParams=puff['gaussianParams']
            self.mean_image=puff['mean_image']
            self.gaussianFit=puff['gaussianFit']
            try:
                self.color=puff.color #(255,0,0,255)
            except:
                pass
            return None
        self.trace=None
        self.kinetics=dict()
        
        #######################################################################
        #############          FIND (x,y) ORIGIN       ########################
        #######################################################################
        '''
        For debugging, use the following code:
        self=g.m.puffAnalyzer.puffs.getPuff()
        from plugins.detect_puffs.threshold_cluster import *
        [(t0,t1),(x0,x1),(y0,y1)]=self.bounds
        mt,mx,my=self.puffs.normalized_window.image.shape
        '''
        
        bb=self.bounds
        before=bb[0][0]-self.udc['paddingT_pre']
        after=bb[0][1]+self.udc['paddingT_post']
        if before<0: before=0
        if after>=mt: after=mt-1;
        self.kinetics['before']=before
        self.kinetics['after']=after
        self.sisterPuffs=[] # the length of this list will show how many gaussians to fit 
        for idx,cluster in enumerate(self.clusters.bounds):
            if np.any(np.intersect1d(np.arange(cluster[0,0],cluster[1,0]),np.arange(t0,t1))):
                if np.any(np.intersect1d(np.arange(cluster[0,1],cluster[1,1]),np.arange(x0,x1))):
                    if np.any(np.intersect1d(np.arange(cluster[0,2],cluster[1,2]),np.arange(y0,y1))):
                        if idx != self.starting_idx:
                            self.sisterPuffs.append(idx)
        I=self.puffs.normalized_window.image[bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1,bb[2][0]:bb[2][1]+1]
        I=np.mean(I,0)
        
        def getFitParams(idx):
            xorigin,yorigin=self.clusters.origins[idx,1:]-np.array([x0,y0])
            sigma=self.clusters.standard_deviations[idx]
            x_lower=xorigin-sigma; x_upper=xorigin+sigma; y_lower=yorigin-sigma; y_upper=yorigin+sigma
            amplitude=np.max(I)/2
            sigma=3
            if self.udc['rotatedfit']:
                sigmax=sigma
                sigmay=sigma
                angle=45
                p0=(xorigin,yorigin,sigmax,sigmay,angle,amplitude)
                #                 xorigin                   yorigin             sigmax, sigmay, angle,    amplitude
                fit_bounds = [(x_lower,x_upper), (y_lower,y_upper),  (2,self.udc['maxSigmaForGaussianFit']), (2,self.udc['maxSigmaForGaussianFit']), (0,90),   (0,np.max(I))]
            else:
                p0=(xorigin,yorigin,sigma,amplitude) 
                #                 xorigin                   yorigin            sigma    amplitude
                fit_bounds = [(x_lower,x_upper), (y_lower,y_upper),  (2,self.udc['maxSigmaForGaussianFit']),    (0,np.max(I))] #[(0.0, 2*self.paddingXY), (0, 2*self.paddingXY),(0,10),(0,10),(0,90),(0,5)]
            return p0, fit_bounds
        p0,fit_bounds=getFitParams(self.starting_idx)
        for puff in self.sisterPuffs:
            sister_p0, sister_fit_bounds=getFitParams(puff)
            p0=p0+sister_p0
            fit_bounds=fit_bounds+sister_fit_bounds
        if self.udc['rotatedfit']:
            p, I_fit, I_fit2= fitRotGaussian(I,p0,fit_bounds,nGaussians=1+len(self.sisterPuffs))
            self.mean_image=I
            self.gaussianFit=I_fit2
            p[0]=p[0]+self.bounds[1][0] #Put back in regular coordinate system.  Add back x
            p[1]=p[1]+self.bounds[2][0] #add back y 
            self.gaussianParams=p
            xorigin,yorigin,sigmax,sigmay,angle,amplitude=self.gaussianParams
        else:
            
            p, I_fit, I_fit2= fitGaussian(I,p0,fit_bounds,nGaussians=1+len(self.sisterPuffs))
            self.mean_image=I
            self.gaussianFit=I_fit2
            p[0]=p[0]+self.bounds[1][0] #Put back in regular coordinate system.  Add back x
            p[1]=p[1]+self.bounds[2][0] #add back y 
            self.gaussianParams=p
            xorigin,yorigin,sigma,amplitude=self.gaussianParams

        if self.udc['rotatedfit']:
            self.kinetics['sigmax']=sigmax; self.kinetics['sigmay']=sigmay; self.kinetics['angle']=angle
        else:
            self.kinetics['sigma']=sigma;
        self.kinetics['x']=xorigin; self.kinetics['y']=yorigin;
        #######################################################################
        #############          FIND PEAK       ########################
        #######################################################################
        if amplitude == 0:
            I_norm = np.zeros(self.gaussianFit.shape)
            I_norm2 = np.zeros(self.gaussianFit.shape)
        else:
            I_norm = I_fit/np.sum(I_fit)
            I_norm2 = I_fit2/np.sum(I_fit2)
            

        #I=self.puffs.highpass_im[before:after+1,bb[1][0]:bb[1][1]+1,bb[2][0]:bb[2][1]+1]
        I=self.puffs.data_window.image[before:after+1,bb[1][0]:bb[1][1]+1,bb[2][0]:bb[2][1]+1]
        #baseline=np.mean(I[:,I_norm2<.01])     
        trace=np.zeros((len(I)))
        x=int(np.floor(xorigin))-self.bounds[1][0]
        y=int(np.floor(yorigin))-self.bounds[2][0]
        roi_width=self.udc['roi_width']
        r=(roi_width-1)/2        
        
        bounds=[x-r,x+r+1,y-r,y+r+1]
        if I[0,x-r:x+r+1,y-r:y+r+1].size==0: # check if roi exceeds the region we cut out of the window
            if bounds[0]<0:
                bounds[0]=0
            if bounds[2]<0:
                bounds[2]=0
            if bounds[1]>I.shape[1]:
                bounds[1]=I.shape[1]
            if bounds[3]>I.shape[2]:
                bounds[3]=I.shape[2]
        for i in np.arange(len(trace)):
            #trace[i]=I[i,x,y]
            trace[i]=np.mean(I[i,bounds[0]:bounds[1],bounds[2]:bounds[3]])
            #trace[i]=2*np.sum((I[i]-baseline)*I_norm)+baseline
            #I'm not sure why the 2 is needed, but this seems to always work out to 1:
            #            from analyze.puffs.gaussianFitting import gaussian
            #            x = np.arange(1000,dtype=float)
            #            y = np.arange(1000,dtype=float)
            #            amplitude=2.8
            #            sigma=7
            #            I=gaussian(x[:,None], y[None,:],100,100,sigma,amplitude)
            #            I_fit=gaussian(x[:,None], y[None,:],100,100,sigma,1)
            #            I_norm=I_fit/np.sum(I_fit)
            #            calc_amp=2*np.sum(I*I_norm)
            #            print('Original amp={}  Calculated amp={}'.format(amplitude,calc_amp))
        self.trace=trace
        self.kinetics['t_start']=t0
        self.kinetics['t_end']=t1
        self.calcRiseFallTimes()
        
    def calcRiseFallTimes(self):
        before=self.kinetics['before']
        t_start=self.kinetics['t_start']-before
        t_end=self.kinetics['t_end']-before
        baseline=self.trace[t_start]
        t_peak=np.argmax(self.trace[t_start:t_end+1])+t_start
        f_peak=self.trace[t_peak]
        if baseline>f_peak:
            baseline=self.trace[t_start]
        amplitude=f_peak-baseline
        thresh20=baseline+amplitude*.2
        thresh50=baseline+amplitude*.5
        thresh80=baseline+amplitude*.8
        tmp=np.argwhere(self.trace>thresh20); tmp=tmp[np.logical_and(tmp>=t_start,tmp<=t_peak)]; 
        if len(tmp)==0: r20=np.nan
        else:  r20=tmp[0]-t_start
        tmp=np.argwhere(self.trace>thresh50); tmp=tmp[np.logical_and(tmp>=t_start,tmp<=t_peak)];
        if len(tmp)==0: r50=np.nan
        else:  r50=tmp[0]-t_start
        tmp=np.argwhere(self.trace>thresh80); tmp=tmp[np.logical_and(tmp>=t_start,tmp<=t_peak)]; 
        if len(tmp)==0: r80=np.nan
        else:  r80=tmp[0]-t_start
        tmp=np.argwhere(self.trace<thresh80); tmp=tmp[tmp>=t_peak]; 
        if len(tmp)==0: f80=np.nan
        else: f80=tmp[0]-t_peak
        
        tmp=np.argwhere(self.trace<thresh50); tmp=tmp[tmp>=t_peak]; 
        if len(tmp)==0: f50=np.nan
        else: f50=tmp[0]-t_peak
        
        tmp=np.argwhere(self.trace<thresh20); tmp=tmp[tmp>=t_peak]; 
        if len(tmp)==0: f20=np.nan
        else: f20=tmp[0]-t_peak
        
        tmp=np.argwhere(self.trace<baseline); tmp=tmp[tmp>=t_peak]; 
        if len(tmp)==0: f0=np.nan
        else: 
            f0=tmp[0]
            if f0<t_end:
                t_end=f0
            f0=f0-t_peak
            
        self.kinetics['amplitude']=amplitude
        self.kinetics['baseline']=baseline
        #self.kinetics['t_end']=t_end+before
        self.kinetics['r20']=r20
        self.kinetics['r50']=r50
        self.kinetics['r80']=r80
        self.kinetics['f20']=f20
        self.kinetics['f50']=f50
        self.kinetics['f80']=f80
        self.kinetics['f0']=f0
        self.kinetics['t_peak']=t_peak+before

            
    def plot(self,figure=None):
        if figure is None:
            figure=pg.plot()
        k=self.kinetics
        baseline=k['baseline']; amplitude=k['amplitude']
        #thresh20=baseline+amplitude*.2
        #thresh50=baseline+amplitude*.5
        #thresh80=baseline+amplitude*.8
        x=np.arange(len(self.trace))+k['before']
        figure.plot(x,self.trace,pen=pg.mkPen(width=2))
        #figure.plot(x,self.fStrong,pen=pg.mkPen('g'))
        #figure.plot(x,self.fWeak,pen=pg.mkPen('r'))
        self.peakLine=figure.addLine(y=baseline,pen=pg.mkPen('y',style=Qt.DashLine))
        self.baselineLine=figure.addLine(y=baseline+amplitude,pen=pg.mkPen('y',style=Qt.DashLine))
        self.startLine=figure.addLine(x=k['t_start'],pen=pg.mkPen('y',style=Qt.DashLine),movable=True,bounds=(self.kinetics['before'],self.kinetics['t_peak']))
        self.endLine=figure.addLine(x=k['t_end'],pen=pg.mkPen('y',style=Qt.DashLine),movable=True, bounds=(self.kinetics['t_peak'],self.kinetics['after']))
        self.startLine.sigDragged.connect(self.changeStartTime)
        self.endLine.sigDragged.connect(self.changeEndTime)
    def changeStartTime(self,line):
        time=line.value()
        time=int(np.round(time))
        if time!=line.value():
            self.startLine.setValue(time)
        oldstart=self.kinetics['t_start']
        self.kinetics['t_start']=time
        if oldstart!=time:
            self.calcRiseFallTimes()
            self.baselineLine.setValue(self.kinetics['baseline'])
            self.peakLine.setValue(self.kinetics['baseline']+self.kinetics['amplitude'])
            self.endLine.setValue(self.kinetics['t_end'])
            self.puffs.puffAnalyzer.drawRedOverlay()
    def changeEndTime(self,line):
        time=line.value()
        time=int(np.round(time))
        if time!=line.value():
            self.endLine.setValue(time)
        oldend=self.kinetics['t_end']
        if oldend!=time:   
            self.kinetics['t_end']=time
            self.puffs.puffAnalyzer.drawRedOverlay()
    