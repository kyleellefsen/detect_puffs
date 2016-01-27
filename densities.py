# -*- coding: utf-8 -*-
"""
Created on Mon Jan 04 14:33:51 2016

@author: Kyle Ellefsen
"""
from process.progress_bar import ProgressBar
import numpy as np
import global_vars as g




def getMask(nt=5,nx=5,ny=5):
    mask=np.zeros((nt,nx,ny))
    center=np.array([(nt-1)/2, (nx-1)/2, (ny-1)/2]).astype(np.int)
    t0,x0,y0=center
    for t in np.arange(nt):
        for x in np.arange(nx):
            for y in np.arange(ny):
                if  ((t-t0)**2) / (t0**2) + ((x-x0)**2) / (x0**2)   +   ((y-y0)**2) / (y0**2) <= 1:
                    mask[t,x,y]=1
    return mask, center


###############################################################################
#                    Density
###############################################################################
def getDensities(Image,maxPuffLen, maxPuffDiameter):
    if maxPuffLen%2==0:
        maxPuffLen+=1
    if maxPuffDiameter%2==0:
        maxPuffDiameter+=1
    nCores = g.m.settings['nCores']
    pxls=np.array(np.where(Image)).T
    block_ends=np.linspace(0,len(pxls),nCores+1).astype(np.int)
    data_blocks=[pxls[block_ends[i]:block_ends[i+1]] for i in np.arange(nCores)]
    args=(Image,maxPuffLen, maxPuffDiameter)
    progress = ProgressBar(calcDensity, data_blocks, args, nCores, msg='Calculating Density')
    if progress.results is None or any(r is None for r in progress.results):
        result=None
        return result
    else:
        result=np.sum(progress.results,0)
    result=np.log(result+1)
    Densities=result
    return Densities
    
def calcDensity(q_results, q_progress, q_status, child_conn, args):
    pxls=child_conn.recv() # unfortunately this step takes a long time
    percent=0  # This is the variable we send back which displays our progress
    status=q_status.get(True) #this blocks the process from running until all processes are launched
    if status=='Stop':
        q_results.put(None) # if the user presses stop, return None
    
    Image, maxPuffLen, maxPuffDiameter=args #unpack all the variables inside the args tuple
    result=np.zeros(Image.shape)
    mt,mx,my=Image.shape
    mask, center=getMask(maxPuffLen, maxPuffDiameter, maxPuffDiameter)
    for i,pxl in enumerate(pxls):
        t,x,y=pxl
        try:
            result[t,x,y]=np.count_nonzero(mask*Image[t-center[0]:t+center[0]+1,x-center[1]:x+center[1]+1,y-center[2]:y+center[2]+1])
        except ValueError:
            t0=t-center[0]
            tf=t+center[0]+1
            x0=x-center[1]
            xf=x+center[1]+1
            y0=y-center[2]
            yf=y+center[2]+1
            mask2=mask
            if t0<0:
                mask2=mask2[center[0]-t:,:,:]
                t0=0
            if x0<0:
                mask2=mask2[:,center[1]-x:,:]
                x0=0
            if y0<0:
                mask2=mask2[:,:,center[2]-y:]
                y0=0
            if tf>mt-1:
                mask2=mask2[:-(tf-mt+1),:,:]
                tf=mt-1
            if xf>mx-1:
                mask2=mask2[:,:-(xf-mx+1),:]
                xf=mx-1
            if yf>my-1:
                mask2=mask2[:,:,:-(yf-my+1)]
                yf=my-1
            result[t,x,y]=np.count_nonzero(mask2*Image[t0:tf,x0:xf,y0:yf])
        
        if percent<int(100*i/len(pxls)):
            percent=int(100*i/len(pxls))
            q_progress.put(percent+2) #I have no idea why the last two percent aren't displayed, but I'm adding 2 so it reaches 100
        if not q_status.empty(): #check if the stop button has been pressed
            stop=q_status.get(False)
            q_results.put(None)
            return                 
    # finally, when we've finished with our calculation, we send back the result
    q_results.put(result)


    
    