# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 13:45:36 2014

@author: Kyle
I found another library which performs a nearly identical function after I 
wrote this.  It is located in https://github.com/ZhuangLab/storm-analysis/blob/master/sa_library/gaussfit.py
"""
import sys
from typing import Callable
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from qtpy import QtCore, QtWidgets
import scipy


def cosd(degrees:float) -> float:
    """Returns the cosine of the input angle in degrees."""
    return np.cos(np.radians(degrees))

def sind(degrees:float) -> float:
    """Returns the sine of the input angle in degrees."""
    return np.sin(np.radians(degrees))

def leastsqbound(func: Callable,
                 x0 : tuple,
                 args: tuple,
                 bounds: tuple,
                 ftol: float = 1.49012e-8,
                 full_output: bool = True) -> np.ndarray:
    """This function is a wrapper around scipy.optimize.least_squares which adds 
    bounds to the fitting. Originally I used 
    https://github.com/jjhelmus/leastsqbound-scipy, but this library broke with
    a numpy change. This uses a similar API to leastsqbound but calls
    scipy.optimize.least_squares.
    
    Returns:
      The fitted x0
    """
    result = scipy.optimize.least_squares(func, x0, args=args, bounds=bounds, ftol=ftol)
    return result.x


#### SYMETRICAL GAUSSIAN
def fitGaussian(image: np.ndarray,
                p0: tuple[float, float, float, float],
                bounds: tuple[tuple[float, float, float, float], tuple[float, float, float, float]],
                nGaussians: int=1,
                display: bool=False):
    """
    Takes an nxm matrix and returns an nxm matrix which is the gaussian fit
    of the first.  p0 is a list of parameters [xorigin, yorigin, sigma,amplitude]
    0-19 should be [-.2889 -.3265 -.3679 -.4263 -.5016 -.6006 ... -.0228 .01913]
    """

    axes = (np.arange(image.shape[0])[:, np.newaxis],
            np.arange(image.shape[1])[np.newaxis, :])
    if display:
        data = Puff3d(image,'Original Data')
    p0_rounded = tuple(round(p, 3) for p in p0)
    if nGaussians == 1:
        p_fit = leastsqbound(_err_single_symmetrical_gaussian, p0_rounded,
                             args=(image, axes), bounds=bounds, ftol=.0000001,
                             full_output=True)
        #xorigin,yorigin,sigmax,sigmay,angle,amplitude=p
        I_fit = gaussian(axes[0], axes[1], *p_fit)
        if not display:
            return p_fit, I_fit, I_fit
        else:
            fitted_data=Puff3d(I_fit,'Fitted')
            return data, fitted_data, p_fit
    elif nGaussians >= 2:
        p_fit = leastsqbound(_err_multiple_symmetric_gaussians, p0_rounded,
                             args=(image, axes), bounds=bounds, ftol=.0000001,
                             full_output=True)
        #xorigin,yorigin,sigmax,sigmay,angle,amplitude=p
        I_fit = gaussian(axes[0], axes[1], *(p_fit[:4]))
        I_fit2 = gaussian2(axes[0], axes[1], *p_fit)
        return p_fit[:4], I_fit, I_fit2
        
def gaussian(x: np.ndarray,
             y: np.ndarray,
             xorigin: float,
             yorigin: float,
             sigma: float,
             amplitude: float):
    """Generates a 2D Gaussian."""
    return amplitude*(np.exp(-(x-xorigin)**2/(2.*sigma**2))*np.exp(-(y-yorigin)**2/(2.*sigma**2)))

def gaussian_1var(p: tuple[float, float, float, float],
                  axes: tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    '''xorigin,yorigin,sigmax,sigmay,angle'''
    xorigin, yorigin, sigma, amplitude = p
    x0 = axes[0]
    x1 = axes[1]
    return amplitude*(np.exp(-(x0-xorigin)**2/(2.*sigma**2))*np.exp(-(x1-yorigin)**2/(2.*sigma**2)))

def _err_single_symmetrical_gaussian(p: tuple[float, float, float, float],
        y: np.ndarray,
        axes: tuple[np.ndarray, np.ndarray]):
    ''' 
    p is a tuple contatining the parameters.  p = (xorigin, yorigin, sigma, amplitude)
    y is the data we are fitting to (the dependent variable)
    axes are the image axes indices.  0 is the x axis, 1 is the y axis
    '''
    remander = y - gaussian_1var(p, axes)
    return remander.ravel()
    
def gaussian2(x: np.ndarray, y: np.ndarray, *args):
    '''xorigin,yorigin,sigma,angle'''
    answer = np.zeros((x.shape[0], y.shape[1]))
    for i in np.arange(0, len(args), 4):
        xorigin, yorigin, sigma, amplitude = args[i:i+4]
        answer += amplitude*(np.exp(-(x-xorigin)**2/(2.*sigma**2))*np.exp(-(y-yorigin)**2/(2.*sigma**2)))
    return answer

def gaussian_1var2(p: tuple[float, float, float, float],
                   axes: tuple[np.ndarray, np.ndarray]):
    '''xorigin, yorigin, sigma, angle'''
    answer = np.zeros((axes[0].shape[0], axes[1].shape[1]))
    for i in np.arange(0, len(p), 4):
        xorigin, yorigin, sigma, amplitude = p[i:i+4]
        x0 = axes[0]-xorigin
        x1 = axes[1]-yorigin
        answer += amplitude * (
            np.exp(-x0**2/(2.*sigma**2)) * 
            np.exp(-x1**2/(2.*sigma**2))
            )
    return answer

def _err_multiple_symmetric_gaussians(p: tuple[float, float, float, float],
         y: np.ndarray,
         axes: tuple[np.ndarray, np.ndarray]):
    """
    p is a tuple contatining the initial parameters.  
      p = (xorigin, yorigin, sigma, amplitude)
    y is the data we are fitting to (the dependent variable)
    axes are the image axes indices.  0 is the x axis, 1 is the y axis
    """
    remander = y - gaussian_1var2(p, axes)
    return remander.ravel()



#### ROTATABLE GAUSSIAN
def fitRotGaussian(image: np.ndarray,
                   p0: tuple,
                   bounds: tuple,
                   nGaussians: int=1, 
                   display: bool=False):
    '''
    Takes an nxm matrix and returns an nxm matrix which is the gaussian fit
     p0 is a list of parameters [xorigin, yorigin, sigmax, sigmay, amplitude]
    '''
    axes = (np.arange(image.shape[0])[:, np.newaxis],
            np.arange(image.shape[1])[np.newaxis, :])
    if display:
        data = Puff3d(image,'Original Data')
    if nGaussians==1:
        p_fit = leastsqbound(err_rot, p0, args=(image, axes), bounds=bounds,ftol=.0000001,full_output=True)
        #xorigin,yorigin,sigmax,sigmay,angle,amplitude=p
        I_fit = gaussian_rot(axes[0], axes[1], *p_fit)
        if not display:
            return p_fit, I_fit, I_fit
        else:
            fitted_data = Puff3d(I_fit,'Fitted')
            return data, fitted_data, p_fit
    elif nGaussians>=2:
        p_fit = leastsqbound(err_rot2, p0,args=(image, axes), bounds = bounds,ftol=.0000001,full_output=True)
        #xorigin,yorigin,sigmax,sigmay,angle,amplitude=p
        p_short = p_fit[:6]
        I_fit = gaussian_rot(axes[0], axes[1], *p_short)
        I_fit2 = gaussian_rot2(axes[0], axes[1], *p_fit)
        return p_fit[:6], I_fit, I_fit2

def gaussian_rot(x, y, xorigin, yorigin, sigmax, sigmay, angle, amplitude):
    '''xorigin,yorigin,sigmax,sigmay,angle'''
    x = x-xorigin
    y = y-yorigin
    x2 = cosd(angle)*x-sind(angle)*y
    y2 = sind(angle)*x+cosd(angle)*y
    return amplitude*(np.exp(-x2**2/(2.*sigmax**2))*np.exp(-y2**2/(2.*sigmay**2)))

def gaussian_rot_1var(p, axes): #INPUT_MAT,xorigin,yorigin,sigma):
    '''xorigin,yorigin,sigmax,sigmay,angle, amplitude'''
    xorigin,yorigin,sigmax,sigmay,angle,amplitude = p
    x = axes[0]-xorigin
    y = axes[1]-yorigin
    x2 = cosd(angle)*x-sind(angle)*y
    y2 = sind(angle)*x+cosd(angle)*y
    return amplitude*(np.exp(-x2**2/(2.*sigmax**2))*np.exp(-y2**2/(2.*sigmay**2)))

def err_rot(p: tuple, y: np.ndarray, axes: tuple[np.ndarray, np.ndarray]):
    remander = y - gaussian_rot_1var(p, axes)
    return remander.ravel()
    
def gaussian_rot2(x, y, *args):
    '''xorigin,yorigin,sigmax,sigmay,angle'''
    answer = np.zeros((x.shape[0], y.shape[1]))
    for i in np.arange(0,len(args),6):
        xorigin,yorigin,sigmax,sigmay,angle,amplitude=args[i:i+6]
        x2 = x-xorigin
        y2 = y-yorigin
        x3 = cosd(angle)*x2-sind(angle)*y2
        y3 = sind(angle)*x2+cosd(angle)*y2
        answer += amplitude*(np.exp(-x3**2/(2.*sigmax**2))*np.exp(-y3**2/(2.*sigmay**2)))
    return answer

def gaussian_rot_1var2(p, axes:tuple[np.ndarray, np.ndarray]): #INPUT_MAT,xorigin,yorigin,sigma):
    '''xorigin,yorigin,sigmax,sigmay,angle, amplitude'''
    answer=np.zeros((axes[0].shape[0], axes[1].shape[1]))
    for i in np.arange(0, len(p), 6):
        xorigin,yorigin,sigmax,sigmay,angle,amplitude = p[i:i+6]
        x = axes[0]-xorigin
        y = axes[1]-yorigin
        x2 = cosd(angle)*x-sind(angle)*y
        y2 = sind(angle)*x+cosd(angle)*y
        answer += amplitude*(np.exp(-x2**2/(2.*sigmax**2))*np.exp(-y2**2/(2.*sigmay**2)))
    return answer

def err_rot2(p: tuple, y: np.ndarray, axes: tuple[np.ndarray, np.ndarray]):
    remander = y - gaussian_rot_1var2(p, axes)
    return remander.ravel()




class Puff3d(gl.GLViewWidget):
    def __init__(self,image,title,parent=None):
        super(Puff3d,self).__init__(parent)
        self.setCameraPosition(distance=50)
        self.grid= gl.GLGridItem()
        #self.grid.scale(2,2,1)
        self.grid.setDepthValue(10) # draw grid after surfaces since they may be translucent
        self.addItem(self.grid)
        #z = ndimage.gaussian_filter(np.random.normal(size=(50,50)), (1,1))
        self.p1 = gl.GLSurfacePlotItem(z=image, shader='shaded', color=(0.5, 0.5, 1, 1))
        self.p1.scale(1, 1, 10.0)
        #self.p1.translate(-volume.shape[1]/2, -volume.shape[2]/2, -volume.shape[0]/4)
        self.addItem(self.p1)
        #self.setMinimumHeight(300)
        self.setWindowTitle(title)
        self.show()

    def updateimage(self,image):
        self.p1.setData(z=image)
        
        
        
def main():
    app = QtWidgets.QApplication([])
    image_size = 30
    x = np.arange(image_size, dtype=float)
    y = np.arange(image_size, dtype=float) 
    p_true = np.array([5.1, 18.5777, 2.2, 3])
    image = gaussian(x[:,None], y[None,:],
                                     xorigin=p_true[0],
                                     yorigin=p_true[1],
                                     sigma=p_true[2],
                                     amplitude=p_true[3])
    rng = np.random.RandomState(30)
    noise = 0.5 * rng.rand(image_size, image_size)
    p0 = (5, 5, 3, 1)
    bounds = [(0, 0, 1, 0), (30, 30, 10, 5)]
    data, fitted_data, p = fitGaussian(image+noise, p0, bounds, display=True)
    print(f"Fitted Values: {p}")
    
    sys.exit(app.exec_())

if __name__=='__main__':
    main()
