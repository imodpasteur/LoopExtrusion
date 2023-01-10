
import numpy as np
from os import walk
import os
import matplotlib.pyplot as plt
import sys
import random
from scipy.optimize import curve_fit
from numpy import sqrt, pi, exp, linspace
from time import process_time
from scipy import integrate
from scipy import ndimage
import pandas as pd
import scipy.stats as stats

from utils import Utils
from models import Models


#folder of first anchor without loop
folder_free0='../Anchor-anchor_simulated_data_subsample/Free/Bead275'
#folder of second anchor without loop
folder_free1='../Anchor-anchor_simulated_data_subsample/Free/Bead324'

#folder of first anchor with loop
folder_loop0='../Anchor-anchor_simulated_data_subsample/Loop/Bead275'
#folder of second anchor with loop
folder_loop1='../Anchor-anchor_simulated_data_subsample/Loop/Bead324'


show=True

#simulated localization errors:
sigma_X=[0.015]
sigma_Y=[0.015]
sigma_Z=[0.03]

#simulated proportions of closed, extrusion and free states:
P_loop=[0.35]
P_extr=[0.15]
P_free=[0.5]


N=10000 #number of anchor-anchor differences simulated










#load list of files of anchor coordinates
free_poly1,free_poly1_indexfile=Utils.loadFolder(folder_free1)
free_poly0,free_poly0_indexfile=Utils.loadFolder(folder_free0)

loop_poly0,loop_poly0_indexfile=Utils.loadFolder(folder_loop0)
loop_poly1,loop_poly1_indexfile=Utils.loadFolder(folder_loop1)





#extract coordinates of anchors according to states {Free, extrusion, closed, all mixed}

free_data0=Utils.exactFrames(free_poly0)
free_data1=Utils.exactFrames(free_poly1)


loop_data0=Utils.exactFramesByLabel(loop_poly0,label=2)
loop_data1=Utils.exactFramesByLabel(loop_poly1,label=2)

extr_data0=Utils.exactFramesByLabel(loop_poly0,label=1)
extr_data1=Utils.exactFramesByLabel(loop_poly1,label=1)




#compute differences of coordinates from 2 anchors coordinates
free_data=Utils.transform_data(free_data0,free_data1,randomAngle=False)
loop_data=Utils.transform_data(loop_data0,loop_data1,shift=.05,randomAngle=False)#shift about 0.05µm to compensate bead size in simulations
extr_data=Utils.transform_data_shiftEndOfExtrusion(extr_data0,extr_data1,.05,60,randomAngle=False)#shift the end of extrusion (60 last time points) to compensate bead size in simulations



for sigma_x,sigma_y,sigma_z in zip(sigma_X,sigma_Y,sigma_Z):

    print('----------------------------------------------------------')
    print('simulated localization precision :  s(X)=',sigma_x,'  s(Y)=',sigma_y,'  s(Z)=',sigma_z,' nm')
    print('----------------------------------------------------------')

    #Be careful, here we multiply var by 2 because anchor to anchor distances contains 2 times localization errors
    var_x=2*sigma_x**2
    var_y=2*sigma_y**2
    var_z=2*sigma_z**2





    loop_data_noisy=Utils.randDataN(loop_data,np.sqrt(var_x),np.sqrt(var_y),np.sqrt(var_z))
    free_data_noisy=Utils.randDataN(free_data,np.sqrt(var_x),np.sqrt(var_y),np.sqrt(var_z))
    extr_data_noisy=Utils.randDataN(extr_data,np.sqrt(var_x),np.sqrt(var_y),np.sqrt(var_z))


    models = Models(var_x,var_y,var_z)



    maxibin=1 #1 µm width
    binstep=0.005 #5 nm step


    bins=np.arange(-maxibin,maxibin+binstep,binstep)
    binsub=bins[0:-1]+binstep/2
    xyz=np.array([binsub]*3).ravel()
    lenvect=len(xyz)//3




    #the last column is class (loop=0;extr=1;free=2)



    free_data_hist=Utils.makeHistogram(free_data_noisy,bins)
    extr_data_hist=Utils.makeHistogram(extr_data_noisy,bins)
    loop_data_hist=Utils.makeHistogram(loop_data_noisy,bins)

    free_data_hist=np.divide(free_data_hist,np.sum(free_data_hist))
    extr_data_hist=np.divide(extr_data_hist,np.sum(extr_data_hist))
    loop_data_hist=np.divide(loop_data_hist,np.sum(loop_data_hist))








    #estimate variance of Free state:
    init_vals = [max(free_data_hist), .2]     # for [amp, cen, wid]
    best_vals, covar = curve_fit(models.gaussian, xyz, free_data_hist, p0=init_vals,bounds=(0, np.inf))
    ampl_free=best_vals[0]
    var_free=abs(best_vals[1])


    #estimate variance of Closed state:
    init_vals = [np.sum(loop_data_hist)/np.sum(models.gaussian(xyz, 1, 0.0001)),0.0001]
    best_vals, covar = curve_fit(models.gaussian, xyz, loop_data_hist, p0=init_vals, maxfev=10000,bounds=(0, np.inf))
    ampl_loop=best_vals[0]
    var_loop=max(best_vals[1],0)







    #precompute Free model
    precomp_free_data = models.gaussian(xyz, 1, var_free)
    precomp_free_data/=np.sum(precomp_free_data)
    def precomp_free(x,amplitude):
        return np.multiply(precomp_free_data,amplitude)

    #precompute loop model
    precomp_loop_data = models.gaussian(xyz, 1, var_loop)
    precomp_loop_data/=np.sum(precomp_loop_data)
    def precomp_loop(x,amplitude):
        return np.multiply(precomp_loop_data,amplitude)


    #precompute extrusion model
    precomp_extr_data = models.gaussian_integrated_fast(np.array(xyz), 1, var_loop, var_free,True)
    precomp_extr_data/=np.sum(precomp_extr_data)
    def precomp_extr(x,amplitude):
        return np.multiply(precomp_extr_data,amplitude)





    #precompute Mixture (Free+Extrusion+ClosedLoop) model
    def precomp_full(x,  ampl_loop,  ampl_free,  ampl_extr):#gaussian model 3d WITHOUT COVARIANCE !!!
        return (np.multiply(precomp_loop_data,ampl_loop)+np.multiply(precomp_extr_data,ampl_extr)+np.multiply(precomp_free_data,ampl_free))




    print('----------------------------------------------------------')
    print('sigma closed estimated = ',np.sqrt(var_loop),' nm')
    print('sigma free estimated = ',np.sqrt(var_free),' nm')
    print('----------------------------------------------------------')


    if show:
        #show distribution of anchor-anchor coordinate differences according to 3 states after fitting them with theoretical models independently
        plt.figure(figsize=(8, 8), dpi=80)
        modelFitted_loop = precomp_loop(xyz, 1)
        plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(loop_data_hist,xyz),'#76FF7B', label='simulation closed')
        plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(modelFitted_loop,xyz),'g--', label='model closed')

        modelFitted_extr = precomp_extr(xyz, 1)
        plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(extr_data_hist,xyz),'#7BC8F6', label='simulation extrusion')
        plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(modelFitted_extr,xyz),'b--', label='model extrusion')

        modelFitted_free = precomp_free(xyz, 1)
        plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(free_data_hist,xyz),'#FFC0CB', label='simulation free')
        plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(modelFitted_free,xyz),'r--', label='model free')

        plt.legend()
        plt.title('independant Gaussian fitting of closed and free distributions ($\sigma_{closed}$ & $\sigma_{free}$  estimation)')
        plt.show()












    #proportion estimation of free, extrusion and closed states:
    for p_loop,p_extr,p_free in zip(P_loop,P_extr,P_free):
        print('----------------------------------------------------------')
        print('Closed proportion simulated : ',p_loop)
        print('Extrusion proportion simulated : ',p_extr)
        print('Free proportion simulated : ',p_free)
        print('----------------------------------------------------------')





        full_data_shuffle=Utils.dataShuffleWithProportions(free_data_noisy,extr_data_noisy,loop_data_noisy,p_free,p_extr,p_loop,N)#in case of loops, we can increase N here
        full_data_shuffle2=full_data_shuffle.copy()
        np.random.shuffle(full_data_shuffle2)



        simuldata=full_data_shuffle2[0:N,:]

        simuldata_full=Utils.makeHistogram(simuldata,bins)
        sumNorm=np.sum(simuldata_full)
        simuldata_full=np.divide(simuldata_full,sumNorm)


        #estimate split closed, extrusion and free states proportions
        init_vals = [.33,.33,.33]     # for [amp, cen, wid]

        best_vals, covar = curve_fit(precomp_full,xyz, simuldata_full,p0=init_vals,bounds=(0, np.inf))

        ampl_loop_fitted  =best_vals[0]
        ampl_free_fitted  =best_vals[1]
        ampl_extr_fitted  =best_vals[2]



        print('----------------------------------------------------------')
        print('Closed proportion estimated : ',ampl_loop_fitted)
        print('Extrusion proportion estimated : ',ampl_extr_fitted)
        print('Free proportion estimated : ',ampl_free_fitted)
        print('----------------------------------------------------------')



        if show:

            model_full=precomp_full(xyz, ampl_loop_fitted, ampl_free_fitted, ampl_extr_fitted)

            #split closed, extrusion and free states from simulated anchor-anchor differences
            simuldata_loop=simuldata[simuldata[:,3]==0]
            simuldata_extr=simuldata[simuldata[:,3]==1]
            simuldata_free=simuldata[simuldata[:,3]==2]

            simuldata_free=Utils.makeHistogram(simuldata_free,bins)
            simuldata_loop=Utils.makeHistogram(simuldata_loop,bins)
            simuldata_extr=Utils.makeHistogram(simuldata_extr,bins)

            simuldata_free=np.divide(simuldata_free,sumNorm)
            simuldata_loop=np.divide(simuldata_loop,sumNorm)
            simuldata_extr=np.divide(simuldata_extr,sumNorm)


            #compute models for closed, extrusion and free states
            model_extr=precomp_extr(xyz, ampl_extr_fitted)
            model_loop=precomp_loop(xyz, ampl_loop_fitted)
            model_free=precomp_free(xyz, ampl_free_fitted)



            plt.figure(figsize=(12, 12), dpi=80)

            plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(simuldata_full,xyz),'k*', label='simulation')
            plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(model_full,xyz),'k--', label='model')


            plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(simuldata_loop,xyz),'g*', label='simulation closed')
            plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(model_loop,xyz),'g--', label='model closed')



            plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(simuldata_free,xyz),'b*', label='simulation free')
            plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(model_free,xyz),'b--', label='model free')




            plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(simuldata_extr,xyz),'r*', label='simulation extr')
            plt.plot(Utils.getMergedXAxis(xyz),Utils.mergecurve(model_extr,xyz),'r--', label='model extr')


            plt.legend()
            plt.title('independant Gaussian fitting of closed and free distribution ($\sigma_{closed}$ & $\sigma_{free}$  estimation)')
            plt.show()











print('proportion estimations finished')
