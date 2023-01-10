

import numpy as np
from os import walk
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



class Utils(object):


    @staticmethod
    def loadFolder(folder):
        #input: folder path including files (of anchor tracks in 3D)
        #output data_full: matrix with [file_number, track_length, [frame,y,y,z,state_label{0:free,1:extruion,2:closed}]] ; state_label only if available
        #output data_full_fileindex: file_index )
        data_full=[]
        data_full_fileindex=[]
        files = []
        molecule=[]
        for (dirpath, dirnames, filenames) in walk(folder):
            #f.extend(filenames)
            for ff in sorted(filenames):
                if (ff.endswith('.txt')):
                    files.append(dirpath+'/'+ff)

        for fid,f in enumerate(files):

            index1=f.rindex('_')
            index2=f.rindex('.txt')

            if fid>=len(molecule):
                molecule.append([])
            data=(np.loadtxt(f, delimiter=' ', skiprows=0))
            #plt.figure(figsize=(5, 4), dpi=80)
            frame=np.transpose(np.arange(np.size(data,0))[np.newaxis])#make frame column
            data=np.concatenate((frame,data),axis=1)#add frame number to data
            data_full.append(data)
            data_full_fileindex.append(int(f[index1+1:index2]))
            #if fid>=20:
            #    print("ERROR: only 1 file loaded")
            #    break
        return data_full,data_full_fileindex

    @staticmethod
    def exactFrames(data,frame_start=0,frame_end=sys.maxsize):
        #input: 4 column matrix of coordinates (frame,x,y,z)
        #output: x,y,z columns between frame_start and frame_end index
        cutted=[]
        for d in data:
            if frame_end>np.size(d,0):
                fe=np.size(d,0)
            else:
                fe=frame_end
            if frame_start<0:
                fs=0
            else:
                fs=frame_start

            cutted.append(d[fs:fe,1:4])

        return cutted


    #
    @staticmethod
    def exactFramesByLabel(data,label):
        #input: 5 column matrix of coordinates (frame,x,y,z,label)
        #output: x,y,z columns whose index (4th column) equals input label
        cutted=[]
        for d in data:
            cutted.append(d[d[:,4]==label,1:4])

        return cutted




    @staticmethod
    def transform_data(distribution0input,distribution1input,shift=0,randomAngle=False):
        #input: 2 anchors coordinates
        #output: difference of 2 anchors
        distribution0=np.vstack(distribution0input)
        distribution1=np.vstack(distribution1input)
        data_tmp=np.subtract(distribution0,distribution1)
        #free_data0=np.subtract(distribution0,distribution0)#center data0
        #relative position of data1
        if randomAngle:

            #get spherical coordinate to change norm and angle randomly, and go back to cartesian coordinate

            #this will be wrong if loc_precision(x),loc_precision(y),loc_precision(z) are different
            #I need to do it because mean != 0, meaning that chromosome is mostly oriented
            p=np.sqrt(np.sum(np.multiply(data_tmp,data_tmp),axis=1))

            #p=np.subtract(p,min(p))

            #theta=np.arccos(np.divide(data_tmp[:,2],p))
            #phi=np.arctan2(data_tmp[:,1],data_tmp[:,0])
            #put random angle on the surface of a sphere:
            new_theta=np.arccos(np.subtract(np.multiply(np.random.rand(np.size(p)),2),1))
            new_phi=np.multiply(np.random.rand(np.size(p)),2.*3.141592)
            if shift>0:
                print('WARNING: here, I removed ',shift,' to the norm')
            p=np.abs(p-shift)

            data_x=np.multiply(p,np.multiply(np.sin(new_theta),np.cos(new_phi)))
            data_y=np.multiply(p,np.multiply(np.sin(new_theta),np.sin(new_phi)))
            data_z=np.multiply(p,np.cos(new_theta))
            data_x=np.transpose(data_x[np.newaxis])
            data_y=np.transpose(data_y[np.newaxis])
            data_z=np.transpose(data_z[np.newaxis])
            data_new=np.concatenate((data_x,data_y,data_z),axis=1)
            #print(np.shape(free_data_randomAngle))
            #angle=atan2(y,x)+pi;



        else:

            #get spherical coordinate to change norm, and go back to cartesian coordinate
            if shift>0:
                p=np.sqrt(np.sum(np.multiply(data_tmp,data_tmp),axis=1))

                devp=np.divide(data_tmp[:,2],p)

                theta=np.arccos(np.divide(data_tmp[:,2],p))


                #print(data_tmp[np.isnan(theta),2],p[np.isnan(theta)])

                phi=np.arctan2(data_tmp[:,1],data_tmp[:,0])

                p=np.abs(p-shift)

                data_x=np.multiply(p,np.multiply(np.sin(theta),np.cos(phi)))
                data_y=np.multiply(p,np.multiply(np.sin(theta),np.sin(phi)))
                data_z=np.multiply(p,np.cos(theta))
                data_x=np.transpose(data_x[np.newaxis])
                data_y=np.transpose(data_y[np.newaxis])
                data_z=np.transpose(data_z[np.newaxis])
                data_new=np.concatenate((data_x,data_y,data_z),axis=1)
            else:
                data_new=data_tmp






        return(data_new)







    @staticmethod
    def transform_data_shiftEndOfExtrusion(distribution0input,distribution1input,shift,startShiftFromEnd,randomAngle=False):
        #input: 2 anchors coordinates
        #startShiftFromEnd=60 indicates we start to shift 60 time points from the end
        #output: difference of 2 anchors

        distribution0=np.array(distribution0input,dtype=object)
        distribution1=np.array(distribution1input,dtype=object)
        data_tmp=np.subtract(distribution0,distribution1)

        #data_tmp=np.subtract(distribution0,distribution1)

        #free_data0=np.subtract(distribution0,distribution0)#center data0
        #relative position of data1
        if randomAngle:

            data_new=[]
            for i,d in enumerate(data_tmp):
                #get spherical coordinate to change norm and angle randomly, and go back to cartesian coordinate

                #this will be wrong if loc_precision(x),loc_precision(y),loc_precision(z) are different
                #I need to do it because mean != 0, meaning that chromosome is mostly oriented
                p=np.sqrt(np.sum(np.multiply(d,d),axis=1))

                start=len(d)-startShiftFromEnd
                end=len(d)
                smoothShift=np.arange(startShiftFromEnd)
                smoothShift=smoothShift*shift/max(smoothShift)
                #print(start,end,startShiftFromEnd,len(smoothShift))
                p[start:end]-=smoothShift

                #if shift>0:
                #    print('WARNING: here, I removed ',shift,' to the norm')
                #p=np.subtract(p,min(p))

                #theta=np.arccos(np.divide(data_tmp[:,2],p))
                #phi=np.arctan2(data_tmp[:,1],data_tmp[:,0])
                #put random angle on the surface of a sphere:
                new_theta=np.arccos(np.subtract(np.multiply(np.random.rand(np.size(p)),2),1))
                new_phi=np.multiply(np.random.rand(np.size(p)),2.*3.141592)

                data_x=np.multiply(p,np.multiply(np.sin(new_theta),np.cos(new_phi)))
                data_y=np.multiply(p,np.multiply(np.sin(new_theta),np.sin(new_phi)))
                data_z=np.multiply(p,np.cos(new_theta))
                data_x=np.transpose(data_x[np.newaxis])
                data_y=np.transpose(data_y[np.newaxis])
                data_z=np.transpose(data_z[np.newaxis])
                data_new.append(np.concatenate((data_x,data_y,data_z),axis=1))
                #print(np.shape(free_data_randomAngle))
                #angle=atan2(y,x)+pi;



        else:

            #get spherical coordinate to change norm, and go back to cartesian coordinate
            if shift>0:
                data_new=[]
                for i,d in enumerate(data_tmp):

                    p=np.abs(np.sqrt(np.sum(np.multiply(d,d),axis=1)))



                    #print('WARNING: here, I removed ',shift,' to the norm')

                    theta=np.arccos(np.divide(d[:,2],p))




                    phi=np.arctan2(d[:,1],d[:,0])
                    #for each simulation

                    start=len(d)-startShiftFromEnd
                    end=len(d)
                    smoothShift=np.arange(startShiftFromEnd)
                    smoothShift=smoothShift*shift/max(smoothShift)
                    p[start:end]-=smoothShift



                    data_x=np.multiply(p,np.multiply(np.sin(theta),np.cos(phi)))
                    data_y=np.multiply(p,np.multiply(np.sin(theta),np.sin(phi)))
                    data_z=np.multiply(p,np.cos(theta))
                    data_x=np.transpose(data_x[np.newaxis])
                    data_y=np.transpose(data_y[np.newaxis])
                    data_z=np.transpose(data_z[np.newaxis])
                    data_new.append(np.concatenate((data_x,data_y,data_z),axis=1))
            else:
                data_new=data_tmp



        data_new=np.vstack(data_new)


        return(data_new)

    @staticmethod
    def randDataN(data,sigmaX,sigmaY,sigmaZ):
        #input: difference of 3D coordinates
        #output: difference of 3D coordinates with additive random gaussian noise
        datarand=np.zeros(np.shape(data))
        datarand[:,0]=np.add(data[:,0],np.random.normal(loc=0.0, scale=sigmaX, size=len(data[:,0])))
        datarand[:,1]=np.add(data[:,1],np.random.normal(loc=0.0, scale=sigmaY, size=len(data[:,1])))
        datarand[:,2]=np.add(data[:,2],np.random.normal(loc=0.0, scale=sigmaZ, size=len(data[:,2])))
        return(datarand)


    @staticmethod
    def dataShuffleWithProportions(free_data,extr_data,loop_data,p_free,p_extr,p_loop,N):

        full_data_shuffle=None

        extr_tmp=extr_data.copy()
        loop_tmp=loop_data.copy()
        free_tmp=free_data.copy()

        np.random.shuffle(extr_tmp)
        np.random.shuffle(loop_tmp)
        np.random.shuffle(free_tmp)


        n_extr=N*p_extr
        n_free=N*p_free
        n_loop=N*p_loop

        n=0
        i=0
        while n<n_extr:
            numberMax=int(min(np.shape(extr_tmp)[0],n_extr-n))

            datawithindex=np.insert(extr_tmp[0:numberMax,:], 3, np.zeros(numberMax)+1, axis=1)
            if full_data_shuffle is None:
                full_data_shuffle=datawithindex
            else:
                full_data_shuffle=np.concatenate((full_data_shuffle,datawithindex),axis=0)#add frame number to data
            n+=numberMax
            i+=1

        n=0
        i=0
        while n<n_loop:
            numberMax=int(min(np.shape(loop_tmp)[0],n_loop-n))
            datawithindex=np.insert(loop_tmp[0:numberMax,:], 3, np.zeros(numberMax)+0, axis=1)
            if full_data_shuffle is None:
                full_data_shuffle=datawithindex
            else:
                full_data_shuffle=np.concatenate((full_data_shuffle,datawithindex),axis=0)#add frame number to data
            n+=numberMax
            i+=1


        n=0
        i=0
        while n<n_free:
            numberMax=int(min(np.shape(free_tmp)[0],n_free-n))
            datawithindex=np.insert(free_tmp[0:numberMax,:], 3, np.zeros(numberMax)+2, axis=1)
            if full_data_shuffle is None:
                full_data_shuffle=datawithindex
            else:
                full_data_shuffle=np.concatenate((full_data_shuffle,datawithindex),axis=0)#add frame number to data
            n+=numberMax
            i+=1

        return full_data_shuffle


    @staticmethod
    def makeHistogram(data,bins):
        #convert 3D distribution to 3 1D histograms concatenated
        tmp=[]
        for dim in [0,1,2]:
            tmp.append(np.histogram(data[:,dim], bins = (bins))[0])
        f=np.asarray(tmp).ravel()
        return(f)


    @staticmethod
    def getMergedXAxis(xyz):
        lenvect=int(len(xyz)//3)
        return(xyz[0:lenvect])




    @staticmethod
    def mergecurve(curve,xyz):
        xaxis=Utils.getMergedXAxis(xyz)
        lenvect=int(len(xyz)//3)
        s=np.sum(curve)
        res=np.zeros(lenvect)
        res+=curve[0:lenvect:]

        res+=curve[lenvect:lenvect+lenvect:]

        res+=curve[2*lenvect:2*lenvect+lenvect:]

        return(res*s/np.sum(res))





















    def toto():
        print('toto')
