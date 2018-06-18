import numpy as np
import tkinter as tk
from tkinter import filedialog
import glob
from tqdm import tqdm
import math
import scipy as sp
import scipy.signal as sps
import pandas as pd

def tt_extractor(filterFlg,upSampleFlg,figFlg):

    """Prompts the user to select a file directory; this directory is used to find *.txt files to process"""
    """REQUIRES THAT ENTIRE FILE HEADER IS COMMENTED OUT"""
    root = tk.Tk() #starts tkinter, a GUI package, as root
    root.withdraw() #retrieves the root window
    file_path = filedialog.askdirectory() #asks the user for folder selection

    """If no file path is selected, the program ends"""
    if file_path == 0: #if the file_path does not exist
        ok_flag = 0 #set ok_flag = 0
        return ok_flag #return ok_flag, ending the function and the program

    """The following initializes various constants used in the program based on experimental condiditons and theoretical constants"""
    fs = 20000 #audio sampling frequency, Hz
    dts = 1000/fs #converting the audio sampling frequency from Hz to ms; time of one sample
    T_record = 500 #time of a single record, ms
    dtt = T_record/1000 #converting time of a single record from ms to s
    N_record = round(T_record/dts) #converting length of a single record to samples
    file_path = file_path + '/*.txt' #adds .txt to the end of the selected file path, specifying what type of file the program should search for
    F = (glob.glob(file_path)) #retrieves all .txt files in the selected directory and saves the names as strings in a list
    ok_flag = np.zeros((len(F),1)) #ok_flag is now a vector of length F filled with zeroes. if two files are found, [0 0]
    T_signal = 5.8 #time of transmission of each chirp in ms
    N_signal = round(T_signal/dts) #length of the transmitted signal in samples
    M = 120 #number of transmissions/receptions in a file
    Fc = 1.2 #kHz, central frequency of the chirp signal
    bandWidth = 0.7 #signal's half bandwidth, kHz

    # filtration of whole record
    df = 1/N_record/dts # kHz, spectral resolution in the frequency domain
    ztn = (range(N_record//2))
    ztn = [x * df for x in ztn]
    ztn = [round(x,3) for x in ztn]
    f = [] # frequencies for signal filtration
    f.append(ztn.index(Fc-bandWidth))
    f.append(ztn.index(Fc+bandWidth))
    f.append(N_record-f[1])
    f.append(N_record-f[0])

    S=8 #number of speakers
    R=8 #number of microphones
    # expected speaker signal locations (in samples) on time axis
    t1Exp = [] #initialize the expected speaker signal location array. each of the folowing append commands add to the end of the (currently empty) array signifying the expected place on the time axis, in samples, of speaker output 
    t1Exp.append(4964)
    t1Exp.append(4164)
    t1Exp.append(8164)
    t1Exp.append(4)
    t1Exp.append(6404)
    t1Exp.append(8004)
    t1Exp.append(1604)
    t1Exp.append(5764)
    t1Exp = [x - 3 for x in t1Exp] #subtracts 3 element-wise from the entire array
    t1Exp = [x - 1 for x in t1Exp] #subtracts 1 element-wise from the entire array
    t1Exp = [x * (fs/40000) for x in t1Exp] #multiplies the entire array element-wise by (fs/40000)
    t1Exp = [round(x) for x in t1Exp] #rounds the entire array element-wise
    t1Exp = [x + 1 for x in t1Exp] #adds 1 to the entire array element-wise

    #upsampling
    N_signal1 = N_signal #unnecessary copying of N_signal
    if upSampleFlg: #if, when initializing the tt_extractor function, the upSampleFlg input existed
        N_signal1 = upSampleFlg * N_signal
        dts = dts/upSampleFlg
        t1Exp = [x - 1 for x in t1Exp]
        t1Exp = [x * upSampleFlg for x in t1Exp]
        t1Exp = [x + 1 for x in t1Exp]

    searchLag = round(7/dts)
    searchLag_S = round(3/dts)

    print(F)

    for fi in range(0, len(F)): #for each .txt file in F, repeat the following

        file_name = F[fi] #set file_name = the string in F representing the filename
        a = np.genfromtxt(file_name) #make a an array = the contents of the file

        tt = np.zeros((64,120))
        s_r = np.zeros((3480,120,8,8))

        with tqdm(total = 7680) as pbar: #initialize our progress bar for progress visualization - 7680 steps in total
            for k in range(0, M): #for each transmission/reception in a file, repeat the following
                T = a[(k*N_record):((k+1)*N_record):1,20] #T, the temperature, is equal to the temperature data in the range of the kth transmission/reception

                for i in range(0, S): #for each speaker in the array, repeat the following
                    x = a[(k*N_record):((k+1)*N_record):1,i] #isolates a single record of speaker data, depending on i, the speaker

                    if upSampleFlg: #if you would like to upsample the incoming data
                        x = sps.resample_poly(x, upSampleFlg, 1) #upsample the record data by upSampleFlg
                    
                    s, t1 = signalOnSpeaker(x,N_signal1,t1Exp[i],searchLag_S) #send the (potentially upsampled) record data, the length of the transmitted signal in samples * upSampleFlg, the expected speaker signal location for that speaker, and searchLag_S to signalOnSpeaker

                    for j in range(0, R): #for each microphone in the array, repeat the following
                        ttExp = getExpt2(i,j,T)
                        t2Exp = np.round(ttExp/dts)+t1
                        x = a[((k)*N_record):((k+1)*N_record):1,(S+j)]

                        if filterFlg: #if, when initializing the tt_extractor function, the filterFlg input existed
                            y = sp.fftpack.fft(x) #returns the fourier transform of the real x array in the form of a complex sequence
                            for p in range(0, f[0]+1): #eliminates the frequencies to be filtered between 0 and f[0]
                                y[p] = 0
                            for p in range(f[1], f[2]+1): #eliminates the frequencies to be filtered between f[1] and f[2]
                                y[p] = 0
                            for p in range(f[3], len(y)): #eliminates the frequencies to be filtered between f[3] and the end of the array
                                y[p] = 0
                            x = sp.fftpack.ifft(y) #returns the inverse fourier transform of the complex y array in the form of a complex (but nearly real) sequence
                            x = [np.real(p) for p in x] #since the output of ifft is not perfectly real, but instead very close to real (Xe-16), here we truncate the imaginary portion of the numbers so they are usable

                        if upSampleFlg: #if, when initializing the tt_extractor function, the upSampleFlg input existed
                            x = sps.resample_poly(x, upSampleFlg, 1) #resample the x array by upSampleFlg/1. If upSampleFlg = 10, the array is upsampled by 10

                        r, t2 = signalOnMic(x,s,t2Exp,searchLag)

                        s_r[:,k,i,j] = r
                        tt[i*R+j,k] = (t2-t1)*dts #write the travel time calculated to the current transmission/reception and speaker/microphone pair
                        pbar.update(1) #update the progress bar one step
                        
    np.save('s_rout', s_r) #save the finished 4D array to a .npy file
    np.savetxt('ttout.txt', tt, delimiter='\t', fmt='%.3f') #save the finished tt array to a .txt file


def getExpt2(i,j,T): #imports the current speaker (i) and microphone (j) as well as the temperature array for the current signal period
    xyR = np.array(([1.4579,39.7960],[33.4493,39.1313],[39.7780,14.6969],[37.8724,-28.0847],[16.7284,-40.0000],[-23.2493,-38.8741],[-40.0000,-14.8975],[-39.1533,27.7008])) #geographic locations of microphones
    xyS = np.array(([0.6579,39.7960],[32.6493,39.1313],[39.7780,15.4969],[ 37.8724,-27.2847],[17.5284,-40.0000],[-22.4493,-38.8741],[-40.0000,-15.6975],[-39.1533,26.9008])) #geographic locations of speakers

    gamma = 1.4 #theoretical constant
    R = 287.058 #theoretical constant
    T = [x + 272.15 for x in T] #converting temperature array from celcius to kelvin
    
    c0 = T #copy temperature array for conversion into speed of sound
    c0 = [x * R * gamma for x in c0] #multiply temperature array element-wise by R and gamma
    c0 = [math.sqrt(x) for x in c0] #square root array element-wise
    c0 = np.mean(c0) #the speed of sound over the interval is the mean of the array

    V0 = np.array([0, 0])

    ell = np.dot((xyS[i,:]-xyR[j,:]),(np.transpose((xyS[i,:]-xyR[j,:]))))
    ell = round(ell,3)
    ell = math.sqrt(ell)

    s = (-xyS[i,:]+xyR[j,:])
    s = [x / ell for x in s]
    s = [round(x) for x in s]
    np.transpose(s)
    tt = ell/(c0+np.dot(V0,s))*1000
    return tt


def signalOnSpeaker(x,N_signal,tExp,searchLag):
    Nx = len(x) #nx = the length of the record data

    if not tExp: #if there is no expected signal location within the record data
        signal_var = []

        for i in range(0,Nx-N_signal+1):
            signal_var = np.append(signal_var,np.sum(np.multiply((x[i:(i+N_signal-1):1]),(x[i:(i+N_signal-1):1]))))
        t1 = np.argmax(signal_var)
        s1 = x[t1:(t1+N_signal-1)]

    else: #if there is an expected signal location within the record data
        signal_var = []

        for i in range(-searchLag,searchLag+1):
            t1 = np.maximum(tExp+i,1)
            t2 = np.minimum(t1+N_signal-1,Nx)
            signal_var = np.append(signal_var,np.sum(np.multiply((x[t1:t2:1]),(x[t1:t2:1]))))
        
        i1 = np.argmax(signal_var)
        i = i1-searchLag-1
        t1 = np.maximum((i+1+tExp),1)
        s1 = x[t1-1:(t1+N_signal-1)]

    return s1, t1


def signalOnMic(x,s,tExp,searchLag):
    N_signal = len(s)
    Nx = len(x)

    if not tExp:
        signal_var = []

        for i in range(0,Nx-N_signal+1):
            signal_var = np.append(signal_var,np.sum(np.multiply((x[i:(i+N_signal-1):1]),s)))
        t1 = np.argmax(signal_var)
        tstart = np.maximum(t1-N_signal,1)
        tfinish = np.minimum(tstart+3*N_signal-1,Nx)
        s1 = x[tstart:tfinish]

    else:
        signal_var = []

        for i in range(-searchLag,searchLag+1):
            t1 = np.int(np.round(np.maximum(tExp+i,1)))
            t2 = np.int(np.round(np.minimum(t1+N_signal-1,Nx)))
            signal_var = np.append(signal_var,np.sum(np.multiply((x[t1:t2+1:1]),s)))

        i1 = np.argmax(signal_var)
        i = i1-searchLag+1
        t1 = np.maximum((i+tExp),1)
        tstart = np.int(np.round(np.maximum(t1-N_signal,1)))
        tfinish = np.int(np.round(np.minimum(tstart+3*N_signal-1,Nx)))
        s1 = x[tstart-1:tfinish]

    return s1, t1


tt_extractor(1,10,1) #function call of tt_extractor, starting the program