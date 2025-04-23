import numpy as np
import scipy
import matplotlib.pyplot as plt
import pywt
from numpy.fft import rfft, rfftfreq
import igrf
from skimage import feature
from scipy.ndimage import gaussian_filter
from skimage.transform import probabilistic_hough_line

########################## DATA EXTRACTION ##########################
# Given the filenames of data files from tuners 1 and 2, we use a short-time FFT to extract the signal intensities as a function of time.
# If there is no interference stronger than the rocket signal, this should correctly isolate the signal due to the rocket.
# Returns: time, amplitude of tuner 1, amplitude of tuner 2

# The freqcutoffstartind, freqcutoffendind parameters are an easy/hacky way to trim the frequency range of the FFT
# If there are interfering signals you wish to exclude
def extract_intensity(fn1,fn2,startind,endind,samprate,freqcutoffstartind,freqcutoffendind,plot=True,xstart=0,xend=None):

    # Read in files
    v = np.fromfile(fn1,dtype=np.complex64)[startind:endind+1]
    v2 = np.fromfile(fn2,dtype=np.complex64)[startind:endind+1]
    
    # Generate time coordinate for the raw signal
    sigtime = np.linspace(0,(len(v)-1)/samprate,len(v))
    
    # Define a window for short time fourier transform
    w = scipy.signal.windows.hamming(64)
    # Define STFT, zero padding fft to get a 512-frequency spectrum despite a 128-sample window
    SFT = scipy.signal.ShortTimeFFT(w, hop=32, fft_mode='centered', fs=samprate, mfft=1024, scale_to='magnitude')
    
    # Perform STFT on v and v2
    Sx = SFT.stft(v)
    Sx2 = SFT.stft(v2)
    
    # Extract frequencies
    freqs = SFT.f


    # freqs = freqs[300:-370]
    # Sx = Sx[300:-370,:]
    # Sx2 = Sx2[300:-370,:]
    freqs = freqs[freqcutoffstartind:freqcutoffendind]
    Sx = Sx[freqcutoffstartind:freqcutoffendind,:]
    Sx2 = Sx2[freqcutoffstartind:freqcutoffendind,:]

    
    # New time coordinate, corresponding to STFT of signal
    time = np.arange(SFT.extent(len(v))[0],SFT.extent(len(v))[1],SFT.delta_t)
    # New sampling period
    sampling_period = np.diff(time).mean()

    # Sum of the squares of the STFT's of both tuners
    Sxmag = np.abs(Sx)**2 + np.abs(Sx2)**2

    # Extracting frequency with maximum mag^2 signal value as a function of time
    bfi = np.argmax(Sxmag,axis=0)
    
    # Extract signal intensities now that we know the right frequency
    ampvec = np.abs(Sx[bfi,np.arange(len(bfi))])  
    ampvec2 = np.abs(Sx2[bfi,np.arange(len(bfi))])


    if plot:
        # Plot signal intensity of tuner 1 with time
        plt.pcolormesh(time[::100],freqs/1000,20*np.log10(np.abs(Sx[:,::100])),vmin=-70,vmax=-15,shading='nearest')
        plt.xlabel('time,seconds')
        plt.ylabel('freq diff from nominal, kHz')
        plt.title('Tuner 1')
        plt.colorbar()
        plt.xlim(xstart,xend)
        plt.show()
    
        # Plot mag^2 and best frequency with time
        plt.pcolormesh(time[::100],freqs/1000,10*np.log10(np.abs(Sxmag[:,::100])),vmin=-70,vmax=-15,shading='nearest')
        plt.xlabel('time,seconds')
        plt.ylabel('freq diff from nominal, kHz')
        plt.title('Mag^2')
        plt.colorbar()
        plt.plot(time,freqs[bfi]/1000,color='red')
        plt.xlim(xstart,xend)
        plt.show()

        # Plot the amplitudes with time
        plt.plot(time,ampvec,label='tuner 1')
        plt.plot(time,ampvec2,label='tuner 2')
    
        plt.xlabel('time,seconds')
        plt.ylabel('amplitude')
        plt.legend()
        plt.xlim(xstart,xend)
        plt.show()

        maxval = np.amax( [np.amax(ampvec[:int(np.round(5/sampling_period))]), np.amax(ampvec2[:int(np.round(5/sampling_period))])] )
        plt.plot(time,ampvec**2,label='tuner 1')
        plt.plot(time,ampvec2**2,label='tuner 2')
        plt.xlim(0,5)
        plt.ylim(0,maxval**2)

        plt.xlabel('time,seconds')
        plt.ylabel('amplitude**2')
        plt.legend()
        plt.show()

        maxval = np.amax( [np.amax(ampvec[-int(np.round(5/sampling_period)):]), np.amax(ampvec2[-int(np.round(5/sampling_period)):])] )
        plt.plot(time,ampvec**2,label='tuner 1')
        plt.plot(time,ampvec2**2,label='tuner 2')
        plt.xlim(time[-1]-5,time[-1])
        plt.ylim(0,maxval**2)
        plt.xlabel('time,seconds')
        plt.ylabel('amplitude**2')
        plt.legend()        
        plt.show()

        maxval = np.amax( [np.amax(ampvec[-int(np.round(5/sampling_period)):]), np.amax(ampvec2[-int(np.round(5/sampling_period)):])] )
        plt.plot(time,ampvec**2,label='tuner 1')
        plt.plot(time,ampvec2**2,label='tuner 2')
        plt.xlim(xstart,xend)
        plt.ylim(0,maxval**2)
        plt.xlabel('time,seconds')
        plt.ylabel('amplitude**2')
        plt.title('custom bounds')
        plt.legend()        
        plt.show()


    return time,ampvec,ampvec2

# Given two signals, assumed to roughly be of the form Env1(t)*sin(phi(t)),Env2(t)*cos(phi(t))
# (with some additive noise), returns roughly: sin(phi(t)),cos(phi(t)),Env1(t),Env2(t)
def normalize_amplitudes(time,ampvec,ampvec2,plot=True):
    # Extract peaks of both signals
    peaks = scipy.signal.find_peaks(ampvec,prominence=np.mean(ampvec)/15)[0]
    peaks2 = scipy.signal.find_peaks(ampvec2,prominence=np.mean(ampvec2)/15)[0]

    # Cubic Spline interpolation of peaks
    peaksscs = scipy.ndimage.gaussian_filter(scipy.interpolate.CubicSpline(time[peaks],ampvec[peaks])(time),300)
    peaks2scs = scipy.ndimage.gaussian_filter(scipy.interpolate.CubicSpline(time[peaks2],ampvec2[peaks2])(time),300)

    # Raw ratio of cubic splines
    rat = scipy.interpolate.CubicSpline(time[peaks],ampvec[peaks])(time)/scipy.interpolate.CubicSpline(time[peaks2],ampvec2[peaks2])(time)

    # 
    ampnorm = (ampvec/peaksscs) / np.sqrt( (ampvec/peaksscs)**2 + (ampvec2/peaks2scs)**2 )
    ampnorm2 = (ampvec2/peaks2scs) / np.sqrt( (ampvec/peaksscs)**2 + (ampvec2/peaks2scs)**2 )

    if plot:
        # Plot peaks + smoothed splines 
        plt.plot(time,ampvec)
        plt.scatter(time[peaks],ampvec[peaks],s=10,color='orange')
        plt.title('peaks, tuner 1')
        plt.show()
        
        plt.plot(time,ampvec2)
        plt.scatter(time[peaks2],ampvec2[peaks2],s=10,color='orange')
        plt.title('peaks, tuner 2')
        plt.show()
    
        plt.scatter(time[peaks],ampvec[peaks],s=10,color='blue')
        plt.scatter(time[peaks2],ampvec2[peaks2],s=10,color='orange')
        plt.plot(time,scipy.interpolate.CubicSpline(time[peaks],ampvec[peaks])(time),label='tuner1')
        plt.plot(time,scipy.interpolate.CubicSpline(time[peaks2],ampvec2[peaks2])(time),label='tuner2')
        
        plt.plot(time,peaksscs,color='blue')
        plt.plot(time,peaks2scs,color='red')
        
        plt.title('Signal Power Peaks')
        plt.xlabel('time')
        plt.ylabel('amplitude, arbitrary units')
        plt.legend()
        plt.show()
    
        # Plot ratio + smoothed ratio
        plt.scatter(time,rat,s=0.1,label='raw')
        plt.plot(time,peaksscs/peaks2scs,color='orange',label='gaussian smoothed')
        
        plt.title('Ratio of amplitude, tuner1 to tuner2')
        plt.legend()
        plt.xlabel('time, seconds')
        plt.ylim(0.,5.)
        plt.show()

        # Plot normalized amplitudes
        plt.plot(time,(ampnorm)**2,label='tuner 1')
        plt.plot(time,(ampnorm2)**2,label='tuner 2')
        plt.title('corrected squared amplitudes')
        plt.legend()
        plt.xlabel('time, seconds')
        plt.xlim(0.,5.)
        plt.show()

        plt.plot(time,(ampnorm)**2,label='tuner 1')
        plt.plot(time,(ampnorm2)**2,label='tuner 2')
        plt.title('corrected squared amplitudes')
        plt.legend()
        plt.xlabel('time, seconds')
        plt.xlim(time[-1]-5,time[-1])
        plt.show()


    return ampnorm,ampnorm2,peaksscs,peaks2scs

# Finds sudden jumps in phase, asssumes they are due to domain changes of trig functions and corrects them
def domaincorrect(vec):
    # Copy input vector
    newvec = np.copy(vec)
    # Find large jumps
    jumps = np.where(np.abs(np.diff(vec))>1.)[0]
    for jump in jumps:
        # Shift either up or down by 2pi
        if np.diff(vec)[jump]>0:
            newvec[jump+1:] -= 2*np.pi
        elif np.diff(vec)[jump]<0:
            newvec[jump+1:] += 2*np.pi
    # Return corrected vector
    return newvec

# Uses a continuous wavelet transform to extract the time-varying phase of a signal, minus a reference
def wavelet_extractphase(time,signal,refsig,plot=True):
    # Make sure signal has mean zero
    chirp = np.copy(signal)-np.mean(signal)
    ######## perform CWT ########
    ################################
    wavelet = "cmor1.5-1.0"
    # logarithmic scale for scales, as suggested by Torrence & Compo:
    #widths = np.geomspace(sr/3, sr/2, num=200)
    widths = np.geomspace(30, 50, num=100)

    # Determine sampling period as an input to CWT
    sampling_period = np.mean(np.diff(time))
    # CWT of signal
    cwtmatr, freqs = pywt.cwt(chirp, widths, wavelet, sampling_period=sampling_period)
    # CWT of ref
    cwtmatr_ref, _ = pywt.cwt(refsig, widths, wavelet, sampling_period=sampling_period)
    ################################
    ################################

    # Find phase difference
    phasediff = np.angle((cwtmatr/cwtmatr_ref))[np.argmax(np.abs(cwtmatr),axis=0),np.arange(len(signal))]
    if plot:
        # Plot the CWT mag + phase diff
        
        plt.pcolormesh(time[::100], freqs, np.abs(cwtmatr)[:,::100])
        plt.xlabel('time')
        plt.ylabel('freq')
        plt.title('magnitude of cwt of recorded signal')
        plt.plot(time,freqs[np.argmax(np.abs(cwtmatr),axis=0)],color='red')
        
        plt.colorbar()
        #plt.ylim(0,4)
        plt.show()
        
        plt.title('phase difference, radians, measured signal from estimated reference')
        plt.pcolormesh(time[::100], freqs, (np.angle(cwtmatr/cwtmatr_ref))[:,::100],vmin=-np.pi,vmax=np.pi)
        plt.colorbar()        
        plt.scatter(time,freqs[np.argmax(np.abs(cwtmatr),axis=0)],color='red',s=0.01)
        plt.xlabel('time')
        plt.ylabel('freq')
    # Return phase diff
    return domaincorrect(phasediff)

# Uses a continuous wavelet transform to extract the time-varying phase of a signal, minus a reference
def wavelet_extractphase_uniformfreq(time,signal,refsig,plot=True):
    # Make sure signal has mean zero
    chirp = np.copy(signal)-np.mean(signal)
    ######## perform CWT ########
    ################################
    wavelet = "cmor1.5-1.0"
    # logarithmic scale for scales, as suggested by Torrence & Compo:
    #widths = np.geomspace(sr/3, sr/2, num=200)
    widths = np.geomspace(30, 50, num=100)

    # Determine sampling period as an input to CWT
    sampling_period = np.mean(np.diff(time))
    # CWT of signal
    cwtmatr, freqs = pywt.cwt(chirp, widths, wavelet, sampling_period=sampling_period)
    # CWT of ref
    cwtmatr_ref, _ = pywt.cwt(refsig, widths, wavelet, sampling_period=sampling_period)

    # Peak frequency of reference signal
    medfreq = np.median(freqs[np.argmax(np.abs(cwtmatr_ref),axis=0)])
    bestind = np.argmin(np.abs(medfreq-freqs))
    ################################
    ################################

    # Find phase difference
    #phasediff = np.angle((cwtmatr/cwtmatr_ref))[np.argmax(np.abs(cwtmatr),axis=0),np.arange(len(signal))]
    phasediff = np.angle((cwtmatr/cwtmatr_ref))[bestind,:]

    if plot:
        # Plot the CWT mag + phase diff
        
        plt.pcolormesh(time[::100], freqs, np.abs(cwtmatr)[:,::100])
        plt.xlabel('time')
        plt.ylabel('freq')
        plt.title('magnitude of cwt of recorded signal')
        plt.plot(time,freqs[np.argmax(np.abs(cwtmatr),axis=0)],color='red')
        plt.plot(time,time*0+freqs[bestind],color='blue')
        
        plt.colorbar()
        #plt.ylim(0,4)
        plt.show()
        
        plt.title('phase difference, radians, measured signal from estimated reference')
        plt.pcolormesh(time[::100], freqs, (np.angle(cwtmatr/cwtmatr_ref))[:,::100],vmin=-np.pi,vmax=np.pi)
        plt.colorbar()        
        plt.plot(time,time*0+freqs[bestind],color='blue')
        plt.xlabel('time')
        plt.ylabel('freq')
        plt.show()
    # Return phase diff
    return domaincorrect(phasediff)


################### FFT ANALYSIS #####################
# Given a function, sample period, and a set of bounds, returns the fourier component (frequency) with the largest amplitude
# Within the listed bounds
def fftpeak(chirp,tsamp,minbound,maxbound,plot=True):
    # FFT of function
    yf = rfft(chirp,len(chirp)*50)
    # Corresponding frequencies
    xf = rfftfreq(len(chirp)*50, tsamp)
    # print(xf)
    if plot:
        plt.scatter(xf/2, np.abs(yf),s=0.01)
        #print(minbound)
        #print(maxbound)
        plt.xlim(minbound/2,maxbound/2)
        #plt.show()
    # Largest amplitude component within bounds
    bestfreq = xf[np.where((minbound<xf) & (maxbound > xf))][np.argmax(np.abs(yf[np.where((minbound<xf) & (maxbound > xf))]))]
    return bestfreq
# Takes in np arrays time, function(time)
# Returns time_1 (new time coord), freqvec(time_1)
# Where freqvec is the largest amplitude FFT component as a function of time_1
def stft_times(time,fin):
    # Window
    w = scipy.signal.windows.hamming(4096)
    # Find STFT
    SFT = scipy.signal.ShortTimeFFT(w, hop=256, fs=1/np.mean(np.diff(time)), mfft=4096*16, scale_to='magnitude')
    Sx = SFT.stft(fin)
    # Extract frequencies, define new time vector
    freqs = SFT.f
    timevec = np.linspace(SFT.extent(len(fin))[0],SFT.extent(len(fin))[1],Sx.shape[1])
    return timevec,freqs[np.argmax(np.abs(Sx),axis=0)]

################# COORDINATE + MAG FIELD FUNCTIONS #####################
# takes in geodetic lat,lon,alt and spits out new cartesian coords centered at venitie
def vcvf(lat,lon,alt):
    veelat = 67.0161
    veelon = -146.4195
    re = 6360.036

    newx = (lat - veelat)*(np.pi/180)*re
    newy = (veelon - lon)*(np.pi/180)*re*np.sin((90-lat)*np.pi/180)
    newz = alt

    return newx,newy,newz

# Goes from VCVF coords back to geodetic
def vcvf_to_lla(x,y,z):
    veelat = 67.0161
    veelon = -146.4195
    re = 6360.036

    alt = z
    lat = (x/((np.pi/180)*re))+veelat
    lon = veelon - y/(re*(np.pi/180)*np.sin((90-lat)*np.pi/180))
    return lat,lon,alt

# 3d linspace between two (x,y,z) coord vectors
# Returns xvec,yvec,zvec
def linspace3d(vec1,vec2,npts):
    xvec = np.linspace(vec1[0],vec2[0],npts)
    yvec = np.linspace(vec1[1],vec2[1],npts)
    zvec = np.linspace(vec1[2],vec2[2],npts)
    return xvec,yvec,zvec
    
# Given a coordinate in VCVF and a date, calculates the vector B field there in VCVF coords using IGRF
def magvec_vcvf(x,y,z,date):
    lat,lon,alt = vcvf_to_lla(x,y,z)
    magstruct = igrf.igrf(date,glat=lat,glon=lon,alt_km=alt)
    Bmag = np.asarray(magstruct.total)[0]
    inc = np.asarray(magstruct.incl)[0]
    decl = np.asarray(magstruct.decl)[0]
    Bhat = np.asarray([-np.cos(inc*np.pi/180)*np.cos(decl*np.pi/180),-np.cos(inc*np.pi/180)*np.sin(decl*np.pi/180),np.sin(inc*np.pi/180)])
    B = Bmag*Bhat
    return B

# Given [x,y,z] start, [x,y,z] end (in VCVF coords), and npts, returns the average value of B// and the length of the path
def path_Bpar(posvec0,posvec1,npts,date,plot=True):
    # Constructing coordinates of path
    xv,yv,zv = linspace3d([posvec0[0],posvec0[1],posvec0[2]],[posvec1[0],posvec1[1],posvec1[2]],npts)
    
    # unit vector pointing along path
    dlv = np.asarray([xv[1]-xv[0],yv[1]-yv[0],zv[1]-zv[0]])
    dlv /= np.linalg.norm(dlv)

    # length of path
    L = np.linalg.norm(posvec1-posvec0)

    # Finding along-path component of B at each point in the path
    bparvec = []
    if plot:
        plt.title('b//')
    for i in range(len(xv)):
        bpar = np.dot(magvec_vcvf(xv[i],yv[i],zv[i],date),dlv)
        bparvec.append(bpar)
        if plot:
            plt.scatter(i,bpar)
    return np.mean(bparvec),L

################### ROGUE SIGNAL DETECTION #####################
class boundary_line:
    def __init__(self, x0,x1,y0,y1, intercept):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.intercept = intercept
def get_line_func(x0,y0,x1,y1):
    slope = (y1 - y0) / (x1 - x0)
    intercept = y0 - slope * x0
    return lambda x: slope * x + intercept
def pixel2real(line, time_arr, freq_arr):
    line.x0 = time_arr[line.x0]
    line.x1 = time_arr[line.x1]
    line.y0 = freq_arr[line.y0]
    line.y1 = freq_arr[line.y1]
    # convert to khz range
    line.y0 = line.y0/1000
    line.y1 = line.y1/1000
    return line


def detect_rogue_crossing(stft, stft2, time_arr, freq_arr, sample_rate, padding_time=5):
    S = np.abs(stft)
    S2 = np.abs(stft2)
    bfi = np.argmax(S**2+S2**2,axis=0)
    pfit = np.polynomial.polynomial.Polynomial.fit(time_arr, freq_arr[bfi]/1000, 3)
    fitted = pfit(time_arr)
    S_dB = 20 * np.log10(S + 1e-12)  # shift to avoid log(0)
    S_blurred = gaussian_filter(S_dB, sigma=25)
    edges = feature.canny(S_blurred, sigma=3)
    lines = probabilistic_hough_line(edges, threshold=10, line_length=3000, line_gap=300)
    filtered_lines = []
    for i in lines:
        x0=i[0][0]
        y0=i[0][1]
        x1=i[1][0]
        y1=i[1][1]
        if np.abs((y1-y0)/(x1-x0))>=0.01:
            filtered_lines.append(i)
    line_list = []
    for i in filtered_lines:
        tmp = boundary_line(i[0][0],i[1][0],i[0][1],i[1][1],(i[0][1]-i[1][1])/(i[0][0]-i[1][0]))
        line_list.append(tmp)
    # min and max is reversed for some reason
    # min_boundary = min(line_list, key=lambda x: x.intercept)
    # max_boundary = max(line_list, key=lambda x: x.intercept)
    # min_boundary = pixel2real(min_boundary, time_arr, freq_arr)
    # max_boundary = pixel2real(max_boundary, time_arr, freq_arr)
    func_list = []
    for i in line_list:
        i = pixel2real(i, time_arr, freq_arr)
        func_list.append(get_line_func(i.x0,i.y0,i.x1,i.y1))
    # plot the boundary over the spectogram
    # real_min_func = get_line_func(min_boundary.x0,min_boundary.y0,min_boundary.x1,min_boundary.y1)
    # real_max_func = get_line_func(max_boundary.x0,max_boundary.y0,max_boundary.x1,max_boundary.y1)
    data_list = []
    for i in func_list:
        data_list.append(i(time_arr))
    # min_data = real_min_func(time_arr)
    # max_data = real_max_func(time_arr)
    # min_diff = fitted-min_data
    # max_diff = fitted-max_data
    crossings = []
    for i in data_list:
        diff = fitted-i
        crossings.append(np.where(np.diff(np.sign(diff)))[0][0])
    # min_crossing = np.where(np.diff(np.sign(min_diff)))[0][0]
    # max_crossing = np.where(np.diff(np.sign(max_diff)))[0][0]
    #padding_idxs = int(padding_time*sample_rate)
    #return real_min_func, real_max_func
    #return np.array([max_crossing-padding_idxs,min_crossing+padding_idxs])
    return np.array([np.min(crossings),np.max(crossings)])

