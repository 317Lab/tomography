# import numpy as np
# from scipy import signal
# import iqtools
# import matplotlib.pyplot as plt
# # filename='381/lav_tuner2_3point125k_25atten_yellowstart.raw'
# # iq_data = np.fromfile(filename, dtype=np.int16).astype(np.float32)
# # iq_complex = iq_data[::2] + 1j * iq_data[1::2]
# # iq_complex.astype(np.complex64).tofile("lav_tuner2_3point125k_25atten_yellowstart.bin")
# filename='lav_tuner2_3point125k_25atten_yellowstart.bin'
# # iq = get_iq_object(filename)
# # iq.read_samples(200*1024)
# # iq.method='fftw'
# # time_grid,freq_grid,power_grid = iq.get_power_spectrogram(lframes = 1024, nframes = 200, sparse=True)
# iqdata = iqtools.GRData(filename, fs = 2.5e6, center=30e6)
# iqdata.read_samples(2000*1024)
# xx, yy, zz = iqdata.get_power_spectrogram(nframes=2000, lframes=1024)
# iqtools.plot_spectrogram(xx, yy, zz, filename='testplt')
 
import argparse
import processing_utilities as util
import subprocess

target_name = 'processed'
# takes flight directory with location subdirectories
parser = argparse.ArgumentParser()
parser.add_argument('unprocessed_directory', type=str)
options = parser.parse_args()
dir = options.unprocessed_directory
print("processing metadata...")
locs = util.get_locs(dir=dir)
receivers = util.get_receivers(locs=locs, dir=dir)
print("setting synchronization delays...")
util.set_delays(receivers)
print("creating directory for processed data...")
util.make_target_direc(dir=dir, name=target_name)
util.make_sinks(receivers=receivers, dir=dir, target_name=target_name)
for i in receivers:
    print("processing "+i.location+" Tuner "+str(i.tuner_num))
    subprocess.run(['python3', 'decsink.py', '--source-name='+i.filepath, '--sink-name='+i.sinkpath,'--startind='+str(i.start_ind)])

