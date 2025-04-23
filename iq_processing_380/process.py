################## Script for automatic IQ file pre-processing - synchronization and decimation + filtering ##################
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

