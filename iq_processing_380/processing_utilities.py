import re
import subprocess
import os

#object to store receiver metadata
class receiver:
    def __init__(self, location: str, tuner_num: int, freq: float, sample_rate: float, if_atten: float, starttime: float, filename: str):
        self.location = location
        self.tuner_num = tuner_num
        self.freq = freq
        self.sample_rate = sample_rate
        self.if_atten = if_atten
        self.starttime = starttime
        self.filename = filename
        delay=0.0
        start_ind = 0
        filepath = ""
        sinkpath = ""


# get locations based on subdirectory names
def get_locs(dir):
    tmp = subprocess.check_output(['ls '+dir],shell=True, text=True).strip()
    return tmp.splitlines()

# generate list of receiver objects
def get_receivers(dir, locs):
    receivers = []
    for loc in locs:
        filenames = []
        tmp = subprocess.check_output(['ls '+dir+"/"+loc],shell=True, text=True).strip()
        lines = tmp.splitlines()

        # get the data filenames
        for i, element in enumerate(lines):
            if "150M" in element:
                filenames.append(element)

        # get the last word of the directory, corresponds to the tuner location
        tuner_location = loc
        # the following is AI generated because I have no patience for regex
        for name in filenames:
            freq_match = re.search(r'(\d+)([kMG]?)[bB]?', name)
            unit_multipliers = {'k': 1e3, 'M': 1e6, 'G': 1e9}
            freq = int(freq_match.group(1)) * unit_multipliers.get(freq_match.group(2), 1)

            # Initialize values
            tuner_num = 0
            sample_rate = 0.0
            if_atten = 0.0
            starttime = 0.0

            # Special handling for tuner#
            tuner_match = re.search(r'tuner(\d+)', name)
            if tuner_match:
                tuner_num = int(tuner_match.group(1))

            # Extract parameters that use '='
            pairs = re.findall(r'([a-zA-Z_]+)=(\d+(?:\.\d+)?)', name)
            for key, value in pairs:
                value = float(value) if '.' in value else int(value)
                key = key.strip('_').lower()

                if key == 'samprate':
                    sample_rate = value
                elif key == 'ifatten':
                    if_atten = value
                elif key == 'starttime':
                    starttime = value

            # Create and append the receiver object
            receivers.append(receiver(
                location=tuner_location,
                tuner_num=tuner_num,
                freq=freq,
                sample_rate=sample_rate,
                if_atten=if_atten,
                starttime=starttime,
                filename=name
            ))
            for i in receivers:
                i.filepath = dir + "/" + i.location + "/" + i.filename
    return receivers

def get_earliest_start(receivers):
    return max(receivers, key=lambda r: r.starttime).starttime

def set_delays(receivers):
    first_start = get_earliest_start(receivers)
    for i in receivers:
        i.delay = first_start - i.starttime
        i.start_ind = int(i.delay * i.sample_rate)

def make_target_direc(dir, name):
    # Create the target directory if it doesn't exist
    try:
        os.mkdir(dir+'/'+name)
        print("Processed directory created at: " + dir+'/'+name)
    except FileExistsError:
        print("Processed directory already exists at: " + dir+'/'+name)

def make_sinks(receivers, dir, target_name):
    for i in receivers:
        filepath = os.path.join(dir, target_name, f"{i.location}{i.tuner_num}_processed.bin")
        # ensure target exists
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w'):
            pass
        i.sinkpath = filepath
