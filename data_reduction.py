"""
NAM: Data Reduction (module)
DES: contains the Data analysis methods for the MDAP pipeline.
"""

# ----------------------------------------------------------------------------------
# ----- Import ---------------------------------------------------------------------
# ----------------------------------------------------------------------------------

from astropy.io import ascii
import numpy as np
import os
import tkinter as tk
from tkinter import filedialog
import concurrent.futures
import pickle
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------
# ----- File Selection class -------------------------------------------------------
# ----------------------------------------------------------------------------------

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    RED = "\033[1;31m"
    BLUE = "\033[1;34m"
    CYAN = "\033[1;36m"
    GREEN = "\033[0;32m"
    RESET = "\033[0;0m"
    BOLD = "\033[;1m"
    REVERSE = "\033[;7m"

class BinaryConvert:
    """
    Contains modules for the BINARY CONVERSION and Artifact Removal (Cleaning)
    """
    def __init__(self):
        self.path_loc = []
        print(f"{bcolors.BOLD}Welcome!{bcolors.ENDC}")
        print(f"Start by defining a path to the folder with the MKID readout data using {bcolors.RED}select_files(){bcolors.ENDC} method")

    def select_files(self):
        root = tk.Tk()
        root.withdraw()

        file_path = filedialog.askopenfilenames()
        return file_path


    def average_sweeps(self, filename):
        """
        filename: file index.

        Averages the up and down sweep values.
        This is to be used internally.

        """

        data = ascii.read(filename)
        data_av = (data['col3'] + data['col4']) / 2.  # averaging the two sweeps
        print('File {} has been averaged'.format(filename))
        return data_av


    def save_as_pickle(self, fsp_files, filename):

        # TODO: Add antlog and obs, time shift files as well into Xarray format.

        """
        fsp_files: the list of MKID frequency shift raw files to be used.
        filename: name of the file to be created.

        Averages the up and down frequency sweep values values.
        """

        average_values = []

        with concurrent.futures.ProcessPoolExecutor() as executor:
            for prime in executor.map(self.average_sweeps, fsp_files):
                average_values.append(prime)



        with open(str(filename)+'.pickle', 'wb') as handle:
            pickle.dump(average_values, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print(f"File {filename} has been created successfully!")

    def load_pickle(self, filename):
        """
        filename: name of the pickle file to be loaded.

        Loads the pickle file.
        """

        with open(str(filename), 'rb') as handle:
            data = pickle.load(handle)

        return data

    def plot(self, data, id):
        """
        data: the averaged data from fsp.
        id: the index of the data to be plotted

        plots the data for visualization.
        """

        plt.plot(data[id-1])
        plt.show()