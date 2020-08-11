"""
NAM: Data Reduction (module)
DES: contains the Data analysis methods for the MDAP pipeline.
"""

# ----------------------------------------------------------------------------------
# ----- Import ---------------------------------------------------------------------
# ----------------------------------------------------------------------------------
import Utilities as ut

from astropy.io import ascii
import numpy as np
import os
import tkinter as tk
from tkinter import filedialog
import concurrent.futures
import pickle
import matplotlib.pyplot as plt

from scipy.signal import argrelextrema, argrelmin, argrelmax
import heapq

from astropy.stats import sigma_clip

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
        self.mkid_exclude = []
        self.idMKIDs_use = self.num_of_pixels()
        self.len_chunks = []

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


    def save_fsp_pickle(self, fsp_files, filename):

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

    def save_as_pickle(self, data, filename):
        """
        fsp_files: the list of MKID frequency shift raw files to be used.
        filename: name of the file to be created.

        Averages the up and down frequency sweep values values.
        """

        with open(str(filename)+'.pickle', 'wb') as handle:
            pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print(f"File {filename} has been created successfully!")

    def load_pickle(self, filename):
        """
        filename: name of the pickle file to be loaded.

        Loads the pickle file.
        """

        with open(str(filename), 'rb') as handle:
            data = pickle.load(handle)

        self.num_of_pixels(len(data)) # automatically creates num pixels
        return data

    def num_of_pixels(self, num=101):
        """
        num: the number of pixels in the dataset

        create a list of pixels to be used internally, does not return anything.
        If not defined, exclusively, will use the length of data as number of pixels.
        """
        return np.arange(1, num+1, 1)



    def plot(self, data, id):
        """
        data: the averaged data from fsp.
        id: the index of the data to be plotted

        plots the data for visualization.
        """
        # TODO: make a window with plot where one can decide whether to use it or discard the TOD.

        plt.title(f"TOD of MKID{id}")
        plt.plot(data[id-1], ',')
        plt.show()

    def pixels_to_exclude(self, ids):
        """
        ids: ids of MKIDs to be excluded
        takes in ids and identifies bad mkids such that it is not used for data analysis stages.
        """

        try:
            if type(ids) is list :
                self.mkid_exclude = ids
            else:
                raise TypeError
        except TypeError:
            print("the ids should be a list")

# ------------------  FREQUENCY IDENTIFICATION CODE--------------------

    def identify_frequencies(self, data, minima_order=4000, exclude=True, plot=False, save_file=True):
        """
        data: is the TOD of frequency shift including the Hot-Rod reading.

        return: file with frequency shift values, can save it as well as well.
        """
        leave = []
        frequencies = []

        if exclude:
            leave = self.mkid_exclude
        else:
            leave = []

        for mkidid in range(len(data)):
            if mkidid + 1 in leave:
                print(mkidid + 1, 'excluded')
                frequencies.append([mkidid + 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

            else:
                l1 = data[mkidid]

                minima_index, minima = self.minima_func(l1, minima_order)  # 400 is chosen as an optimum order of argrelextrema to get
                # just one minima at the load while recognising signal minimas near the load

                fon, fonindex, fload, floadindex = self.indentify(minima_index, minima)

                if plot:
                    plt.plot(l1)
                    plt.plot(minima_index, minima, 'o', color='purple')
                    plt.plot(fonindex, fon, 'o', color='red', label='f_brightest')
                    plt.plot(floadindex, fload, 'o', color='green', label='f_load')
                    plt.legend()
                    plt.title("MKID{:03}".format(mkidid + 1))
                    # plt.savefig('MKID{:03}frequency_detected.png'.format(mkidid + 1))
                    plt.show()
                    plt.close()

                f_off = np.amax(l1)
                f_off_index = [i for i, j in enumerate(l1) if j == f_off]
                f_base = self.fbase_load_subtractor(l1, fload, floadindex)
                f = [mkidid + 1, fon[0], fonindex[0], fload[0], floadindex[0], fload[1], floadindex[1], f_off,
                     f_off_index[0], f_base[0], f_base[1]]
                frequencies.append(f)

        if save_file:
            np.savetxt('frequencies.txt', frequencies, fmt='%s')
            print("the Frequencies have been saved into a file titled \"Frequencies.txt\"")
        return frequencies


    def fbase_load_subtractor(self, data, fload, floadindex):
        fbase1 = np.average(data[floadindex[0] + 2000:floadindex[0] + 3000])
        fbase2 = np.average(data[floadindex[1] - 3000:floadindex[1] - 2000])
        fbase = [fbase1 - fload[0], fbase2 - fload[1]]
        # print( fbase, 'FBASE')
        return fbase

    def minima_func(self, data, order):
        minima_index = argrelextrema(data, np.less, order=order)
        minima = [data[i] for i in minima_index]
        # print minima[0], minima_index[0]

        return minima_index[0], minima[0]

    def indentify(self, minima_index, minima):
        fload = [heapq.nsmallest(1, minima)[-1], heapq.nsmallest(2, minima)[-1]]
        # print fload, 'fload'
        floadindex = [minima_index[i] for i, x in enumerate(minima) if x in fload]

        fon = [heapq.nsmallest(3, minima)[-1]]
        fonindex = [minima_index[i] for i, x in enumerate(minima) if x in fon]

        # print fon, fonindex, fload, floadindex

        return fon, fonindex, fload, floadindex

    def average(self, data, value):
        l = len(data)
        print
        'total number of data points are', l
        ############Reducing the data points by taking average-----------------------
        ml = divmod(l, value)
        print
        'to make an average over 100 data point,', ml[1], 'data from the end are rejected'
        data = data[:-ml[1]]
        x = np.array(data)
        nx = np.mean(x.reshape(-1, value), axis=1)  # averages it out over 100 points
        print
        'the total data points after averaging are', len(nx)

        return nx


# ------------------ CLEANING THE DATASET --------------------


    def outlier_remove(self, data):
        """
        data: The data from which the outlier is to be removed.

        Automatically detects and removes outliers.
        """
        # TODO: Outliers

        pass

    def changepoint_remove(self, data, patch, jump_at):
        """
        data: The data from which the outlier is to be removed.

        Automatically detects and removes outliers.
        """
        # TODO: Actual Changepoint detection using an algorithm

        for i in range(len(data)):
            f = data[i]

            mean_left = np.mean(f[jump_at - patch:jump_at])
            mean_right = np.mean(f[jump_at: jump_at + patch])
            diff_mean = mean_left - mean_right

            f[jump_at:-1] = f[jump_at:-1] + diff_mean
            data[i] = f

        return data

# ------------------ Decorrelation --------------------

    def clip(self, data, starting_index, ending_index):
        """
        Clips the data for decorrealtion
        """
        res =[]
        for idMKID in self.idMKIDs_use:
            if idMKID not in self.mkid_exclude:
                raw_data = data[idMKID-1]
                raw_data = np.array(raw_data, dtype=float).transpose()
                dat_FS = raw_data[starting_index:ending_index]
                res.append(dat_FS)
        res = np.array(res)  # Very important
        return res

    def mean_center(self, data):
        """
        removes the mean from each TOD
        """

        for i in range(len(data)):
            m = np.mean(data[i])
            data[i] = data[i] - m
        return data

    def pca_decor(self, data, n_components=2):
        """
        Decorrelates using PCA function
        """

        res = ut.work_pca(data, n_components)
        res = np.array(res)
        return res

    def flatten(self, pca_data, index, save=False):
        flat = []
        for i in range(len(pca_data)):
            flat.append(pca_data[i][:, index])
        flat = np.array(flat)

        # TODO: see if flat is top be saved or not
        return flat

    def get_sigmaclipped_mask(self, data, avg_points=32, sigma=4, maxiters=5, axis = 1):
        """
        Returns mask
        """

        running_avged = []
        for i in range(len(data)):
            running_avged.append(np.convolve(data[i], np.ones((avg_points,)) / avg_points, mode='same'))
        print(np.shape(running_avged))
        # running_avged = normalize(running_avged)

        filtered_data = sigma_clip(running_avged, sigma=sigma, maxiters=maxiters, axis=axis)
        mask = np.ma.getmask(filtered_data)
        return mask

    def make_chunk_matrix(self, mask, num_chunks):
        """
        returns the chunk matrix
        """
        num_chunks = num_chunks
        len_chunks = int(mask.shape[1] / num_chunks)  # TODO: make it take fraction as well
        self.len_chunks = len_chunks
        len_remaining = mask.shape[1] % num_chunks  # see if there are left over data

        if len_remaining != 0:
            chunk_matrix = np.zeros((len(mask), num_chunks + 1))
            for i in range(len(mask)):
                for j in range(num_chunks):
                    dot_val = np.prod(~mask[i][(j * len_chunks):(len_chunks * (j + 1))])
                    chunk_matrix[i][j] = (~dot_val + 2)
                chunk_matrix[i][num_chunks] = (~(np.prod(~mask[i][(-len_chunks):-1])) + 2)
            return chunk_matrix
        elif len_remaining == 0:
            chunk_matrix = np.zeros((len(mask), num_chunks))
            for i in range(len(mask)):
                for j in range(num_chunks):
                    dot_val = np.prod(~mask[i][(j * len_chunks):(len_chunks * (j + 1))])
                    chunk_matrix[i][j] = (~dot_val + 2)
            return chunk_matrix


    def chunkpca_decor(self, mask, chunk_matrix, data):
        """
        returns decorrelated data using ChunkPCA
        """
        t_chunk = chunk_matrix.T
        len_chunks = self.len_chunks
        for i in range(len(t_chunk)):  # range(len(t_chunk))
            pooled_data = []
            full = []
            empty = []
            for j in range(len(t_chunk[i])):  # LET'S CALL IT THE "POOLING METHOD"
                if t_chunk[i][j] == 0:
                    val = data[j][(i * len_chunks):(len_chunks * (i + 1))]
                    pooled_data.append(val - np.mean(val))  # mean centering the pooled_data
                    full.append(j)
                elif t_chunk[i][j] == 1:
                    empty.append(j)
            print("chunk {} pool data length is {}".format(i, np.shape(pooled_data)))
            result_masked_pca = ut.work_pca(pooled_data, n_components=3)  # pca on the pool
            result_masked_pca = np.array(result_masked_pca)

            ##### TODO: make sure the empty doesn't have too much of the pixels, since then the pca won't work.
            # print("unused pixels in chunk {}".format(i + 1))
            # print(empty)

            ############# Finding the baseline using the most correlated components from the pooled_data #############

            correlated = []

            for q in empty:
                val = data[q][(i * len_chunks):(len_chunks * (i + 1))]
                correlation_matrix = np.corrcoef(np.vstack((pooled_data, val)))
                fpr = []
                for w, t in enumerate(correlation_matrix[-1]):
                    if t > 0:
                        fpr.append(w)
                correlated.append([q, fpr[:-1]])

            ############# Putting the missing pixels values with proper data, baseline #############

            for k in empty:
                val = data[k][(i * len_chunks):(len_chunks * (i + 1))]
                masked_val = np.ma.masked_array(data[k][(i * len_chunks):(len_chunks * (i + 1))],
                                                mask[k][(i * len_chunks):(len_chunks * (i + 1))])
                val = val - np.mean(masked_val)
                for x in correlated:
                    if x[0] == k:
                        empty_baseline_pool = []
                        for p in x[1]:
                            if p not in empty:
                                got = result_masked_pca[p][:, 1]
                                empty_baseline_pool.append(got)
                        baseline = np.mean(empty_baseline_pool, axis=0)
                dat = val - baseline
                average = np.stack([dat, baseline, val], axis=1)
                result_masked_pca = np.insert(result_masked_pca, k, average, 0)

            # print("Final", np.shape(result_masked_pca))

            ##### JOINING ALL THE CHUNKS ######

            if i == 0:
                final = result_masked_pca
            else:
                final = np.concatenate([final, result_masked_pca], axis=1)
        return final

# ------------------ Calibration --------------------

    def load_frequency(self):
        """
        loads and returns frequency file
        """

        try:
            freq = np.loadtxt("frequencies.txt")
        except FileNotFoundError:
            print("File \'frequencies.txt\' not found! create it using identify_frequencies, with save=True")

        return freq

    def find_source_frequency(self, data, minimaorder):
        """
        data: The decorrelated data, with the source being the lowest point in the observation

        returns the frequency at the brightest point of the source observation(lowest)
        """
        f_bright = []
        # find the lowest point.

        # TODO: Make this parallel

        for i in range(len(data)):
            l1 = data[i]
            minima_index, minima = self.minima_func(l1, minimaorder)
            #     print(minima_index, minima)
            f_brightest = heapq.nsmallest(1, minima)[-1]
            f_bright.append([i + 1, f_brightest])
        print(f_bright)
        #find the matching pixel id.
        pix_bright = []
        index = []

        for idMKID in self.idMKIDs_use:
            if idMKID not in self.mkid_exclude:
                index.append(idMKID)

        for i in range(len(f_bright)):
            pix_bright.append([index[i], f_bright[i][1]])

        return pix_bright

    def antenna_temperature(self, frequency_table, pix_fbright, t_atm, plot=False):
        """
        frequency_table: the frequency file generated using "identify_frequencies"
        frequency_brightest: ferq and f_brightest created using find_source_frequency
        t_atm: the atmospheric temperature that day

        returns the pixel-id and antenna temperature.
        """
        t_a = []
        f_load_av = []

        for k in range(len(pix_fbright)):
            index = pix_fbright[k][0]
            #     print(index)
            j = int(index - 1)
            f_load = np.mean([frequency_table[j][-1], frequency_table[j][-2]])
            f_load_av.append(f_load)
            if f_load != 0:
                print(j, f_load)
                print(k)
                t_a.append([index, -t_atm * (pix_fbright[k][1] / f_load)])

        t_a = np.array(t_a)

        if plot:
            plt.rcParams['figure.figsize'] = [10, 6]

            plt.title("T_a* plot")
            plt.plot(t_a[:, 0], t_a[:, 1], "o", c='r')
            plt.hlines(np.average(t_a[:, 1]), xmin=0, xmax=np.amax(t_a[:, 0]), label="average :{: .2f} K".format(np.average(t_a[:, 1])))
            plt.xlabel("pixels")
            plt.ylabel("Antenna temperature (K)")
            plt.legend()
            plt.show()
        return t_a



    def _planck(self, nu, T):
        h = 6.626e-34
        # c = 3.0e+8
        k = 1.38e-23

        a = (h * nu) / k
        b = (np.exp((h * nu) / (k * T)) - 1)
        intensity = a / b
        return intensity

    def _main_beam_efficiency(self, ta_star, nu, t_source, theta_eq, theta_pol, theta_mb):

        #     mb_eff = Ta_star / ( T_source * (1 - np.exp(-np.log(2)*((theta_eq * theta_pol)/theta_mb))))

        mb_eff = (0.5 * ta_star) / ((self._planck(nu, t_source) - self._planck(nu, 2.28)) * (
                    1 - np.exp(-np.log(2) * ((theta_eq * theta_pol) / theta_mb ** 2))))

        return mb_eff

    def get_main_beam_efficiency(self, ta_star, nu, t_source, theta_eq, theta_pol, theta_mb, plot=False):
        eta_mb = []

        for i in range(len(ta_star)):
            eta_mb.append([ta_star[i][0], self._main_beam_efficiency(ta_star[i][1], nu, t_source, theta_eq, theta_pol, theta_mb)])

        eta_mb = np.array(eta_mb)
        eta_mb[:, 1] = eta_mb[:, 1] * 100

        if plot:
            plt.rcParams['figure.figsize'] = [10, 6]

            plt.title("$\eta_{\mathrm{mb}}$ plot")
            plt.plot(eta_mb[:, 0],eta_mb[:, 1], "o")
            plt.hlines(np.average(eta_mb[:, 1]), xmin=0, xmax=np.amax(eta_mb[:, 0]), label="average :{: .2f} %".format(np.average(eta_mb[:, 1])))
            plt.xlabel("pixels")
            plt.ylabel("Main Beam efficiency")
            plt.legend()
            plt.show()

        return eta_mb