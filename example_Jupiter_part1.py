"""
To describe the workflow using some methods

** Intended to be used line by line
"""
import matplotlib.pyplot as plt # to be removed later.
import numpy.ma as ma
import numpy as np
from astropy import wcs


import data_reduction as dr
import image_synthesis as im

working_dir = "/home/pranshu/Downloads/work/Jupiter_2021/"

bc = dr.BinaryConvert()  # Call the object
# mkid_files = bc.select_files() # Select the MKID readout files

# # # to average and save the fsp files to a single binary file.
# bc.save_fsp_pickle(mkid_files, working_dir + "Jupiter_0601", mode="single_col") #mode can be "average" of "single_col"

data = bc.load_pickle(working_dir + "Jupiter_0601.pickle")
# data = bc.load_pickle('Mars_0601.pickle')
data = np.array(data)
print(len(data))
print(data.shape)

sample_rate = 256
# for i in range(len(data)):
#     bc.plot(data, i+1)

exclude = [9, 12, 23, 25, 28, 35, 39, 41, 49, 58] # excluded pixels were selected using plotted data

bc.pixels_to_exclude(exclude) # to exclude these pixels from further calculations.

# freq = bc.identify_frequencies(data, 4000, exclude=True, plot=False, save_file=True, file_path=working_dir)

# IMP: set the clipping length properly
clipped = bc.clip(data, sample_rate * 70, -sample_rate * 75) # clip for decorrelation

mean_sub = bc.mean_center(clipped)

# print(len(mean_sub))  # Important, use only the good and mean subtracted data.
# plt.plot(mean_sub[5], ',')
# plt.show()

# for i in range(10):
#     plt.plot(mean_sub[i], alpha = 0.5)
#     plt.show()


res_pca = bc.pca_decor(mean_sub, n_components=3)

# for i in range(10):
#     plt.plot(res_pca[i][:, 1], alpha = 0.5)
#     plt.show()

flat = bc.flatten(res_pca, 0)
# IMP: reduce avg_points if the masked data is too sparse
mask = bc.get_sigmaclipped_mask(flat, avg_points=32, sigma=4, maxiters=4, axis=1)

# # ----------- TO GET THE SIGMA CLIPPED MASK PLOT --------------------
pix = 4
# plt.plot(running_avged[pix])
plt.plot(ma.masked_array(mean_sub[pix], ~mask[pix]), 'o', label='Masked data')
plt.plot(ma.masked_array(mean_sub[pix], mask[pix]), label='Unmasked data', alpha = 0.5)
plt.legend()
plt.title('Sigma clipping at 3-sigma')
plt.show()


# -------------------------------------------------------------------

chunk_matrix = bc.make_chunk_matrix(mask, num_chunks=5600)

# plt.imshow(chunk_matrix)
# plt.colorbar()
# plt.show()

# ---------  CHUNK PCA METHOD CHECKING ------------------------------
#
final = bc.chunkpca_decor(mask, chunk_matrix, mean_sub)
#
plt.rcParams['figure.figsize'] = [16, 8]
N = 32
pix = 4
plt.plot(final[pix][:,2])# - np.mean(final[pix][:,0]))
plt.plot(final[pix][:,1])
plt.plot(ma.masked_array(mean_sub[pix][(0):(70000)], mask[pix][(0):(70000)]), alpha=0.6)
plt.show()

plt.plot(np.convolve(final[pix][:,0], np.ones((N,))/N, mode='same'))
plt.show()
#
final_flat = bc.flatten(final, 0)
bc.save_as_pickle(final_flat, working_dir + "final_flat_selected_by_beamsize")

cleaned = bc.load_pickle(working_dir + "final_flat_selected_by_beamsize")

plt.plot(cleaned[0])
plt.show()
# # ---------------- Calibration ---------------
# # # plt.plot(cleaned[0], )
# # # plt.show()
