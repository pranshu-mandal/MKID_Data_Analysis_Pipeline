"""
To describe the workflow using some methods

** Intended to be used line by line
"""
import matplotlib.pyplot as plt # to be removed later.
import numpy.ma as ma
import numpy as np

import data_reduction as dr

bc = dr.BinaryConvert()  # Call the object
# mkid_files = bc.select_files() # Select the MKID readout files
# bc.save_fsp_pickle(mkid_files, "Mars_0601.pickle") # to average and save the fsp files to a single binary file.

# data = bc.load_pickle("Mars_0601.pickle")

#### To check which pixels are good #####
# for i in range(len(data)):
#     bc.plot(data, i+1)
#########################################

# exclude = [48, 71, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 89, 90, 91, 92, 93, 94, 95, 98, 99]
# bc.pixels_to_exclude(exclude) # to exclude these pixels from further calculations.

# freq = bc.identify_frequencies(data, 4000, exclude=True, plot=False, save_file=True)

# clipped = bc.clip(data, 16 * 255, -16 * 160) # clip for decorrelation
#
# jump_removed = bc.changepoint_remove(clipped, 1000, 17173)
#
# mean_sub = bc.mean_center(clipped)

# print(len(mean_sub))  # Important, use only the good and mean subtracted data.
# plt.plot(mean_sub[2], ',')
# plt.show()

# res_pca = bc.pca_decor(mean_sub, n_components=3)

# for i in range(10):
#     plt.plot(res_pca[i][:, 1], alpha = 0.5)

# plt.plot(results[0][:,2], alpha=0.5)
# plt.plot(results[0][:,1])
# plt.show()

# plt.plot(results[0][:,0])
# plt.show()

# flat = bc.flatten(res_pca, 0)
# mask = bc.get_sigmaclipped_mask(flat, avg_points=32, sigma=4, maxiters=4, axis=1)


# ----------- TO GET THE SIGMA CLIPPED MASK PLOT --------------------
# pix = 4
# # plt.plot(running_avged[pix])
# plt.plot(ma.masked_array(mean_sub[pix], ~mask[pix]), 'o', label='Masked data')
# plt.plot(ma.masked_array(mean_sub[pix], mask[pix]), label='Unmasked data', alpha = 0.5)
# plt.legend()
# plt.title('Sigma clipping at 3-sigma')
# plt.show()

# -------------------------------------------------------------------

# chunk_matrix = bc.make_chunk_matrix(mask, num_chunks=400)

# plt.imshow(chunk_matrix)
# plt.colorbar()
# plt.show()

# ---------  CHUNK PCA METHOD CHECKING ------------------------------
#
# final = bc.chunkpca_decor(mask, chunk_matrix, mean_sub)
#
# plt.rcParams['figure.figsize'] = [16, 8]
# N = 32
# pix = 4
# plt.plot(final[pix][:,2])# - np.mean(final[pix][:,0]))
# plt.plot(final[pix][:,1])
# plt.plot(ma.masked_array(mean_sub[pix][(0):(70000)], mask[pix][(0):(70000)]), alpha=0.6)
# plt.show()
#
# plt.plot(np.convolve(final[pix][:,0], np.ones((N,))/N, mode='same'))
# plt.show()
#
# final_flat = bc.flatten(final, 0)
# bc.save_as_pickle(final_flat, "final_flat")

cleaned = bc.load_pickle("final_flat.pickle")

plt.plot(cleaned[0], )
plt.show()

freq = np.loadtxt("frequencies.txt")