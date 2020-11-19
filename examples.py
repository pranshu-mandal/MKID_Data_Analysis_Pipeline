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

bc = dr.BinaryConvert()  # Call the object
# mkid_files = bc.select_files() # Select the MKID readout files
# bc.save_fsp_pickle(mkid_files, "Mars_0601") # to average and save the fsp files to a single binary file.

data = bc.load_pickle("Mars_0601.pickle")

#### To check which pixels are good #####
# for i in range(len(data)):
#     bc.plot(data, i+1)
#########################################

exclude = [48, 71, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 89, 90, 91, 92, 93, 94, 95, 98, 99]
# exclude = [7, 13, 16, 38, 48, 61, 71, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 89, 90, 91, 92, 93, 94, 95, 97, 98, 99]
bc.pixels_to_exclude(exclude) # to exclude these pixels from further calculations.

# freq = bc.identify_frequencies(data, 4000, exclude=True, plot=False, save_file=True)
#
# clipped = bc.clip(data, 16 * 255, -16 * 160) # clip for decorrelation
#
# jump_removed = bc.changepoint_remove(clipped, 1000, 17173)
#
# mean_sub = bc.mean_center(clipped)
#
# print(len(mean_sub))  # Important, use only the good and mean subtracted data.
# plt.plot(mean_sub[2], ',')
# plt.show()
#
# res_pca = bc.pca_decor(mean_sub, n_components=3)
#
# # for i in range(10):
# #     plt.plot(res_pca[i][:, 1], alpha = 0.5)
#
# # plt.plot(results[0][:,2], alpha=0.5)
# # plt.plot(results[0][:,1])
# # plt.show()
#
# # plt.plot(results[0][:,0])
# # plt.show()
#
# flat = bc.flatten(res_pca, 0)
# mask = bc.get_sigmaclipped_mask(flat, avg_points=32, sigma=4, maxiters=4, axis=1)
#
#
# # ----------- TO GET THE SIGMA CLIPPED MASK PLOT --------------------
# # pix = 4
# # # plt.plot(running_avged[pix])
# # plt.plot(ma.masked_array(mean_sub[pix], ~mask[pix]), 'o', label='Masked data')
# # plt.plot(ma.masked_array(mean_sub[pix], mask[pix]), label='Unmasked data', alpha = 0.5)
# # plt.legend()
# # plt.title('Sigma clipping at 3-sigma')
# # plt.show()
#
# # -------------------------------------------------------------------
#
# chunk_matrix = bc.make_chunk_matrix(mask, num_chunks=400)
#
# # plt.imshow(chunk_matrix)
# # plt.colorbar()
# # plt.show()
#
# # ---------  CHUNK PCA METHOD CHECKING ------------------------------
# #
# final = bc.chunkpca_decor(mask, chunk_matrix, mean_sub)
# #
# # plt.rcParams['figure.figsize'] = [16, 8]
# # N = 32
# # pix = 4
# # plt.plot(final[pix][:,2])# - np.mean(final[pix][:,0]))
# # plt.plot(final[pix][:,1])
# # plt.plot(ma.masked_array(mean_sub[pix][(0):(70000)], mask[pix][(0):(70000)]), alpha=0.6)
# # plt.show()
# #
# # plt.plot(np.convolve(final[pix][:,0], np.ones((N,))/N, mode='same'))
# # plt.show()
# #
# final_flat = bc.flatten(final, 0)
# bc.save_as_pickle(final_flat, "final_flat_selected_by_beamsize")
#
cleaned = bc.load_pickle("final_flat.pickle")
# # ---------------- Calibration ---------------
# # # plt.plot(cleaned[0], )
# # # plt.show()
# #
# freq = bc.load_frequency()
# pix_bright = bc.find_source_frequency(cleaned, minimaorder=4000)
# np.savetxt('pix_bright.txt', pix_bright)
#
#
# antenna_temp = bc.antenna_temperature(freq, pix_bright, 282, plot=True)
# print(antenna_temp)
#
# eta_mb = bc.get_main_beam_efficiency(antenna_temp, 1e+11, 195, 15.3, 15.3, 17.7, plot=True)
#
# eta_mb = bc.get_main_beam_efficiency(antenna_temp, 1e+11, 195, 15.3, 15.3, 17, plot=True)

mkids_used = bc.mkids_used()


path_to_data = '/home/pranshu/Downloads/Data/rNROq_20180522093543_3600.MKID/data/'

""" So let's include a bit more user freindly way to get to the controls inside the mapping methods
Things changed from marsmap to venus map are
path_to_data:
filename:
mkids_to_exclude
mkids_end (maybe add an mkids start)
"""


settings = {
    'path': '/home/pranshu/git/Data_handling/',
    'filename': "rNROq_20180522093543_3600.MKID001.f0o8w5_0_29061.avrn8.txt",
    'position_log': "20180601054644",
    'mkids_to_exclude': exclude,
    'map_type': 'skymap', # 'skymap' or 'beammap'
    'mkids_end': 101,
    'beamsize_fwhm': 16, # In arcsecs.     (Default: 16 ")
    'kernelsize_fwhm': 4, # In arcsecs.
    'save_directory': "/home/pranshu/Downloads/Data/rNROq_20180522093543_3600.MKID/data/",
    'scan_size': [4, 4],#'auto', #[4, 4], # auto option is there. ----can be explicitely specified as [2, 2] for 2'X2' scan.
    #This is exclusively for Temperature callibratio using chopper-wheel method
    'chopper_wheel_frequency_table': 'frequencies.txt', # specify None value if not being used
    'Process' : 'beam_center_and_size', # 'combined_map'
    'position_crop':[64*64, -5*64], #The numbers are the beginning and end of the manual crop
    'intensity_offset': -456, #(-12*64)#(19*5)+1 #The intensity offset value to position the peaks correctly.212-v-347
    'show_scatter_plot': False, #scatter plot with the
    'mkids_used': mkids_used,
    'plot_type': 'tmb', # 'ta_star', 'freq_shift', tmb
    'return_type': 'beamcenter'
}

beam = im.BeamCharcteristics()

# ------------------- fIND bEAM SIZE ---------------------
# beamsize = []
# for i in range(len(cleaned)):
#     world_coord = beam.work_fits(data = (cleaned[i]), settings=settings)
#     beamsize.append([mkids_used[i], world_coord])
# beamsize = np.array(beamsize)
# np.savetxt('beamsize.txt', beamsize)
# ----------------------------------------------------------------
beamsize = np.loadtxt('beamsize.txt')
print(beamsize)

# plt.plot(cleaned[15])
# plt.show()

# plt.plot(beamsize[:, 0], beamsize[:, 1], 'o')
# plt.show()
# selecting the good beamsizes
beamaverage = np.mean(beamsize[:, 1])

new_sizes = []
threshhold = 3
for i in range(len(beamsize)):
    if beamsize[i][1] < beamaverage + threshhold:
        if beamsize[i][1] > beamaverage - threshhold:
            new_sizes.append([beamsize[i][0], beamsize[i][1]])
print(len(new_sizes))
new_sizes = np.array(new_sizes)
# plt.plot(new_sizes[:, 0], new_sizes[:, 1], 'o', c='green')
# plt.hlines(np.mean(new_sizes[:, 1]), 0, np.amax(new_sizes[:, 0]), label='average: {: .2f} \"' .format(np.mean(new_sizes[:, 1])))
# plt.title('Beamsize of elements of the camera')
# plt.xlabel('elements')
# plt.ylabel('beamsize (arcsec)')
# plt.ylim(0, 20)
# plt.legend()
# plt.show()

filtered_elements = []
for i in range(len(beamsize)):
    if beamsize[i][0] not in new_sizes[:, 0]:
        filtered_elements.append(beamsize[i][0])

print(filtered_elements)

freq = bc.load_frequency()
pix_bright = np.loadtxt('pix_bright.txt')
# antenna_temp = bc.antenna_temperature(freq, pix_bright, 282, exclude=mkids_used, plot=False)
# np.savetxt('beamsize_filtered.txt', new_sizes)
# eta_mb = bc.get_main_beam_efficiency(antenna_temp, 1e+11, 195, 15.3, 15.3, new_sizes, plot=False)

antenna_temp = bc.antenna_temperature(freq, pix_bright, 282, exclude=[], plot=False)
# np.savetxt('beamsize_filtered.txt', new_sizes)
eta_mb = bc.get_main_beam_efficiency(antenna_temp, 1e+11, 195, 15.3, 15.3, beamsize, plot=False)

# tmb = bc.aperture_efficiency(eta_mb, new_sizes,plot=True)


# ------------------------ Composite RA-Dec map in Tmb scale ------------------
marsmap = im.CompositeMKIDmapper()

offsets = np.loadtxt("RA_DEC_center_offsets.txt") #"center_offsets.txt" or "RA_DEC_center_offsets.txt"
cygrid_cube, target_header = marsmap.work_fits(data = cleaned, frequencies=freq, eta_mb=eta_mb, mkid_ids=mkids_used, offsets=offsets, settings=settings, mkid_end=1, roomtemp=282)
target_wcs = wcs.WCS(target_header)


axis_extent = 110

xaxmin = -axis_extent
xaxmax = axis_extent
yaxmin = -axis_extent
yaxmax = axis_extent

# print yaxmax
# print yaxmin


##########################################################

wmcut = target_wcs[40:130, 40:130]
datacut = cygrid_cube[40:130, 40:130]

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection=wmcut.celestial)
cax = fig.add_axes([0.905, 0.12, 0.03, 0.75])
im = ax.imshow(
    datacut,
    origin='lower', interpolation='nearest', vmin = 0
    )#, cmap='gist_ncar'
lon, lat = ax.coords
lon.set_axislabel('R.A. [deg]')
lat.set_axislabel('Dec [deg]')
ax.set_title('Map in $T_{\mathrm{mb}}$ (K) scale | Mars (06/01)')
fig.colorbar(im, cax=cax, orientation='vertical')
# CS = plt.contour(cygrid_db, levels, colors=('white', 'r', 'green', (1, 1, 0), '#afeeee', '0.5'))#, extent=extent)
# plt.savefig('./results/Mars_RA_DEC_projection_zoomed.pdf')
plt.show()





# marsmap = im.SingleMKIDmapper()
#
# for i in range(2):  # len(cleaned)
#     # i = 2
#
#     cygrid_cube, header = marsmap.work_fits(data=(cleaned[i]), pixel = mkids_used[i], settings=settings)
#
#     xaxmin = -120
#     xaxmax = 120
#     yaxmin = -120
#     yaxmax = 120
#
#     # print yaxmax
#     # print yaxmin
#
#     extent = (xaxmin, xaxmax, yaxmin, yaxmax)
#
#     levels = [-15, -10, -3]
#
#     cygrid_cube = np.array(cygrid_cube)
#
#     ### Correcting for the dB scale of negative numbers ######
#
#     #     for i in range(len(cygrid_cube)):
#     #         for j in range(len(cygrid_cube[i])):
#     #             if cygrid_cube[i][j] < 0:
#     #                 cygrid_cube[i][j] == 0.0001
#
#     #     ##########################################################
#
#     #     cygrid_cube= 10 * np.log10(cygrid_cube)
#
#     plt.figure(figsize=(10, 10))
#     plt.imshow(cygrid_cube, origin='lower', interpolation='nearest', cmap='gist_ncar', extent=extent, vmin=0.1, vmax=1)
#     plt.colorbar(fraction=0.046, pad=0.04)
#     CS = plt.contour(cygrid_cube, levels, colors=('white', 'r', 'green', (1, 1, 0), '#afeeee', '0.5'), extent=extent)
#     plt.clabel(CS, fontsize=9, inline=1)
#     #     ax.set_aspect('equal')
#     plt.title('Beammap-Mars, {}, Bw-{}", Kw-{}".'.format(mkids_used[i], 16, 4))
#     plt.ylabel('dEl (arcsec)')
#     plt.xlabel('dAz (arcsec)')
#     #     plt.savefig('./results/Mars_MKID_{}.pdf'.format(mkid_exclude[i]))
#     plt.show()
