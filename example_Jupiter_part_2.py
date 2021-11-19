import matplotlib.pyplot as plt # to be removed later.
import numpy.ma as ma
import numpy as np
from astropy import wcs


import data_reduction as dr
import image_synthesis as im

working_dir = "/home/pranshu/Downloads/work/Jupiter_2021/"

bc = dr.BinaryConvert()  # Call the object

# get the correct list of MKID_IDs -----------
original_data = bc.load_pickle(working_dir + 'Jupiter_0601.pickle')
exclude = [9, 12, 23, 25, 28, 35, 39, 41, 49, 58] # excluded pixels were selected using plotted data
bc.pixels_to_exclude(exclude)
mkids_used = bc.mkids_used()
print("MKIDs used :", mkids_used)
original_data = []
# --------------------------------------------

cleaned = bc.load_pickle(working_dir + "final_flat_selected_by_beamsize.pickle")
print(len(cleaned))
# plt.plot(cleaned[0])
# plt.show()



# Making the Beammaps ------------------------------
# Setting up the mapping parameters

path_to_data = '/home/pranshu/Downloads/work/Jupiter_2021/'

""" So let's include a bit more user freindly way to get to the controls inside the mapping methods
Things changed from marsmap to venus map are
path_to_data:
filename:
mkids_to_exclude
mkids_end (maybe add an mkids start)
"""

settings = {
    'path': '/home/pranshu/Downloads/work/Jupiter_2021/',
    'filename': "rNROq_20180522093543_3600.MKID001.f0o8w5_0_29061.avrn8.txt",
    'position_log': "20210608043832",
    'mkids_to_exclude': [21, 38, 39, 71, 77, 85, 89, 90, 91, 92],  # [21, 38, 39, 71, 77],
    'map_type': 'beammap',  # 'skymap' or 'beammap'
    'mkids_end': 4,
    'beamsize_fwhm': 16,  # In arcsecs.     (Default: 16 ")
    'kernelsize_fwhm': 4,  # In arcsecs.
    'save_directory': "/home/pranshu/Downloads/work/Jupiter_2021/",
    'scan_size': [6, 6],
    'sampling_rate': 256,
    # 'auto', #[4, 4], # auto option is there. ----can be explicitely specified as [2, 2] for 2'X2' scan.
    # This is exclusively for Temperature callibratio using chopper-wheel method
    'chopper_wheel_frequency_table': 'frequencies.txt',  # specify None value if not being used
    'Process': 'beam_center_and_size',  # 'combined_map'
    'position_crop': [52 * 256, -5 * 256],  # The numbers are the beginning and end of the manual crop
    'intensity_offset': -250,  # (-12*64)#(19*5)+1 #The intensity offset value to position the peaks correctly.212-v-347
    'show_scatter_plot': False,  # scatter plot with the
}

marsmap = im.SingleMKIDmapper()

# marsmap.map(data = results_pca_comp[0][:,0], settings=settings)
for i in range(len(cleaned)):  # len(final)
    # i = 2

    cygrid_cube, target_wcs = marsmap.work_fits(data=(cleaned[i]), settings=settings)

    xaxmin = -120
    xaxmax = 120
    yaxmin = -120
    yaxmax = 120

    # print yaxmax
    # print yaxmin

    extent = (xaxmin, xaxmax, yaxmin, yaxmax)

    levels = [-15, -10, -3]

    cygrid_cube = np.array(cygrid_cube)

    ### Correcting for the dB scale of negative numbers ######

    #     for i in range(len(cygrid_cube)):
    #         for j in range(len(cygrid_cube[i])):
    #             if cygrid_cube[i][j] < 0:
    #                 cygrid_cube[i][j] == 0.0001

    #     ##########################################################

    #     cygrid_cube= 10 * np.log10(cygrid_cube)

    plt.figure(figsize=(10, 10))
    plt.imshow(cygrid_cube, origin='lower', interpolation='nearest', cmap='gist_ncar', extent=extent, vmin=0.0, vmax=1)
    plt.colorbar(fraction=0.046, pad=0.04)
    CS = plt.contour(cygrid_cube, levels, colors=('white', 'r', 'green', (1, 1, 0), '#afeeee', '0.5'), extent=extent)
    plt.clabel(CS, fontsize=9, inline=1)
    #     ax.set_aspect('equal')
    plt.title('Beammap-Jupiter, {}, Bw-{}", Kw-{}".'.format(mkids_used[i], 16, 4))
    plt.ylabel('dEl (arcsec)')
    plt.xlabel('dAz (arcsec)')
    plt.savefig('/home/pranshu/Downloads/work/Jupiter_2021/results/beammaps/Jupiter_MKID_{}.png'.format(mkids_used[i]))
    # plt.show()




