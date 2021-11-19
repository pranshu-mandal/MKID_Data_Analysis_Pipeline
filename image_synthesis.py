import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import cygrid

# self made python files
import param
import imaging

Param = param.Param.copy()


class SingleMKIDmapper():

    def map(self, data, pixel, settings):
        """meant to perform a check on the settings parameters, and if anything is found wrong, popup an error message"""
        if settings['path'] == None:
            print('ERROR: Path to the folder is not specified')
        if settings['position_log'] == None:
            print('ERROR: position_log is not specified')
        if settings['mkids_to_exclude'] == None:
            print('No MKIDs are being excluded')
        self.work_fits(data, settings)

    def work_fits(self, data, settings):
        """
        redirects the data to the respective functions.
        :param data: the 1D TOD of the given pixel
        :param settings: has other arguments like antenna log and stuff.
        :return: single pixel cygrid datacube.
        """
        path_to_folder = settings["path"]
        # filename = settings["filename"]
        position_log = path_to_folder + settings["position_log"]
        mkids_exclude = settings["mkids_to_exclude"]  # mkids to exclude
        mkid_end = settings["mkids_end"]

        return self.run(f=data, position_log=position_log, settings=settings)

    def run(self, f, position_log, settings):
        """
        adjusts the cropping of the two values. calls imaging intensity, and calls write_map function.
        :param f: the 1D tod from the work_fits.
        :param position_log: the antenna log file name
        :param settings: has other information regarding the scan.
        :return: calls the write_map function.
        """

        # Cropping the transition scan values
        manual_selection = settings["position_crop"]  # the indexes of the position log cropping
        offset = settings["intensity_offset"]  # the value to allign the position log with the correct position log.
        maptype = settings["map_type"]
        crop = []
        sampling_rate = settings["sampling_rate"]

        if settings["scan_size"] == "auto":
            print("scan is auto")
            crop = 0
        elif len(settings["scan_size"]) == 2:
            crop = settings["scan_size"]
        else:
            print(('\x1b[0;31;40m'
                   + 'scan size must be a list of 2 parameters, X-scanwidth, and Y-scanheight in arcsec Could not crop' +
                   '\n Uncropeed data imported'
                   '\x1b[0m'))

        print(crop, 'crop1')

        ra, dec, index_cropped = imaging.position_import(antlog=position_log, maptype=maptype,
                                                         manual_selection=manual_selection,
                                                         crop=crop, sampling_rate=sampling_rate)
        intensity = imaging.intensity_import(f, index_cropped=index_cropped, offset=offset, normalize=True)


        # if settings['plot_type'] == 'frequency_shift':
        #     pass
        # elif settings['plot_type'] == 'tmb':
        #     intensity = imaging.antenna_temp(intensity, t_atm=settings['t_atm'])



        if settings['show_scatter_plot'] == True:
            plt.scatter(ra, dec, c=intensity)
            plt.title('')
            plt.show()
        else:
            pass

        return self.write_map(intensity, ra, dec, settings)

    def write_map(self, intensity, ra, dec, settings):
        """
        uses all the information provided for intensity, and positions it returns the data cube from cygrid.
        :param intensity: the cropped 1D TOD for the given pixel.
        :param ra: ra position
        :param dec: dec position.
        :param settings: other information.
        :return: the cygrid data cube.
        """
        overall_maxra = np.amax(ra)
        overall_minra = np.amin(ra)
        overall_maxdec = np.amax(dec)
        overall_mindec = np.amin(dec)
        map_width, map_height = (overall_maxra - overall_minra), (
                overall_maxdec - overall_mindec)  # in RA-dec

        #         map_width, map_height = (overall_maxra - overall_minra) * 3600., (
        #                 overall_maxdec - overall_mindec) * 3600.  # in arcsecs

        # print map_width, map_height

        num_mkids = 101  # len(timewithra[0, :]) - 1

        if settings["beamsize_fwhm"] == None:
            beamsize_fwhm = 16.  # arcsec;
        else:
            beamsize_fwhm = settings["beamsize_fwhm"]

        '''IF THE gridder.grid() gives nan as output from the intensity, the BEAMSIZE / KERNELSIZE is too small
         in the calculation.'''

        target_header = imaging.setup_header((((overall_minra + (overall_maxra - overall_minra) / 2)),
                                              (overall_mindec + (overall_maxdec - overall_mindec) / 2)),
                                             ((overall_maxdec - overall_mindec), (overall_maxdec - overall_mindec)),
                                             beamsize_fwhm / 3600., )  # (overall_maxra - overall_minra) replaced with height to account for removing the transition positioning
        target_header['NAXIS3'] = 1  # dummy spectral axis

        target_wcs = WCS(target_header)

        # print 'Header:', target_header

        # Defining the gridder

        gridder = cygrid.WcsGrid(target_header)

        kernelsize_fwhm = 4.
        kernelsize_fwhm /= 3600.  # need to convert to degree
        kernelsize_sigma = kernelsize_fwhm / np.sqrt(8 * np.log(2))
        support_radius = 4. * kernelsize_sigma
        healpix_reso = kernelsize_sigma / 2.

        gridder.set_kernel(
            'gauss1d',
            (kernelsize_sigma,),
            support_radius,
            healpix_reso,
        )

        tmp_map = np.array(intensity)
        tmp_map_f = tmp_map.flatten()
        xcoords = np.array(ra)
        ycoords = np.array(dec)

        xcoords = np.array(xcoords)
        ycoords = np.array(ycoords)

        gridder.grid(xcoords.flatten(), ycoords.flatten(), tmp_map_f[:, None])

        cygrid_cube = gridder.get_datacube().squeeze()
        target_wcs = gridder.get_wcs()

        return cygrid_cube, target_header


class BeamCharcteristics():
    """
    This is similar to the SIngleMKIDMapper, but instead of returning the cygrid datacube, it returns the center positions for a given mkid.
    """

    def map(self, data, settings):
        """meant to perform a check on the settings parameters, and if anything is found wrong, popup an error message"""
        if settings['path'] == None:
            print('ERROR: Path to the folder is not specified')
        if settings['position_log'] == None:
            print('ERROR: position_log is not specified')
        if settings['mkids_to_exclude'] == None:
            print('No MKIDs are being excluded')
        self.work_fits(data, settings)

    def work_fits(self, data, settings):
        """
        redirects the data to the respective functions.
        :param data: the 1D TOD of the given pixel
        :param settings: has other arguments like antenna log and stuff.
        :return: single pixel cygrid datacube.
        """
        path_to_folder = settings["path"]
        # filename = settings["filename"]
        position_log = path_to_folder + settings["position_log"]
        mkids_exclude = settings["mkids_to_exclude"]  # mkids to exclude
        mkid_end = settings["mkids_end"]

        return self.run(f=data, position_log=position_log, settings=settings)

    def run(self, f, position_log, settings):
        """
        adjusts the cropping of the two values. calls imaging intensity, and calls write_map function.
        :param f: the 1D tod from the work_fits.
        :param position_log: the antenna log file name
        :param settings: has other information regarding the scan.
        :return: calls the write_map function.
        """

        # Cropping the transition scan values
        manual_selection = settings["position_crop"]  # the indexes of the position log cropping
        offset = settings["intensity_offset"]  # the value to allign the position log with the correct position log.
        crop = []

        if settings["scan_size"] == "auto":
            print("scan is auto")
            crop = 0
        elif len(settings["scan_size"]) == 2:
            crop = settings["scan_size"]
        else:
            print(('\x1b[0;31;40m'
                   + 'scan size must be a list of 2 parameters, X-scanwidth, and Y-scanheight in arcsec Could not crop' +
                   '\n Uncropeed data imported'
                   '\x1b[0m'))

        print(crop, 'crop1')

        ra, dec, index_cropped = imaging.position_import(antlog=position_log, maptype=settings['map_type'], manual_selection=manual_selection,
                                                         crop=crop)

        intensity = imaging.intensity_import(f, index_cropped=index_cropped, offset=offset, normalize=True)

        if settings['show_scatter_plot'] == True:
            plt.scatter(ra, dec, c=intensity)
            plt.title('')
            plt.show()
        else:
            pass

        return self.write_map(intensity, ra, dec, settings)

    def write_map(self, intensity, ra, dec, settings):
        """
        uses all the information provided for intensity, and positions it returns the data cube from cygrid.
        :param intensity: the cropped 1D TOD for the given pixel.
        :param ra: ra position
        :param dec: dec position.
        :param settings: other information.
        :return: the cygrid data cube.
        """
        overall_maxra = np.amax(ra)
        overall_minra = np.amin(ra)
        overall_maxdec = np.amax(dec)
        overall_mindec = np.amin(dec)
        map_width, map_height = (overall_maxra - overall_minra), (
                overall_maxdec - overall_mindec)  # in RA-dec

        #         map_width, map_height = (overall_maxra - overall_minra) * 3600., (
        #                 overall_maxdec - overall_mindec) * 3600.  # in arcsecs

        # print map_width, map_height

        num_mkids = 101  # len(timewithra[0, :]) - 1

        if settings["beamsize_fwhm"] == None:
            beamsize_fwhm = 16.  # arcsec;
        else:
            beamsize_fwhm = settings["beamsize_fwhm"]

        '''IF THE gridder.grid() gives nan as output from the intensity, the BEAMSIZE / KERNELSIZE is too small
         in the calculation.'''

        target_header = imaging.setup_header((((overall_minra + (overall_maxra - overall_minra) / 2)),
                                              (overall_mindec + (overall_maxdec - overall_mindec) / 2)),
                                             ((overall_maxdec - overall_mindec), (overall_maxdec - overall_mindec)),
                                             beamsize_fwhm / 3600., )  # (overall_maxra - overall_minra) replaced with height to account for removing the transition positioning
        target_header['NAXIS3'] = 1  # dummy spectral axis

        target_wcs = WCS(target_header)

        # print 'Header:', target_header

        # Defining the gridder

        gridder = cygrid.WcsGrid(target_header)

        kernelsize_fwhm = 4.
        kernelsize_fwhm /= 3600.  # need to convert to degree
        kernelsize_sigma = kernelsize_fwhm / np.sqrt(8 * np.log(2))
        support_radius = 4. * kernelsize_sigma
        healpix_reso = kernelsize_sigma / 2.

        gridder.set_kernel(
            'gauss1d',
            (kernelsize_sigma,),
            support_radius,
            healpix_reso,
        )

        tmp_map = np.array(intensity)
        tmp_map_f = tmp_map.flatten()
        xcoords = np.array(ra)
        ycoords = np.array(dec)

        xcoords = np.array(xcoords)
        ycoords = np.array(ycoords)

        gridder.grid(xcoords.flatten(), ycoords.flatten(), tmp_map_f[:, None])

        cygrid_cube = gridder.get_datacube().squeeze()

        # Get the db conversion and then plot the contour

        center = []

        cygrid_cube = 10 * np.log10(cygrid_cube)  # dB conversion

        xaxmin = -120
        xaxmax = 120
        yaxmin = -120
        yaxmax = 120

        extent = (xaxmin, xaxmax, yaxmin, yaxmax)

        levels = [-3]  # defining the contour level
        plt.imshow(cygrid_cube)
        CS = plt.contour(cygrid_cube, levels, colors=('r', (1, 1, 0), '#afeeee', '0.5'), extent=extent)  # , extent=extent)
        dat0 = CS.allsegs[0][0]
        CS_center = [np.mean(dat0[:, 0]), np.mean(dat0[:, 1])]
        CS_beamsize = [np.amax(dat0[:, 0]) - np.amin(dat0[:, 0]), np.amax(dat0[:, 1]) - np.amin(dat0[:, 1])]
        CS_beamsize_average = np.mean(CS_beamsize)
        # CS_plot = plt.text(np.mean(dat0[:, 0]), np.mean(dat0[:, 1]), "mkid {}, bw = {}".format(mkid_exclude[i] , CS_beamsize_average), size=6)
        center.append([CS_center[0], CS_center[1], CS_beamsize_average])
        # plt.show()
        plt.close()

        print("BEAM CENTER : ", center)

        w = WCS(target_header, naxis=2)
        world = w.wcs_pix2world(np.array([[center[0][0], center[0][1]]], dtype=np.float64), 0)


        return CS_beamsize_average

        i

        #
        # target_wcs = gridder.get_wcs()

        # return cygrid_cube, target_wcs

class CompositeMKIDmapper():

    def __init__(self):
        self.freq=[]
        self.roomtemp=[]
        self.etamb = []


    def map(self, data, mkid_ids, offsets, settings):
        """meant to perform a check on the settings parameters, and if anything is found wrong, popup an error message"""
        if settings['path'] == None:
            print('ERROR: Path to the folder is not specified')
        if settings['position_log'] == None:
            print('ERROR: position_log is not specified')
        if settings['mkids_to_exclude'] == None:
            print('No MKIDs are being excluded')
        self.work_fits(data, settings)

    def work_fits(self, data, frequencies, eta_mb, mkid_ids, offsets, settings, mkid_end, roomtemp=282):
        """
        works with the offsets table to properly adjust the coordinates of the map.
        :param data: the 2D cleaned data matrix.
        :param mkid_ids: the list of MKID ids that are there in the 2d data in the same order.
        :param offsets: the offset table which contains the MKID id, offset_x, offset_y relative to the central pixel.
        :param settings: cotains other information like beam width, position log filename etc.
        :return: the composite cygrid data cube.
        """
        self.freq = frequencies
        self.roomtemp = roomtemp
        self.etamb = eta_mb




        path_to_folder = settings["path"]
        filename = settings["filename"]
        position_log = path_to_folder + settings["position_log"]

        three_values = []

        for i in range(mkid_end): #len(offsets)
            for j, val in enumerate(mkid_ids):
                if offsets[i][0] == val:
                    three_values.append(self.run(j, f=data[j], off_x = offsets[i][1], off_y = offsets[i][2],
                                                 position_log=position_log, settings=settings))


        three_values = np.array(three_values)
        intensity_composite = three_values[:, 0]
        ra_composite = three_values[:, 1]
        dec_composite = three_values[:, 2]

        if settings['show_scatter_plot'] == True:
            print(" You should be careful scatter plotting so much points in python. It may crash, if you still want it, change the source code")
            # plt.scatter(ra_composite, dec_composite, c=intensity_composite)
            # plt.title('')
            # plt.show()
        else:
            pass

        return self.composite_write_map(intensity_composite, ra_composite, dec_composite, settings)


    def run(self, j, f, off_x, off_y, position_log, settings):

        # Cropping the transition scan values
        manual_selection = settings["position_crop"]
        intensity_offset = settings["intensity_offset"]
        map_type = settings["map_type"]
        crop = []

        if settings["scan_size"] == "auto":
            print("scan is auto")
            crop = 0
        elif len(settings["scan_size"]) == 2:
            crop = settings["scan_size"]
        else:
            print(('\x1b[0;31;40m'
                   + 'scan size must be a list of 2 parameters, X-scanwidth, and Y-scanheight in arcsec Could not crop' +
                   '\n Uncropeed data imported'
                   '\x1b[0m'))

        print(crop, 'crop1')

#         ra, dec, index_cropped = imaging.position_import(antlog=position_log, manual_selection=manual_selection,
#                                 crop=crop)

        ra, dec, index_cropped = imaging.pos_import_offset(antlog=position_log, maptype= map_type, manual_selection=manual_selection,  crop=crop, offset_x=off_x, offset_y=off_y)
# #         position_offset = 3 #position_offset
#         ra = ra - (off_x/3600.)
#         dec = dec - (off_y/3600.)

        intensity = np.array(imaging.intensity_import(f, index_cropped=index_cropped, offset=intensity_offset, normalize=False))

        if settings['plot_type'] == 'tmb':
            t_a = []
            f_load_av = []

            # for k in range(len(self.etamb)):

            index = self.etamb[j][0]
            
            #     print(index)
            l = int(index - 1)
            f_load = np.mean([self.freq[l][-1], self.freq[l][-2]])
            f_load_av.append(f_load)
            if f_load != 0:
                # print(j, f_load)
                # print(k)
                t_a.append(-self.roomtemp * (-intensity / f_load))

            t_a = np.array(t_a)
            plt.close()
            plt.plot(t_a)
            plt.show()
            t_mb = t_a/(self.etamb[j][1]/100)
            intensity= t_mb

        ######################### Changed For Tastar calculation ################################
        # f_load = 7.5
        # t_atm = 300 #K
        # intensity = t_atm * (intensity/f_load)
        ######################### Changed For Tastar calculation ################################


        # plt.scatter(ra, dec, c=intensity)
        # plt.title('')
        # plt.show()

        return intensity.flatten(), ra, dec

    def composite_write_map(self, intensity, ra, dec, settings):
        """
                uses all the information provided for intensity, and positions it returns the data cube from cygrid.
                :param intensity: the cropped 1D TOD for the given pixel.
                :param ra: ra position
                :param dec: dec position.
                :param settings: other information.
                :return: the cygrid data cube.
                """
        overall_maxra = np.amax(ra)
        overall_minra = np.amin(ra)
        overall_maxdec = np.amax(dec)
        overall_mindec = np.amin(dec)
        map_width, map_height = (overall_maxra - overall_minra) * 3600., (
                overall_maxdec - overall_mindec) * 3600.  # in arcsecs

        # print map_width, map_height

        num_mkids = 101  # len(timewithra[0, :]) - 1

        if settings["beamsize_fwhm"] == None:
            beamsize_fwhm = 16.  # arcsec;
        else:
            beamsize_fwhm = settings["beamsize_fwhm"]

        '''IF THE gridder.grid() gives nan as output from the intensity, the BEAMSIZE / KERNELSIZE is too small
         in the calculation.'''

        target_header = imaging.setup_header((((overall_minra + (overall_maxra - overall_minra) / 2)),
                                              (overall_mindec + (overall_maxdec - overall_mindec) / 2)),
                                             ((overall_maxdec - overall_mindec), (overall_maxdec - overall_mindec)),
                                             beamsize_fwhm / 3600., )  # (overall_maxra - overall_minra) replaced with height to account for removing the transition positioning
        target_header['NAXIS3'] = 1  # dummy spectral axis

        target_wcs = WCS(target_header)

        # print 'Header:', target_header

        # Defining the gridder

        gridder = cygrid.WcsGrid(target_header)

        kernelsize_fwhm = 4.
        kernelsize_fwhm /= 3600.  # need to convert to degree
        kernelsize_sigma = kernelsize_fwhm / np.sqrt(8 * np.log(2))
        support_radius = 4. * kernelsize_sigma
        healpix_reso = kernelsize_sigma / 2.

        gridder.set_kernel(
            'gauss1d',
            (kernelsize_sigma,),
            support_radius,
            healpix_reso,
        )

        tmp_map = np.array(intensity)
        tmp_map_f = tmp_map.flatten()
        xcoords = np.array(ra)
        ycoords = np.array(dec)

        xcoords = np.array(xcoords)
        ycoords = np.array(ycoords)

        gridder.grid(xcoords.flatten(), ycoords.flatten(), tmp_map_f[:, None])

        cygrid_cube = gridder.get_datacube().squeeze()
        target_wcs = gridder.get_wcs()

        return cygrid_cube, target_header
        