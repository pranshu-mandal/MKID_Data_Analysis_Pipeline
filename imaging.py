#!python
# -*- coding: utf-8 -*-

import numpy as np
import datetime
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.interpolate import interp1d


def pos_import_offset(antlog, maptype, manual_selection, crop, offset_x, offset_y):
    """
    Calculates several coordinate related values from the antenna log. Converted to python from the code written in
    C sent by Nagai-san.

    :param antlog: antenna log file from the telescope
    :return: RA-Dec, dAz-dEl, etc.

    Improvements to be made: Commandline selection of the output, maybe by an "if".
    Rotation due to change in elevation calculation.
    """

    LAT_NRO = 0.62729799791922308

    with open(antlog) as f:
        fdata = [line.rstrip() for line in f]  # this splits each row
        #	if "AZEL" in line:
        cdata = [line.split() for line in fdata]  # this splits each column
        a = 1  # this is the number of row
        #	b=2 # this is the number of column in the above row
        # c=len(fdata)
        c = len(fdata)
        # print len(fdata)
        delout = []
        dazout = []
        elout = []
        time = []

        # Observing in AZ-EL------------------
        realazelf0 = []
        realazelf1 = []
        progazelf0 = []
        progazelf1 = []

        # observing in RA-DEC-----------------
        realradecf0 = []
        realradecf1 = []
        realradecval0 = []
        realradecval1 = []

        #########   the start of loop for offsets
        # ofc = len(lens_offsetdeg_x)
        ofb = 0
        dazimuth = []
        dazcol = []

        while (a != (c - 1)):
            a += 1

            # print a
            #		g= fdata[a] #in case of extracting a single row
            n = cdata[a]  # this is splitting the a`th row into the 'n' array

            #	eval(str) converts a string to a number which can be calculated with
            adprog00 = (eval(n[1]) / 180) * np.pi  # RA value
            adprog01 = (eval(n[2]) / 180) * np.pi  # DEC value
            aeprog0 = (eval(n[3]) / 180) * np.pi  # AZ value
            aeprog1 = (eval(n[4]) / 180) * np.pi  # El value
            aeprog00 = (eval(n[9]) / 180) * np.pi  # AZ center indication value(not including offset)
            aeprog01 = (eval(n[10]) / 180) * np.pi  # El center indication value(not including offset)
            aereal0 = (eval(n[5]) / 180) * np.pi  # AZ encoder value(collimator value)
            aereal1 = (eval(n[6]) / 180) * np.pi  # El encoder value(collimator value)

            #	Beam offsets

            dofsx = ((0 / 3600.0) / 180) * np.pi  # x beam offset is -40"
            dofsy = ((0 / 3600.0) / 180) * np.pi  # y beam offset is 40"

            #	Rotaion by elevation

            rxval = -aereal1

            #	Scan offset + beam offset

            dEl = dofsx * np.sin(rxval) + dofsy * np.cos(rxval) + (aereal1 - aeprog01)
            dAz = dofsx * np.cos(rxval) - dofsy * np.sin(rxval) + ((aereal0 - aeprog00)) * np.cos(aereal1)

            dEl = dEl - ((580 / 3600.0) / 180.) * np.pi
            #
            dAz = dAz + ((40 / 3600.0) / 180.) * np.pi
            #	Add offsets

            aereal1 = aeprog01 + dEl
            aereal0 = aeprog00 + dAz / np.cos(aereal1)

            #	Calculate parallactic angle to convert dAZEL to dRADEC
            para = - np.arcsin(np.sin(aereal0) * np.cos(LAT_NRO) / np.cos(adprog01))
            dRa = -np.cos(para) * dAz + np.sin(para) * dEl
            dDec = np.sin(para) * dAz + np.cos(para) * dEl

            #	Add dRa DDec

            adreal1 = adprog01 + dDec
            adreal0 = adprog00 + dRa / np.cos(adreal1)

            #	converting back to proper units

            progazel0 = (aeprog0 / np.pi) * 180
            progazel1 = (aeprog1 / np.pi) * 180
            realazel0 = (aereal0 / np.pi) * 180
            realazel1 = (aereal1 / np.pi) * 180
            dAz = (dAz / np.pi) * 180  # * 3600
            dEl = (dEl / np.pi) * 180  # * 3600
            progradec0 = (adprog00 / np.pi) * 180
            progradec1 = (adprog01 / np.pi) * 180
            realradec0 = (adreal0 / np.pi) * 180
            realradec1 = (adreal1 / np.pi) * 180
            #		appending the calculated values into arrays
            delout.append(dEl)                 ########## Changing this for now
            dazout.append(dAz)
            # elout.append(realazel1)
            realazelf0.append(realazel0)
            realazelf1.append(realazel1)
            progazelf0.append(progazel0)
            progazelf1.append(progazel1)
            realradecf0.append(realradec0)
            realradecf1.append(realradec1)
            time1 = eval(n[0])  # extract time from log
            time1 = repr(time1)  # converts time1 from float to string which is essential for strptime()
            dt = datetime.datetime.strptime(time1, '%y%m%d%H%M%S.%f')
            time.append(dt.time())
            # time.append(time1)
            ofb += 1

        # print ("Antenna log starting time is {}".format(time[0]))
        # print ("Antenna log ending time is {}".format(time[len(time) - 1]))



        # for i in range(len(time)):
        #     print time[i].time()

        # print time,realradecf0, realradecf1

        #

        ################################# plot the antenna log #####################

    # print realazelf0


    # plt.xlabel('Az')
    # plt.ylabel('El')
    # plt.title('Venus observation Az-El position in the sky')
    # dazout = dazout * 3600
    # delout = delout * 3600
    # plt.scatter(dazout, delout)
    # plt.title('dEl-dAz with beam-offset (-60,60)')
    # plt.xlabel('dAz(deg)')
    # plt.ylabel('dEl(deg)')
    # # plt.scatter(realradecf0, realradecf1)
    # plt.show()

    if maptype == "beammap":
        realradecf0 = dazout                                      ########## CHANGED THIS  ! aND IT WORKED
        realradecf1 = delout
    elif maptype == "skymap":
        pass
    else:
        print("Please check *map_type* parameter in the settings")




    ############*OLD CODE*################# Define the starttime and average the position, and take time every second
    # start_time = datetime.datetime.strptime("05:49:00", "%H:%M:%S").time()
    #
    # # print ("Antenna log starting time is {}".format(time[0]))
    # # print ("Antenna log ending time is {}".format(time[len(time) - 1]))
    # begin_index = [i for i, j in enumerate(time) if j == start_time]
    # # print [j for i, j in enumerate(time) if j == start_time]
    # # print begin_index[0]

    # realradecf1 = realradecf1[begin_index[0]:]
    # realradecf0 = realradecf0[begin_index[0]:]
    # time = time[begin_index[0]:]

    avtime = time[::10]
    # print 'Antlog recorded ', len(avtime), 's'



    realradecf0 = realradecf0[190:-540]  # For MARSx4
    realradecf1 = realradecf1[190:-540]

    # realradecf0 = realradecf0[41:-38]  # FOr venus 5/22
    # realradecf1 = realradecf1[41:-38]

    # ## INTERPOLATING RA DEC SEPARATELY #######


    sampling_rate = 16*4
    b = len(realradecf0)
    rem = b % 5   # 5 samples for half a second. ( 10 per second)



    #trimming the end of the position log such that the interpolation can be done properly
    if rem !=0:
        realradecf0 = realradecf0[:-rem] # samples more than half a second are discarded
        realradecf1 = realradecf1[:-rem]


    # The length of ineterpolation
    t_obs = (len(realradecf0) / 5) # number of seconds observed mutliplied by 2
    interpol_length = t_obs * (sampling_rate / 2) # total number of interpolation points

    # Interpolating
    x_num = np.linspace(0,len(realradecf0), num=len(realradecf0), endpoint=True)
    x_interp= np.linspace(0,len(realradecf0), num=int(t_obs)*int((sampling_rate/2)), endpoint=True)

    f = interp1d(x_num, realradecf0, "linear")
    ra_interp = f(x_interp)

    x_num = np.linspace(0,len(realradecf1), num=len(realradecf1), endpoint=True)
    x_interp= np.linspace(0,len(realradecf1), num=int(t_obs)*int((sampling_rate/2)), endpoint=True)

    f = interp1d(x_num, realradecf1, "linear")
    dec_interp = f(x_interp)

    print(len(realradecf1), len(ra_interp))








    ######## AVERAGING THE POSITION
    # realradecf0 = np.array(realradecf0)
    b = len(realradecf0)
    rem = b % 10

    realradecf0 = realradecf0[:-rem]  # working with the thrid column of the intensity table
    realradecf0 = np.mean(np.reshape(realradecf0, (-1, 10)), axis=1)  # averages it out over 32 points

    realradecf1 = realradecf1[:-rem]  # working with the thrid column of the intensity table
    realradecf1 = np.mean(np.reshape(realradecf1, (-1, 10)), axis=1)  # averages it out over 32 points

    ############## NOT AVERAGING BUT TAKING EVERY 10TH VALUE ######################

    # realradecf0 = realradecf0[::10]
    # realradecf1 = realradecf1[::10]

    ###############################################################################




    # realradecf0 = realradecf0[19:-54]  # For MARSx4
    # realradecf1 = realradecf1[19:-54]

    # # realradecf0 = realradecf0[41:-38]  # FOr venus 5/22
    # # realradecf1 = realradecf1[41:-38]

    # xl = len(realradecf0)
    # # print 'the total data points after averaging are', xl, 'hence', xl, ' seconds from position import'
    # # plt.plot(realradecf0, realradecf1)
    # # plt.show()

    # ################## OLD INTERPOLATING 16 POSITIONS BETWEEN EACH SECOND #################
    # realradecf0 = np.array(realradecf0, dtype='float64')
    # realradecf1 = np.array(realradecf1, dtype='float64')



    # realradecf0intpol = []
    # realradecf1intpol = []

    # sampling_rate = 16*4 # 32 for venus

    # for i in range(xl - 1):
    #     if realradecf0[i] < realradecf0[i + 1]:
    #         posf0 = np.linspace(realradecf0[i], realradecf0[i + 1], sampling_rate) #the last element is removed to avoid overlap with the previous point
    #         # print 'POSF0 ---------------------', len(posf0), i
    #         realradecf0intpol.append(
    #             [np.interp(posf0, [realradecf0[i], realradecf0[i + 1]], [realradecf0[i], realradecf0[i + 1]])])
    #         # print 'realradecf0----------------', realradecf0
    #         posf1 = np.linspace(realradecf1[i], realradecf1[i + 1], sampling_rate)
    #         realradecf1intpol.append(
    #             [np.interp(posf1, [realradecf1[i], realradecf1[i + 1]], [realradecf1[i], realradecf1[i + 1]])])
    #     else:
    #         posf0 = np.linspace(realradecf0[i], realradecf0[i + 1], sampling_rate)
    #         # print 'POSF0 ---------------------', len(posf0), i
    #         realradecf0intpol.append(
    #             [np.interp(posf0, [realradecf0[i + 1], realradecf0[i]], [realradecf0[i + 1], realradecf0[i]])])
    #         # print 'realradecf0----------------', realradecf0
    #         posf1 = np.linspace(realradecf1[i + 1], realradecf1[i], sampling_rate)
    #         realradecf1intpol.append(
    #             [np.interp(posf1, [realradecf1[i + 1], realradecf1[i]], [realradecf1[i + 1], realradecf1[i]])])

    ###################################################################################

    # print realradecf1intpol[0], 'realradecf1intpol'

    realradecf0intpol = np.array(ra_interp)
    realradecf1intpol = np.array(dec_interp)

    realradecf0intpol = realradecf0intpol.flatten()
    realradecf1intpol = realradecf1intpol.flatten()

    xl = len(realradecf0intpol)
    print(('Position data points imported:', xl))

    # plt.plot(realradecf0intpol, realradecf1intpol)
    # plt.show()

    # plt.scatter(realradecf0intpol, realradecf1intpol)
    # plt.show()

    """Here we crop the data, first we manually crop, then we use the crop method
        """
    realradecf0intpol = realradecf0intpol[manual_selection[0]:manual_selection[1]]
    realradecf1intpol = realradecf1intpol[manual_selection[0]:manual_selection[1]]

    overall_maxra = np.amax(realradecf0intpol)
    overall_minra = np.amin(realradecf0intpol)
    overall_maxdec = np.amax(realradecf1intpol)
    overall_mindec = np.amin(realradecf1intpol)
    mapcenter = (((overall_minra + (overall_maxra - overall_minra) / 2)),
                 (overall_mindec + (overall_maxdec - overall_mindec) / 2))
    mapsize = ((overall_maxdec - overall_mindec), (overall_maxdec - overall_mindec))

    j=0

    realradecf0intpol_masked = []
    realradecf1intpol_masked = []
    index_cropped = []

    """NEW Cropping method"""

    if crop == 0:
        xboundary = [mapcenter[0] - (mapsize[0]/2), mapcenter[0] + (mapsize[0]/2)]
        yboundary = [mapcenter[1] - (mapsize[1]/2), mapcenter[1] + (mapsize[1]/2)]
    elif len(crop) ==2:
        coords = [realradecf0intpol, realradecf1intpol]
        mapwidth = crop[0] / 60.
        mapheight = crop[1]/ 60.
        xboundary = [mapcenter[0] - (mapwidth / 2), mapcenter[0] + (mapwidth / 2)]
        yboundary = [mapcenter[1] - (mapheight / 2), mapcenter[1] + (mapheight / 2)]
        xmask = ma.masked_outside(coords[0], xboundary[0], xboundary[1])
        ymask = ma.masked_outside(coords[1], yboundary[0], yboundary[1])

        realradecf0intpol_masked = xmask
        realradecf1intpol_masked = ymask

        index_cropped = [j for j,i in enumerate(xmask*ymask)]# & coords[1][j] == realradecf1intpol_masked[j]]


        # print(len(xmask*ymask),'MAMAMAMAMAMA')

        print(len(realradecf0intpol_masked))
        # print(len(index_cropped),'YMASK')

        index =np.arange(0,len(coords[0]),1)


    # # with intensity as intens:
    
    # # uncropped_intensity = intensity_import_old(intensity)
    # # cropped_intensity = intensity_import(intensity, index_cropped=index_cropped) # FIXME: get the intensity if possible
    
    # plt.scatter(realradecf0intpol, realradecf1intpol, c = 'grey', alpha=0.4)
    # plt.scatter(realradecf0intpol_masked, realradecf1intpol_masked, c='green', alpha = 0.7)
    # plt.vlines(x=xboundary[0], ymin=yboundary[0], ymax=yboundary[1])
    # plt.vlines(x=xboundary[1], ymin=yboundary[0], ymax=yboundary[1])
    # plt.hlines(y=yboundary[0], xmin=xboundary[0], xmax=xboundary[1])
    # plt.hlines(y=yboundary[1], xmin=xboundary[0], xmax=xboundary[1])
    # plt.Circle(mapcenter, 0.5 / 60 , color='r')
    # # plt.savefig("cropping_scan_area.pdf")
    # plt.title('manual selection {}'.format(manual_selection))
    # plt.show()



    xl = len(realradecf0intpol_masked)
    print(('Position data points imported:', xl))
    # return realradecf0, realradecf1

    # index_cropped = [i for i, j in enumerate(realradecf1intpol)]

    # print(np.shape(realradecf0intpol_masked))
#     return realradecf0intpol_masked.compressed(), realradecf1intpol_masked.compressed(), index_cropped

    # return np.array(realradecf0intpol_masked), np.array(realradecf1intpol_masked) , index_cropped    # For no-offset beammaps
    # return np.array(realradecf0intpol_masked) + (offset_x/3600.), np.array(realradecf1intpol_masked) - (offset_y/3600.), index_cropped    # Worked for Beam maps
    # return np.array(realradecf0intpol_masked) - (offset_x), np.array(realradecf1intpol_masked) - (offset_y), index_cropped  # removed the degree conversion for RA-Dec


    if maptype == "beammap":
        return np.array(realradecf0intpol_masked) + (offset_x/3600.), np.array(realradecf1intpol_masked) - (offset_y/3600.), index_cropped    # Worked for Beam maps
    elif maptype == "skymap":
        return np.array(realradecf0intpol_masked) - (offset_x), np.array(realradecf1intpol_masked) - (offset_y), index_cropped  # removed the degree conversion for RA-Dec
    else:
        print("Please check *map_type* parameter in the settings")




def position_import(antlog, maptype, manual_selection, crop, sampling_rate):

    """
    Calculates several coordinate related values from the antenna log. Converted to python from the code written in
    C sent by Nagai-san.

    :param antlog: antenna log file from the telescope
    sampling_rate: the number of samples in the intensity per second
    :return: RA-Dec, dAz-dEl, etc.

    Improvements to be made: Commandline selection of the output, maybe by an "if".
    Rotation due to change in elevation calculation.
    """

    LAT_NRO = 0.62729799791922308

    with open(antlog) as f:
        fdata = [line.rstrip() for line in f]  # this splits each row
        #	if "AZEL" in line:
        cdata = [line.split() for line in fdata]  # this splits each column
        a = 1  # this is the number of row
        #	b=2 # this is the number of column in the above row
        # c=len(fdata)
        c = len(fdata)
        # print len(fdata)
        delout = []
        dazout = []
        elout = []
        time = []

        # Observing in AZ-EL------------------
        realazelf0 = []
        realazelf1 = []
        progazelf0 = []
        progazelf1 = []

        # observing in RA-DEC-----------------
        realradecf0 = []
        realradecf1 = []
        realradecval0 = []
        realradecval1 = []

        #########   the start of loop for offsets
        # ofc = len(lens_offsetdeg_x)
        ofb = 0
        dazimuth = []
        dazcol = []

        while (a != (c - 1)):
            a += 1

            # print a
            #		g= fdata[a] #in case of extracting a single row
            n = cdata[a]  # this is splitting the a`th row into the 'n' array

            #	eval(str) converts a string to a number which can be calculated with
            adprog00 = (eval(n[1]) / 180) * np.pi  # RA value
            adprog01 = (eval(n[2]) / 180) * np.pi  # DEC value
            aeprog0 = (eval(n[3]) / 180) * np.pi  # AZ value
            aeprog1 = (eval(n[4]) / 180) * np.pi  # El value
            aeprog00 = (eval(n[9]) / 180) * np.pi  # AZ center indication value(not including offset)
            aeprog01 = (eval(n[10]) / 180) * np.pi  # El center indication value(not including offset)
            aereal0 = (eval(n[5]) / 180) * np.pi  # AZ encoder value(collimator value)
            aereal1 = (eval(n[6]) / 180) * np.pi  # El encoder value(collimator value)

            #	Beam offsets

            dofsx = ((0 / 3600.0) / 180) * np.pi  # x beam offset is -40"
            dofsy = ((0 / 3600.0) / 180) * np.pi  # y beam offset is 40"

            #	Rotaion by elevation

            rxval = -aereal1

            #	Scan offset + beam offset

            dEl = dofsx * np.sin(rxval) + dofsy * np.cos(rxval) + (aereal1 - aeprog01)
            dAz = dofsx * np.cos(rxval) - dofsy * np.sin(rxval) + ((aereal0 - aeprog00)) * np.cos(aereal1)

            dEl = dEl - ((580 / 3600.0) / 180.) * np.pi
            #
            dAz = dAz + ((40 / 3600.0) / 180.) * np.pi
            #	Add offsets

            aereal1 = aeprog01 + dEl
            aereal0 = aeprog00 + dAz / np.cos(aereal1)

            #	Calculate parallactic angle to convert dAZEL to dRADEC
            para = - np.arcsin(np.sin(aereal0) * np.cos(LAT_NRO) / np.cos(adprog01))
            dRa = -np.cos(para) * dAz + np.sin(para) * dEl
            dDec = np.sin(para) * dAz + np.cos(para) * dEl

            #	Add dRa DDec

            adreal1 = adprog01 + dDec
            adreal0 = adprog00 + dRa / np.cos(adreal1)

            #	converting back to proper units

            progazel0 = (aeprog0 / np.pi) * 180
            progazel1 = (aeprog1 / np.pi) * 180
            realazel0 = (aereal0 / np.pi) * 180
            realazel1 = (aereal1 / np.pi) * 180
            dAz = (dAz / np.pi) * 180  # * 3600
            dEl = (dEl / np.pi) * 180  # * 3600
            progradec0 = (adprog00 / np.pi) * 180
            progradec1 = (adprog01 / np.pi) * 180
            realradec0 = (adreal0 / np.pi) * 180
            realradec1 = (adreal1 / np.pi) * 180
            #		appending the calculated values into arrays
            dazout.append(dAz)                 ########## Changing this for now
            delout.append(dEl)
            # elout.append(realazel1)
            realazelf0.append(realazel0)
            realazelf1.append(realazel1)
            progazelf0.append(progazel0)
            progazelf1.append(progazel1)
            realradecf0.append(realradec0)
            realradecf1.append(realradec1)
            time1 = eval(n[0])  # extract time from log
            time1 = repr(time1)  # converts time1 from float to string which is essential for strptime()
            dt = datetime.datetime.strptime(time1, '%y%m%d%H%M%S.%f')
            time.append(dt.time())
            # time.append(time1)
            ofb += 1

        # print ("Antenna log starting time is {}".format(time[0]))
        # print ("Antenna log ending time is {}".format(time[len(time) - 1]))



        # for i in range(len(time)):
        #     print time[i].time()

        # print time,realradecf0, realradecf1

        #




    if maptype == "beammap":
        realradecf0 = dazout                                      ########## CHANGED THIS  ! aND IT WORKED
        realradecf1 = delout
    elif maptype == "skymap":
        pass
    else:
        print("Please check *map_type* parameter in the settings")




    ############*OLD CODE*################# Define the starttime and average the position, and take time every second
    # start_time = datetime.datetime.strptime("05:49:00", "%H:%M:%S").time()
    #
    # # print ("Antenna log starting time is {}".format(time[0]))
    # # print ("Antenna log ending time is {}".format(time[len(time) - 1]))
    # begin_index = [i for i, j in enumerate(time) if j == start_time]
    # # print [j for i, j in enumerate(time) if j == start_time]
    # # print begin_index[0]

    # realradecf1 = realradecf1[begin_index[0]:]
    # realradecf0 = realradecf0[begin_index[0]:]
    # time = time[begin_index[0]:]

    avtime = time[::10]
    # print 'Antlog recorded ', len(avtime), 's'



    realradecf0 = realradecf0[190:-540]  # For MARSx4
    realradecf1 = realradecf1[190:-540]

    # realradecf0 = realradecf0[41:-38]  # FOr venus 5/22
    # realradecf1 = realradecf1[41:-38]

    # ## INTERPOLATING RA DEC SEPARATELY #######


    sampling_rate = sampling_rate
    b = len(realradecf0)
    rem = b % 5   # 5 samples for half a second. ( 10 per second)



    #trimming the end of the position log such that the interpolation can be done properly
    if rem !=0:
        realradecf0 = realradecf0[:-rem] # samples more than half a second are discarded
        realradecf1 = realradecf1[:-rem]


    # The length of ineterpolation
    t_obs = (len(realradecf0) / 5) # number of seconds observed mutliplied by 2
    interpol_length = t_obs * (sampling_rate / 2) # total number of interpolation points

    # Interpolating
    x_num = np.linspace(0,len(realradecf0), num=len(realradecf0), endpoint=True)
    x_interp= np.linspace(0,len(realradecf0), num=int(t_obs)*int((sampling_rate/2)), endpoint=True)

    f = interp1d(x_num, realradecf0, "linear")
    ra_interp = f(x_interp)

    x_num = np.linspace(0,len(realradecf1), num=len(realradecf1), endpoint=True)
    x_interp= np.linspace(0,len(realradecf1), num=int(t_obs)*int((sampling_rate/2)), endpoint=True)

    f = interp1d(x_num, realradecf1, "linear")
    dec_interp = f(x_interp)

    print(len(realradecf1), len(ra_interp))

    ################# PLOTTING THE INTERPOLATED VALUES FOR THESIS ##################

    # plt.plot(realradecf0, realradecf1, 'o-', label='Anetnna-log coordinates')
    # plt.plot(ra_interp, dec_interp, '.-', label='interpolated coordinates')
    # plt.title('Position Interpolation')
    # plt.xlabel('dAz(deg)')
    # plt.ylabel('dEl(deg)')
    # plt.legend()
    # plt.show()








    ######## AVERAGING THE POSITION
    # realradecf0 = np.array(realradecf0)
    b = len(realradecf0)
    rem = b % 10

    realradecf0 = realradecf0[:-rem]  # working with the thrid column of the intensity table
    realradecf0 = np.mean(np.reshape(realradecf0, (-1, 10)), axis=1)  # averages it out over 32 points

    realradecf1 = realradecf1[:-rem]  # working with the thrid column of the intensity table
    realradecf1 = np.mean(np.reshape(realradecf1, (-1, 10)), axis=1)  # averages it out over 32 points

    ############## NOT AVERAGING BUT TAKING EVERY 10TH VALUE ######################

    # realradecf0 = realradecf0[::10]
    # realradecf1 = realradecf1[::10]

    ###############################################################################




    # realradecf0 = realradecf0[19:-54]  # For MARSx4
    # realradecf1 = realradecf1[19:-54]

    # # realradecf0 = realradecf0[41:-38]  # FOr venus 5/22
    # # realradecf1 = realradecf1[41:-38]

    # xl = len(realradecf0)
    # # print 'the total data points after averaging are', xl, 'hence', xl, ' seconds from position import'
    # # plt.plot(realradecf0, realradecf1)
    # # plt.show()



    # print realradecf1intpol[0], 'realradecf1intpol'

    realradecf0intpol = np.array(ra_interp)
    realradecf1intpol = np.array(dec_interp)

    realradecf0intpol = realradecf0intpol.flatten()
    realradecf1intpol = realradecf1intpol.flatten()

    xl = len(realradecf0intpol)
    print(('Position data points imported:', xl))

    # plt.plot(realradecf0intpol, realradecf1intpol)
    # plt.show()

    # plt.scatter(realradecf0intpol, realradecf1intpol)
    # plt.show()

    """Here we crop the data, first we manually crop, then we use the crop method
        """
    realradecf0intpol = realradecf0intpol[manual_selection[0]:manual_selection[1]]
    realradecf1intpol = realradecf1intpol[manual_selection[0]:manual_selection[1]]

    overall_maxra = np.amax(realradecf0intpol)
    overall_minra = np.amin(realradecf0intpol)
    overall_maxdec = np.amax(realradecf1intpol)
    overall_mindec = np.amin(realradecf1intpol)
    mapcenter = (((overall_minra + (overall_maxra - overall_minra) / 2)),
                 (overall_mindec + (overall_maxdec - overall_mindec) / 2))
    mapsize = ((overall_maxdec - overall_mindec), (overall_maxdec - overall_mindec))

    j=0

    realradecf0intpol_masked = []
    realradecf1intpol_masked = []
    index_cropped = []

    """NEW Cropping method"""

    if crop == 0:
        xboundary = [mapcenter[0] - (mapsize[0]/2), mapcenter[0] + (mapsize[0]/2)]
        yboundary = [mapcenter[1] - (mapsize[1]/2), mapcenter[1] + (mapsize[1]/2)]
    elif len(crop) ==2:
        coords = [realradecf0intpol, realradecf1intpol]
        mapwidth = crop[0] / 60.
        mapheight = crop[1]/ 60.
        xboundary = [mapcenter[0] - (mapwidth / 2), mapcenter[0] + (mapwidth / 2)]
        yboundary = [mapcenter[1] - (mapheight / 2), mapcenter[1] + (mapheight / 2)]
        xmask = ma.masked_outside(coords[0], xboundary[0], xboundary[1])
        ymask = ma.masked_outside(coords[1], yboundary[0], yboundary[1])

        realradecf0intpol_masked = xmask
        realradecf1intpol_masked = ymask

        index_cropped = [j for j,i in enumerate(xmask*ymask)]# & coords[1][j] == realradecf1intpol_masked[j]]


        # print(len(xmask*ymask),'MAMAMAMAMAMA')
        #
        # print(len(realradecf0intpol_masked))
        # print(len(index_cropped),'YMASK')

        index =np.arange(0,len(coords[0]),1)


    # # with intensity as intens:
    
    # # uncropped_intensity = intensity_import_old(intensity)
    # # cropped_intensity = intensity_import(intensity, index_cropped=index_cropped) # FIXME: get the intensity if possible

    ################# PLOTTING THE Mask FOR THESIS ##################


    

    # plt.scatter(realradecf0intpol, realradecf1intpol, c = 'grey', alpha=0.4)
    # plt.scatter(realradecf0intpol_masked*3600, realradecf1intpol_masked*3600, c='green', alpha = 0.7, label='Scanning Region')
    # plt.vlines(x=xboundary[0], ymin=yboundary[0], ymax=yboundary[1])
    # plt.vlines(x=xboundary[1], ymin=yboundary[0], ymax=yboundary[1])
    # plt.hlines(y=yboundary[0], xmin=xboundary[0], xmax=xboundary[1])
    # plt.hlines(y=yboundary[1], xmin=xboundary[0], xmax=xboundary[1])
    # plt.Circle(mapcenter, 0.5 , color='r', label= 'MapCenter')
    # # plt.savefig("cropping_scan_area.pdf")
    # plt.title('manual selection {}'.format(manual_selection))
    # plt.xlabel('dEl (deg)')
    # plt.ylabel('dAz (deg)')
    # plt.legend()
    # plt.show()



    xl = len(realradecf0intpol_masked)
    print(('Position data points imported:', xl))
    # return realradecf0, realradecf1

    # index_cropped = [i for i, j in enumerate(realradecf1intpol)]
    return realradecf0intpol_masked, realradecf1intpol_masked, index_cropped

def intensity_import(intensity_log, index_cropped, offset = 0, invert = True, normalize = True):

    rdata = intensity_log
    ######### OLD chunk to read intensity from a txt file ####################
    # rdata = np.loadtxt(intensity_log)
    # # row = [line.rstrip() for line in intensity_log]
    # # rdata = [line.split() for line in row]
    # # rdata = np.array(rdata).astype(np.float)
    # rdata = rdata[:, 2]
    ##########################################################################

    print(rdata)

    b = len(rdata)
    # c = len(row)
    ########## CHANGED FOR VENUS, NEED TO UPDATE FOR EVERYDATA POINTS

    rem = b % 64 #making sure there are 16 intensity datapoints per second

    # print rem, 'rem'

    if rem != 0:
        y = rdata[:-rem]  # working with the thrid column of the intensity table
    else:
        y = rdata

    offset = offset
    top_crop = (12*64)+ offset
    bottom_crop= - (10*64) + offset
    y = y[top_crop:bottom_crop]  # 52, 35   [best position average 170:-6] [best position for 10th value 163:-13]
    # # print(len(a), 'LENGTH OF A')


    '''Now for the sake of taking the dB, we subtract the value from f0 and find df'''
    ymin = np.amin(y)

    if invert==True:
        y = -1 * y # intensity inversion
    else:
        pass

    if normalize==True:
        ymin = np.amin(y)
        ymax = np.amax(y)
        y = y/ymax
#         y = (y - ymin) / (ymax - ymin)
    else: pass

    # Cropping

    y = [y[i] for i in index_cropped]


    # plt.plot(y)
    # plt.title('ymax{}, ymin{}'.format(ymax, ymin))
    # plt.show()

    xl = len(y)
    print(('Intensity data points imported:', xl))

    # return 100* np.log10(y)
    return y

def intensity_import_old(intensity_log):

    """
    Reads frequency offset values of MKIDs. So far written for PCA data in the form of
    single file per MKID.

    :param intensity_log: filename of the frequency offset values.
    :return: cropped, normalized, dB value of the frequency offset

    Improvements to be made: include "if" to give log or non log values, or maybe just log them later.
    The cropping method needs to be improved for combined map image synthesis.
    """

    # Import Data

    row = [line.rstrip() for line in intensity_log]
    rdata = [line.split() for line in row]
    rdata = np.array(rdata).astype(np.float)
    rdata = rdata[2]


    b = len(rdata)
    c = len(row)

    rem = b % 32

    # print rdata

    if rem != 0:
        y = rdata[:-rem]  # working with the thrid column of the intensity table
    else:
        y = rdata
    print((y, 'y'))
    #################################################################################
    # y = np.mean(np.reshape(y, (-1, 16)), axis=1)  # averages it out over 16 points
    #################################################################################
    # xl = len(y)
    # print 'the total data points after averaging are', xl, 'hence', xl, ' seconds from intensity import'


    # plt.plot(y)
    # plt.xlabel('Time (sec)')
    # plt.ylabel('Freq (GHz)')
    # plt.title('Intensity(averaged) Plot')
    # plt.show()

    y = y[172:-4]  # 52, 35   [best position average 170:-6] [best position for 10th value 163:-13]
    # print len(y)
    xl = len(y)
    print(('Intensity data points imported:', xl))

    '''Now for the sake of taking the dB, we subtract the value from f0 and find df'''
    ymin = np.amin(y)

    # print ymin, 'ymin'


    # y = 3869.01726091842 - y
    # print y
    # cyb0 = 0
    # for i in range(len(y)-1):
    #     if y[i] < 0.0:
    #         y[i] = 0.00001

    # print cyb0, 'less than 0 values'
    y = -1 * y

    ymin = np.amin(y)

    # plt.plot(y)
    #
    # # y = y + abs(2*ymin)
    # plt.plot(y)
    # plt.show()
    ymin = np.amin(y)
    # print ymin, 'ymin'

    # print np.count_nonzero(y), 'non zero y'
    ymin = np.amin(y)

    # print ymin, 'ymin'

    ymax = np.max(y)
    y = (y - ymin) / (ymax - ymin)

    y = y + abs(ymin)

    ymin = np.amin(y)

    print(y)

    # plt.plot(y)
    # plt.show()

    # print ymin, 'ymin'


    # print np.count_nonzero(y), 'Non zero y'
    # y = y - (ymin + 0.001)
    # print np.bincount(y, 0.0), 'Non zero y'
    y = 10 * np.log10(y)

    # plt.plot(y)
    # plt.show()

    ymin = np.amin(y)

    ##################### Normalized #################
    # ymin = np.amin(y)
    #
    # print ymin, 'ymin'
    # ymax = np.max(y)
    # y = (y - ymin) / (ymax - ymin)
    # print ymin, 'ymin after dB'
    for i in range(len(y) - 1):
        if y[i] == ymin:  # ymin is -inf
            y[i] = -14
    # plt.plot(y)
    # plt.show()

    # print y
    return y



    # Since every second has 32 data points, it has been averaged into the intensity for each second

def setup_header(mapcenter, mapsize, beamsize_fwhm):
    '''
    Produce a FITS header that contains the target field.

    :param mapcenter: calculated from the position_log
    :param mapsize: calculated from the position_log.
    :param beamsize_fwhm: defined in the gridding section of mapping.py.
    :return: FITS header for the target field.
    '''

    # define target grid (via fits according to WCS convention)
    # a good size is a third of the FWHM of the PSF (avoids aliasing)

    pixsize = beamsize_fwhm / 5.
    dnaxis1 = int(mapsize[0] / pixsize)
    dnaxis2 = int(mapsize[1] / pixsize)

    header = {
        'NAXIS': 3,
        'NAXIS1': dnaxis1,
        'NAXIS2': dnaxis2,
        'NAXIS3': 1,  # need a dummy spectral axes
        'CTYPE1': 'RA---SIN',
        'CTYPE2': 'DEC--SIN',
        'CUNIT1': 'deg',
        'CUNIT2': 'deg',
        'CDELT1': -pixsize,
        'CDELT2': pixsize,
        'CRPIX1': (dnaxis1 + 1) / 2.,
        'CRPIX2': (dnaxis2 + 1) / 2.,
        'CRVAL1': mapcenter[0],
        'CRVAl2': mapcenter[1],
    }

    # self._header(header=header)

    return header
def position_import_offset(antlogfile, manual_selection, crop, offset_x, offset_y): #Exclusively used for CombinedMKIDMap class, adding offset values to position coordinates
    """
    Same as the position_import function with added offset values. The offset values are read from a file. Calculated by
    Beammapcenter class of mapping.py and beammap_position_size_plot.py.

    :param antlogfile:  antenna log file from the telescope.
    :param x: the index number of the MKID.
    :return: offset added RA-Dec, dAz-dEl, etc.
    """

    """Now, let us move on to calculate the RA-Dec from the antenna log"""

    LAT_NRO = 0.62729799791922308

    # offsets = x

    with open(antlogfile) as antlog:
        fdata = [line.rstrip() for line in antlog]  # this splits each row
        #	if "AZEL" in line:
        cdata = [line.split() for line in fdata]  # this splits each column
        a = 1  # this is the number of row
        #	b=2 # this is the number of column in the above row
        # c=len(fdata)
        c = len(fdata)
        # print len(fdata)
        delout = []
        dazout = []
        elout = []
        time = []

        # Observing in AZ-EL------------------
        realazelf0 = []
        realazelf1 = []
        progazelf0 = []
        progazelf1 = []

        # observing in RA-DEC-----------------
        realradecf0 = []
        realradecf1 = []
        realradecval0 = []
        realradecval1 = []

        #########   the start of loop for offsets
        # ofc = len(lens_offsetdeg_x)
        ofb = 0
        dazimuth = []
        dazcol = []

        while (a != (c - 1)):
            a += 1

            # print a
            #		g= fdata[a] #in case of extracting a single row
            n = cdata[a]  # this is splitting the a`th row into the 'n' array

            #	eval(str) converts a string to a number which can be calculated with
            adprog00 = (eval(n[1]) / 180) * np.pi  # RA value
            adprog01 = (eval(n[2]) / 180) * np.pi  # DEC value
            aeprog0 = (eval(n[3]) / 180) * np.pi  # AZ value
            aeprog1 = (eval(n[4]) / 180) * np.pi  # El value
            aeprog00 = (eval(n[9]) / 180) * np.pi  # AZ center indication value(not including offset)
            aeprog01 = (eval(n[10]) / 180) * np.pi  # El center indication value(not including offset)
            aereal0 = (eval(n[5]) / 180) * np.pi  # AZ encoder value(collimator value)
            aereal1 = (eval(n[6]) / 180) * np.pi  # El encoder value(collimator value)

            #	Beam offsets

            #	Beam offsets

            dofsx = ((0 / 3600.0) / 180) * np.pi  # x beam offset is -40"
            dofsy = ((0 / 3600.0) / 180) * np.pi  # y beam offset is 40"

            #	Rotaion by elevation

            rxval = -aereal1

            #	Scan offset + beam offset

            dEl = dofsx * np.sin(rxval) + dofsy * np.cos(rxval) + (aereal1 - aeprog01)
            dAz = dofsx * np.cos(rxval) - dofsy * np.sin(rxval) + ((aereal0 - aeprog00)) * np.cos(aereal1)

            dEl = dEl - ((580 / 3600.0) / 180.) * np.pi
            #
            dAz = dAz + ((40 / 3600.0) / 180.) * np.pi

            ##########################OLD Code############################

            # # dofsx = (offsets[x, 0] / 3600.0)
            # # dofsy = (offsets[x, 1] / 3600.0)
            # #
            # # dofsxdif = ((offset_x / 3600.0) / 180.) * np.pi  # x beam offset is -40"
            # # dofsydif = ((offset_y / 3600.0) / 180.) * np.pi  # y beam offset is 40"
            # #
            # dofsx = ((0 / 3600.0) / 180.) * np.pi  # x beam offset is -40"
            # dofsy = ((0 / 3600.0) / 180.) * np.pi  # y beam offset is 40"
            #
            # #	Rotaion by elevation
            #
            # rxval = -aereal1
            #
            # #	Scan offset + beam offset
            #
            # dEl = dofsx * np.sin(rxval) + dofsy * np.cos(rxval) + (aereal1 - aeprog01)
            # dAz = dofsx * np.cos(rxval) - dofsy * np.sin(rxval) + ((aereal0 - aeprog00)) * np.cos(aereal1)
            #
            # # dEl = dEl - ((580/ 3600.0) / 180.) * np.pi - dofsydif
            # # #
            # # dAz = dAz + ((20/ 3600.0) / 180.) * np.pi + dofsxdif
            ###################################################################

            # print dEl, dAz


            #	Add offsets

            aereal1 = aeprog01 + dEl
            aereal0 = aeprog00 + dAz / np.cos(aereal1)

            #	Calculate parallactic angle to convert dAZEL to dRADEC
            para = - np.arcsin(np.sin(aereal0) * np.cos(LAT_NRO) / np.cos(adprog01))
            dRa = -np.cos(para) * dAz + np.sin(para) * dEl
            dDec = np.sin(para) * dAz + np.cos(para) * dEl

            #	Add dRa DDec

            adreal1 = adprog01 + dDec
            adreal0 = adprog00 + dRa / np.cos(adreal1)

            #	converting back to proper units

            progazel0 = (aeprog0 / np.pi) * 180
            progazel1 = (aeprog1 / np.pi) * 180
            realazel0 = (aereal0 / np.pi) * 180
            realazel1 = (aereal1 / np.pi) * 180
            dAz = (dAz / np.pi) * 180# * 3600
            dEl = (dEl / np.pi) * 180# * 3600
            progradec0 = (adprog00 / np.pi) * 180
            progradec1 = (adprog01 / np.pi) * 180
            realradec0 = (adreal0 / np.pi) * 180
            realradec1 = (adreal1 / np.pi) * 180
            #		appending the calculated values into arrays
            delout.append(dEl)
            dazout.append(dAz)
            # elout.append(realazel1)
            realazelf0.append(realazel0)
            realazelf1.append(realazel1)
            progazelf0.append(progazel0)
            progazelf1.append(progazel1)
            realradecf0.append(realradec0)
            realradecf1.append(realradec1)
            time1 = eval(n[0])  # extract time from log
            time1 = repr(time1)  # converts time1 from float to string which is essential for strptime()
            dt = datetime.datetime.strptime(time1, '%y%m%d%H%M%S.%f')
            time.append(dt.time())
            # time.append(time1)
            ofb += 1

        # print ("Antenna log starting time is {}".format(time[0]))
        # print ("Antenna log ending time is {}".format(time[len(time) - 1]))



        # for i in range(len(time)):
        #     print time[i].time()

        # print time,realradecf0, realradecf1

        #

        ################################# plot the antenna log #####################

    # print realazelf0


    # plt.xlabel('Az')
    # plt.ylabel('El')
    # plt.title('Venus observation Az-El position in the sky')
    # dazout = dazout * 3600
    # delout = delout * 3600
    # plt.scatter(dazout, delout)
    # plt.title('dEl-dAz with beam-offset (-60,60)')
    # plt.xlabel('dAz(deg)')
    # plt.ylabel('dEl(deg)')
    # # plt.scatter(realradecf0, realradecf1)
    # plt.show()

    realradecf0 = dazout
    realradecf1 = delout
        ################################ Define the starttime and average the position, and take time every second
    start_time = datetime.datetime.strptime("05:49:00", "%H:%M:%S").time()

    # print ("Antenna log starting time is {}".format(time[0]))
    # print ("Antenna log ending time is {}".format(time[len(time) - 1]))
    begin_index = [i for i, j in enumerate(time) if j == start_time]
    # print [j for i, j in enumerate(time) if j == start_time]
    # print begin_index[0]

    realradecf1 = realradecf1[begin_index[0]:]
    realradecf0 = realradecf0[begin_index[0]:]
    time = time[begin_index[0]:]

    avtime = time[::10]
    # print 'Antlog recorded ', len(avtime), 's'

    ######## AVERAGING THE POSITION
    # realradecf0 = np.array(realradecf0)
    b = len(realradecf0)
    rem = b % 10

    realradecf0 = realradecf0[:-rem]  # working with the thrid column of the intensity table
    realradecf0 = np.mean(np.reshape(realradecf0, (-1, 10)), axis=1)  # averages it out over 32 points

    realradecf1 = realradecf1[:-rem]  # working with the thrid column of the intensity table
    realradecf1 = np.mean(np.reshape(realradecf1, (-1, 10)), axis=1)  # averages it out over 32 points


    ############## NOT AVERAGING BUT TAKING EVERY 10TH VALUE ######################

    # realradecf0 = realradecf0[::10]
    # realradecf1 = realradecf1[::10]

    ###############################################################################






    realradecf0 = realradecf0[19:-54] # [29:24]
    realradecf1 = realradecf1[19:-54]
    xl = len(realradecf0)
    # print 'the total data points after averaging are', xl, 'hence', xl, ' seconds from position import'
    # plt.plot(realradecf0, realradecf1)
    # plt.show()

    ################## INTERPOLATING 16 POSITIONS BETWEEN EACH SECOND #################
    realradecf0 = np.array(realradecf0, dtype='float64')
    realradecf1 = np.array(realradecf1, dtype='float64')

    realradecf0intpol = []
    realradecf1intpol = []

    for i in range(xl-1):
        if realradecf0[i] < realradecf0[i + 1]:
            posf0 = np.linspace(realradecf0[i], realradecf0[i + 1], 16)
            # print 'POSF0 ---------------------', len(posf0), i
            realradecf0intpol.append([np.interp(posf0, [realradecf0[i], realradecf0[i + 1]], [realradecf0[i], realradecf0[i + 1]])])
            # print 'realradecf0----------------', realradecf0
            posf1 = np.linspace(realradecf1[i], realradecf1[i + 1], 16)
            realradecf1intpol.append([np.interp(posf1, [realradecf1[i], realradecf1[i + 1]], [realradecf1[i], realradecf1[i + 1]])])
        else:
            posf0 = np.linspace(realradecf0[i], realradecf0[i + 1], 16)
            # print 'POSF0 ---------------------', len(posf0), i
            realradecf0intpol.append(
                [np.interp(posf0, [realradecf0[i + 1], realradecf0[i]], [realradecf0[i + 1], realradecf0[i]])])
            # print 'realradecf0----------------', realradecf0
            posf1 = np.linspace(realradecf1[i + 1], realradecf1[i], 16)
            realradecf1intpol.append(
                [np.interp(posf1, [realradecf1[i + 1], realradecf1[i]], [realradecf1[i + 1], realradecf1[i]])])

    ###################################################################################

    # print realradecf1intpol[0], 'realradecf1intpol'

    realradecf0intpol = np.array(realradecf0intpol)
    realradecf1intpol = np.array(realradecf1intpol)



    realradecf0intpol = realradecf0intpol.flatten()
    realradecf1intpol = realradecf1intpol.flatten()


    xl = len(realradecf0intpol)
    print(('Position points imported:', xl))

    """Here we crop the data, first we manually crop, then we use the crop method
            """
    realradecf0intpol = realradecf0intpol[manual_selection[0]:manual_selection[1]]
    realradecf1intpol = realradecf1intpol[manual_selection[0]:manual_selection[1]]

    overall_maxra = np.amax(realradecf0intpol)
    overall_minra = np.amin(realradecf0intpol)
    overall_maxdec = np.amax(realradecf1intpol)
    overall_mindec = np.amin(realradecf1intpol)
    mapcenter = (((overall_minra + (overall_maxra - overall_minra) / 2)),
                 (overall_mindec + (overall_maxdec - overall_mindec) / 2))
    mapsize = ((overall_maxdec - overall_mindec), (overall_maxdec - overall_mindec))

    j = 0

    realradecf0intpol_masked = []
    realradecf1intpol_masked = []
    index_cropped = []

    """NEW Cropping method"""

    if crop == 0:
        xboundary = [mapcenter[0] - (mapsize[0] / 2), mapcenter[0] + (mapsize[0] / 2)]
        yboundary = [mapcenter[1] - (mapsize[1] / 2), mapcenter[1] + (mapsize[1] / 2)]
    elif len(crop) == 2:
        coords = [realradecf0intpol, realradecf1intpol]
        mapwidth = crop[0] / 60.
        mapheight = crop[1] / 60.
        xboundary = [mapcenter[0] - (mapwidth / 2), mapcenter[0] + (mapwidth / 2)]
        yboundary = [mapcenter[1] - (mapheight / 2), mapcenter[1] + (mapheight / 2)]
        xmask = ma.masked_outside(coords[0], xboundary[0], xboundary[1])
        ymask = ma.masked_outside(coords[1], yboundary[0], yboundary[1])

        realradecf0intpol_masked = xmask
        realradecf1intpol_masked = ymask

        index_cropped = [j for j, i in enumerate(xmask * ymask)]  # & coords[1][j] == realradecf1intpol_masked[j]]

        print(len(xmask * ymask), 'MAMAMAMAMAMA')

        print(len(realradecf0intpol_masked))
        print(len(index_cropped), 'YMASK')

        index = np.arange(0, len(coords[0]), 1)

    # with intensity as intens:
    
    # uncropped_intensity = intensity_import_old(intensity)
    # cropped_intensity = intensity_import(intensity, index_cropped=index_cropped) # FIXME: get the intensity if possible
    
    plt.scatter(realradecf0intpol, realradecf1intpol, c = 'grey', alpha=0.4)
    plt.scatter(realradecf0intpol_masked, realradecf1intpol_masked, c='green', alpha = 0.7)
    plt.vlines(x=xboundary[0], ymin=yboundary[0], ymax=yboundary[1])
    plt.vlines(x=xboundary[1], ymin=yboundary[0], ymax=yboundary[1])
    plt.hlines(y=yboundary[0], xmin=xboundary[0], xmax=xboundary[1])
    plt.hlines(y=yboundary[1], xmin=xboundary[0], xmax=xboundary[1])
    plt.Circle(mapcenter, 0.5 / 60 , color='r')
    # plt.savefig("cropping_scan_area.pdf")
    plt.title('manual selection {}'.format(manual_selection))
    plt.show()

    xl = len(realradecf0intpol_masked)
    print(('Position data points imported:', xl))
    # return realradecf0, realradecf1

#     index_cropped = [i for i, j in enumerate(realradecf1intpol)]
    return realradecf0intpol_masked, realradecf1intpol_masked, index_cropped

#     return realradecf0intpol_masked + (offset_x/3600.), realradecf1intpol_masked + (offset_y/3600.), index_cropped

def position_import_offset_old(antlogfile, x): #Exclusively used for CombinedMKIDMap class, adding offset values to position coordinates
    """
    Same as the position_import function with added offset values. The offset values are read from a file. Calculated by
    Beammapcenter class of mapping.py and beammap_position_size_plot.py.

    :param antlogfile:  antenna log file from the telescope.
    :param x: the index number of the MKID.
    :return: offset added RA-Dec, dAz-dEl, etc.
    """

    """Now, let us move on to calculate the RA-Dec from the antenna log"""

    LAT_NRO = 0.62729799791922308

    offsets = np.loadtxt("/home/pranshu/Downloads/mappingpackage/data/42nd_pixel_set_abs_0.txt")

    with open(antlogfile) as antlog:
        fdata = [line.rstrip() for line in antlog]  # this splits each row
        #	if "AZEL" in line:
        cdata = [line.split() for line in fdata]  # this splits each column
        a = 1  # this is the number of row
        #	b=2 # this is the number of column in the above row
        # c=len(fdata)
        c = len(fdata)
        # print len(fdata)
        delout = []
        dazout = []
        elout = []
        time = []

        # Observing in AZ-EL------------------
        realazelf0 = []
        realazelf1 = []
        progazelf0 = []
        progazelf1 = []

        # observing in RA-DEC-----------------
        realradecf0 = []
        realradecf1 = []
        realradecval0 = []
        realradecval1 = []

        #########   the start of loop for offsets
        # ofc = len(lens_offsetdeg_x)
        ofb = 0
        dazimuth = []
        dazcol = []

        while (a != (c - 1)):
            a += 1

            # print a
            #		g= fdata[a] #in case of extracting a single row
            n = cdata[a]  # this is splitting the a`th row into the 'n' array

            #	eval(str) converts a string to a number which can be calculated with
            adprog00 = (eval(n[1]) / 180) * np.pi  # RA value
            adprog01 = (eval(n[2]) / 180) * np.pi  # DEC value
            aeprog0 = (eval(n[3]) / 180) * np.pi  # AZ value
            aeprog1 = (eval(n[4]) / 180) * np.pi  # El value
            aeprog00 = (eval(n[9]) / 180) * np.pi  # AZ center indication value(not including offset)
            aeprog01 = (eval(n[10]) / 180) * np.pi  # El center indication value(not including offset)
            aereal0 = (eval(n[5]) / 180) * np.pi  # AZ encoder value(collimator value)
            aereal1 = (eval(n[6]) / 180) * np.pi  # El encoder value(collimator value)

            #	Beam offsets

            # dofsx = (offsets[x, 0] / 3600.0)
            # dofsy = (offsets[x, 1] / 3600.0)
            #
            dofsxdif = ((offsets[x, 0] / 3600.0) / 180.) * np.pi  # x beam offset is -40"
            dofsydif = ((offsets[x, 1] / 3600.0) / 180.) * np.pi  # y beam offset is 40"
            #
            dofsx = ((0 / 3600.0) / 180.) * np.pi  # x beam offset is -40"
            dofsy = ((0 / 3600.0) / 180.) * np.pi  # y beam offset is 40"

            #	Rotaion by elevation

            rxval = -aereal1

            #	Scan offset + beam offset

            dEl = dofsx * np.sin(rxval) + dofsy * np.cos(rxval) + (aereal1 - aeprog01)
            dAz = dofsx * np.cos(rxval) - dofsy * np.sin(rxval) + ((aereal0 - aeprog00)) * np.cos(aereal1)

            dEl = dEl - ((710/ 3600.0) / 180.) * np.pi - dofsydif
            #
            dAz = dAz + ((20/ 3600.0) / 180.) * np.pi + dofsxdif

            # print dEl, dAz


            #	Add offsets

            aereal1 = aeprog01 + dEl
            aereal0 = aeprog00 + dAz / np.cos(aereal1)

            #	Calculate parallactic angle to convert dAZEL to dRADEC
            para = - np.arcsin(np.sin(aereal0) * np.cos(LAT_NRO) / np.cos(adprog01))
            dRa = -np.cos(para) * dAz + np.sin(para) * dEl
            dDec = np.sin(para) * dAz + np.cos(para) * dEl

            #	Add dRa DDec

            adreal1 = adprog01 + dDec
            adreal0 = adprog00 + dRa / np.cos(adreal1)

            #	converting back to proper units

            progazel0 = (aeprog0 / np.pi) * 180
            progazel1 = (aeprog1 / np.pi) * 180
            realazel0 = (aereal0 / np.pi) * 180
            realazel1 = (aereal1 / np.pi) * 180
            dAz = (dAz / np.pi) * 180# * 3600
            dEl = (dEl / np.pi) * 180# * 3600
            progradec0 = (adprog00 / np.pi) * 180
            progradec1 = (adprog01 / np.pi) * 180
            realradec0 = (adreal0 / np.pi) * 180
            realradec1 = (adreal1 / np.pi) * 180
            #		appending the calculated values into arrays
            delout.append(dEl)
            dazout.append(dAz)
            # elout.append(realazel1)
            realazelf0.append(realazel0)
            realazelf1.append(realazel1)
            progazelf0.append(progazel0)
            progazelf1.append(progazel1)
            realradecf0.append(realradec0)
            realradecf1.append(realradec1)
            time1 = eval(n[0])  # extract time from log
            time1 = repr(time1)  # converts time1 from float to string which is essential for strptime()
            dt = datetime.datetime.strptime(time1, '%y%m%d%H%M%S.%f')
            time.append(dt.time())
            # time.append(time1)
            ofb += 1

        # print ("Antenna log starting time is {}".format(time[0]))
        # print ("Antenna log ending time is {}".format(time[len(time) - 1]))



        # for i in range(len(time)):
        #     print time[i].time()

        # print time,realradecf0, realradecf1

        #

        ################################# plot the antenna log #####################

    # print realazelf0


    # plt.xlabel('Az')
    # plt.ylabel('El')
    # plt.title('Venus observation Az-El position in the sky')
    # dazout = dazout * 3600
    # delout = delout * 3600
    # plt.scatter(dazout, delout)
    # plt.title('dEl-dAz with beam-offset (-60,60)')
    # plt.xlabel('dAz(deg)')
    # plt.ylabel('dEl(deg)')
    # # plt.scatter(realradecf0, realradecf1)
    # plt.show()

    realradecf0 = dazout
    realradecf1 = delout
        ################################ Define the starttime and average the position, and take time every second
    start_time = datetime.datetime.strptime("05:49:00", "%H:%M:%S").time()

    # print ("Antenna log starting time is {}".format(time[0]))
    # print ("Antenna log ending time is {}".format(time[len(time) - 1]))
    begin_index = [i for i, j in enumerate(time) if j == start_time]
    # print [j for i, j in enumerate(time) if j == start_time]
    # print begin_index[0]

    realradecf1 = realradecf1[begin_index[0]:]
    realradecf0 = realradecf0[begin_index[0]:]
    time = time[begin_index[0]:]

    avtime = time[::10]
    # print 'Antlog recorded ', len(avtime), 's'

    ######## AVERAGING THE POSITION
    # realradecf0 = np.array(realradecf0)
    b = len(realradecf0)
    rem = b % 10

    realradecf0 = realradecf0[:-rem]  # working with the thrid column of the intensity table
    realradecf0 = np.mean(np.reshape(realradecf0, (-1, 10)), axis=1)  # averages it out over 32 points

    realradecf1 = realradecf1[:-rem]  # working with the thrid column of the intensity table
    realradecf1 = np.mean(np.reshape(realradecf1, (-1, 10)), axis=1)  # averages it out over 32 points


    ############## NOT AVERAGING BUT TAKING EVERY 10TH VALUE ######################

    # realradecf0 = realradecf0[::10]
    # realradecf1 = realradecf1[::10]

    ###############################################################################






    realradecf0 = realradecf0[19:-54] # [29:24]
    realradecf1 = realradecf1[19:-54]
    xl = len(realradecf0)
    # print 'the total data points after averaging are', xl, 'hence', xl, ' seconds from position import'
    # plt.plot(realradecf0, realradecf1)
    # plt.show()

    ################## INTERPOLATING 16 POSITIONS BETWEEN EACH SECOND #################
    realradecf0 = np.array(realradecf0, dtype='float64')
    realradecf1 = np.array(realradecf1, dtype='float64')

    realradecf0intpol = []
    realradecf1intpol = []

    for i in range(xl-1):
        if realradecf0[i] < realradecf0[i + 1]:
            posf0 = np.linspace(realradecf0[i], realradecf0[i + 1], 16)
            # print 'POSF0 ---------------------', len(posf0), i
            realradecf0intpol.append([np.interp(posf0, [realradecf0[i], realradecf0[i + 1]], [realradecf0[i], realradecf0[i + 1]])])
            # print 'realradecf0----------------', realradecf0
            posf1 = np.linspace(realradecf1[i], realradecf1[i + 1], 16)
            realradecf1intpol.append([np.interp(posf1, [realradecf1[i], realradecf1[i + 1]], [realradecf1[i], realradecf1[i + 1]])])
        else:
            posf0 = np.linspace(realradecf0[i], realradecf0[i + 1], 16)
            # print 'POSF0 ---------------------', len(posf0), i
            realradecf0intpol.append(
                [np.interp(posf0, [realradecf0[i + 1], realradecf0[i]], [realradecf0[i + 1], realradecf0[i]])])
            # print 'realradecf0----------------', realradecf0
            posf1 = np.linspace(realradecf1[i + 1], realradecf1[i], 16)
            realradecf1intpol.append(
                [np.interp(posf1, [realradecf1[i + 1], realradecf1[i]], [realradecf1[i + 1], realradecf1[i]])])

    ###################################################################################

    # print realradecf1intpol[0], 'realradecf1intpol'

    realradecf0intpol = np.array(realradecf0intpol)
    realradecf1intpol = np.array(realradecf1intpol)



    realradecf0intpol = realradecf0intpol.flatten()
    realradecf1intpol = realradecf1intpol.flatten()


    xl = len(realradecf0intpol)
    print(('Position points imported:', xl))

    # plt.plot(realradecf0intpol, realradecf1intpol)
    # plt.show()

    # plt.scatter(realradecf0intpol, realradecf1intpol)
    # plt.show()



    #
    # # return realradecf0, realradecf1
    #
    # # realradecf0intpol = np.array(realradecf0intpol)
    # # realradecf1intpol = np.array(realradecf1intpol)
    # # return realradecf0intpol - (offsets[x, 0]), realradecf1intpol - (offsets[x, 1])
    #
    # print ((np.min(realradecf0intpol)* 3600) / 180 ) * np.pi, 'dAz'
    # print ((np.min(realradecf1intpol)* 3600) / 180 ) * np.pi, 'dEl'
    # print offsets[x, 0], offsets[x, 1]

    return realradecf0intpol, realradecf1intpol

def antenna_temp(intensity, t_atm, frequency_table, eta_mb, pixel):

    f_load_av = []
    index = []
    etamb = []

    for i, k in range(enumerate(eta_mb)):
        if i+1 == pixel:
            index = i+1
            etamb = k
    #     print(index)
    j = int(index - 1)
    f_load = np.mean([frequency_table[j][-1], frequency_table[j][-2]])
    ta_star = t_atm * ((intensity) / (f_load))
    t_mb = ta_star/etamb
    return ta_star, t_mb


def temperature_calculator(intensity_log, index_cropped, frequencies, mkidid): #Exclusively used for chopperwheel method temperature calculation

    print(frequencies)
    row = [line.rstrip() for line in intensity_log]
    print((row, 'rdata'))
    rdata = [line.split() for line in row]
    rdata = np.array(rdata).astype(np.float)

    rdata = rdata[:, 2]

    # print rdata

    b = len(rdata)
    c = len(row)
    ########## CHANGED FOR VENUS, NEED TO UPDATE FOR EVERYDATA POINTS

    rem = b % 32  # making sure there are 16 intensity datapoints per second

    # print rem, 'rem'

    if rem != 0:
        y = rdata[:-rem]  # working with the thrid column of the intensity table
    else:
        y = rdata

    offset = 66
    y = y[
        1410 + offset:-1278 + offset]  # 52, 35   [best position average 170:-6] [best position for 10th value 163:-13]
    # print len(y)

    # Cropping

    if index_cropped ==[]:
        pass
    elif len(index_cropped) != 0:
        y = [y[i] for i in index_cropped]

    '''Temperature calculation: Tb = (f_off - f(here 'y'))/(f_load - f_off)'''
    #import frequecies from the frequency file
    if mkidid == frequencies[0]:
        t_atm = 300 # K
        f_brightest = frequencies[1]
        f_off = frequencies[7]
        f_load1 = frequencies[3]
        f_load2 = frequencies[5]

        t_b = t_atm * ((f_off - y) / (((f_load1 + f_load2) / 2.) - f_off))


        ################################### FOR PLOTTING - T_b and frequency shift ###############################
        # fig, ax1 = plt.subplots()
        #
        # color = 'tab:red'
        # ax1.set_xlabel('time (s)')
        # ax1.set_ylabel('- T_b', color=color)
        # ax1.plot(t_b, color=color)
        # ax1.tick_params(axis='y', labelcolor=color)
        #
        # ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        #
        # color = 'tab:blue'
        # ax2.set_ylabel('frequency shift', color=color)  # we already handled the x-label with ax1
        # ax2.plot(y, '.', color=color, alpha=0.4)
        # ax2.tick_params(axis='y', labelcolor=color)
        #
        # fig.tight_layout()  # otherwise the right y-label is slightly clipped
        # plt.title('T_b and frequency shift, MKID {:03}'.format(mkidid))
        # plt.savefig("tbandfreqmkid{:03}.pdf".format(mkidid))
        # # plt.show()
        #
        # # plt.plot(-1. * t_b)
        # # plt.show()
        # plt.close()
        ############################################################################################################

        # the temperature is then converted back to postive
        print(("Intensity data points inported", len(t_b)))
        return (-1. * t_b)
    else:
        print(('ERROR: Problem in frequency file indexing, MKID {} frequencies are not there'.format(mkidid)))

def main_beam_efficiency(Ta_star, T_source, theta_eq, theta_pol, theta_mb):

    eta_mb = Ta_star / ( T_source * (1 - np.exp(-np.log(2)*((theta_eq * theta_pol)/theta_mb))))

def mkids_excluder(intensity,frequencies, mkidid):
    print((mkidid, frequencies))
    if mkidid == frequencies[0]:
        f_brigtest = frequencies[1]
        print(f_brigtest)
        if f_brigtest in intensity:
            print('Beam center in scan area')
        else:
            print('Beam center not in scan area')
            return mkidid
    else:
        print(('ERROR: Problem in frequency file indexing, MKID {} frequencies are not there'.format(mkidid)))


