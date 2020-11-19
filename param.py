#Params
Param = {
'LAT_NRO' : 0.62729799791922308,
'lens_offset_x' : [0, 0.0, 12.99, 12.99, 0.0, -12.99, -12.99, 0.0, 12.99, 25.98, 25.98, 25.98, 12.99, 0.0, -12.99,
                   -25.98, -25.98, -25.98, -12.99, 0.0, 12.99, 25.98, 38.97, 38.97, 38.97, 38.97, 25.98, 12.99, 0.0,
                   -12.99, -25.98, -38.97, -38.97, -38.97, -38.97, -25.98, -12.99, 0.0, 12.99, 25.98, 38.97,
                   51.962500000000006, 51.962500000000006, 51.962500000000006, 51.962500000000006, 51.962500000000006,
                   38.97, 25.98, 12.99, 0.0, -12.99, -25.98, -38.97, -51.962500000000006, -51.962500000000006,
                   -51.962500000000006, -51.962500000000006, -51.962500000000006, -38.97, -25.98, -12.99, 0.0, 12.99,
                   25.98, 38.97, 51.962500000000006, 64.9525, 64.9525, 64.9525, 64.9525, 64.9525, 64.9525,
                   51.962500000000006, 38.97, 25.98, 12.99, 0.0, -12.99, -25.98, -38.97, -51.962500000000006, -64.9525,
                   -64.9525, -64.9525, -64.9525, -64.9525, -64.9525, -51.962500000000006, -38.97, -25.98, -12.99, 25.98,
                   38.97, 51.962500000000006, 77.9425, 77.9425, 77.9425, 51.962500000000006, 38.97, 25.98, -25.98,
                   -38.97, -51.962500000000006, -77.9425, -77.9425, -77.9425, -51.962500000000006, -38.97, -25.98],
'lens_offset_y' : [0, -15.0, -7.5, 7.5, 15.0, 7.5, -7.5, -30.0, -22.5, -15.0, 0.0, 15.0, 22.5, 30.0, 22.5, 15.0, 0.0,
                   -15.0, -22.5, -45.0, -37.5, -30.0, -22.5, -7.5, 7.5, 22.5, 30.0, 37.5, 45.0, 37.5, 30.0, 22.5, 7.5,
                   -7.5, -22.5, -30.0, -37.5, -60.0, -52.5, -45.0, -37.5, -30.0, -15.0, 0.0, 15.0, 30.0, 37.5, 45.0,
                   52.5, 60.0, 52.5, 45.0, 37.5, 30.0, 15.0, 0.0, -15.0, -30.0, -37.5, -45.0, -52.5, -75.0, -67.5,
                   -60.0, -52.5, -45.0, -37.5, -22.5, -7.5, 7.5, 22.5, 37.5, 45.0, 52.5, 60.0, 67.5, 75.0, 67.5, 60.0,
                   52.5, 45.0, 37.5, 22.5, 7.5, -7.5, -22.5, -37.5, -45.0, -52.5, -60.0, -67.5, -75.0, -67.5, -60.0,
                   -15.0, 0.0, 15.0, 60.0, 67.5, 75.0, 75.0, 67.5, 60.0, 15.0, 0.0, -15.0, -60.0, -67.5, -75.0],
'DAQFreq' : 80000.0,
'FSPFreq' : 64.0,
'BSSE_TimeSpanRate' : 10.0,
'LSE_TimeSpanRate' : 0.55,
'mkidBaseShift' : -0.02,
'mkidShiftErr' : 4.0,
'mkidShiftErr_useShift' : 0.05,
'mkidShift_costumS' : 1.0*(64-3.6)/(1221-3.6),
'mkidShift_costumE' : 1.0*(72 - 3.6)/(1221-3.6),
'mkidShiftStep' : 0.001,
#'mkidShiftStep_useShift' : 0.0005,
#'ConvolStep' : 0.02,
'MaskRate' : 0.025,
'mkidForceCount' : 10,
'rejectRate' : 2,

'force' : True,
'use_shift' : True,
'use_costumSE' : True,
'shift_withMask' : True,

'fit_mr' : True,
'mr_s' : 0.02,
'mr_e' : 0.05,
'mr_a' : 0.001,
'rms_s' : 1.0*(90-3.6)/(1221-3.6),
'rms_e' : 1.0*(130 - 3.6)/(1221-3.6),
}
