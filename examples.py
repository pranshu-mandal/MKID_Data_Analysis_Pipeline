"""
To describe the workflow using some methods

** Intended to be used line by line
"""

import data_reduction as dr

bc = dr.BinaryConvert() # Call the object
# mkid_files = bc.select_files() # Select the MKID readout files
# bc.save_as_pickle(mkid_files, "Mars_0601.pickle") # to average and save the fsp files to a single binary file.

data = bc.load_pickle("Mars_0601.pickle")

bc.plot(data, 1)