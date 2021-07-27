# MKID_Data_Analysis_Pipeline
The tools for analyzing MKID camera data

Each function contains the description of its parameters and what they do in their description. Kindly refer to 
the source code for detailed explanation of the code, calculations and the parameters used. The following is a step-by-step guide to data analysis using the new pipeline with basic instructions and explanation.

All the functions are in the file *data_reduction.py* int the class **BinaryConvert**. 

To start using the functions, follow the following steps

    import data_reduction as dr
    import image_synthesis as im
    
    bc = dr.BinaryConvert()  # Call the object

To convert the MKID readout files to binary pickle file format, 

    mkid_files = bc.select_files() # Select the MKID readout files from a folder.
    bc.save_fsp_pickle(mkid_files, "Mars_0601") # Average and save the data to "Mars_0601.pickle" file.

The above step needs to be done only once per data, We shall use the pickle file created to load the data into python.

To load the data we use

    data = bc.load_pickle("Mars_0601.pickle")

Now we shall check which pixels are to be used. Depending upon the noise and the correlation factor of the
pixels, not all pixels are good for noise removal algorithm that we are going to implement. We can do so by 
manually excluding the pixels we don't want to include or by using a noise level calculation and exclusion function.
To do it manually, we can plot each Time Order Data and note down the pixels we wish to exclude. simply run,

    for i in range(len(data)):
    bc.plot(data, i+1)

This will plot the TOD one by one, and the user can note the pixels to be excluded in a list as such,

    exclude = [48, 71, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 89, 90, 91, 92, 93, 94, 95, 98, 99]

Then to exclude these pixels from further calculations,

    bc.pixels_to_exclude(exclude) # to exclude these pixels from further calculations.

The next part is to identify and select the different portions of a typical observation data. While observing the 
hot load, the frequency shift drops significantly, and it is carried out before and after each scan of the source
of interest. The function identify_frequencies can perform this task by providing the required fields. Simply run

    freq = bc.identify_frequencies(data, 4000, exclude=True, plot=False, save_file=True)

The details of the parameters can be found as usual in the source code description of each function. Please tweak the
minima_order(here, the value 4000) to get proper minima detection, using *plot=True* and *save_file=False* at first. Then
change them back as shown above to save this analysis result as **frequencies.txt**. This will be necessary for later steps.

Now that the frequencies are identified, we can crop the observation data. To do so, run the command

    clipped = bc.clip(data, 16 * 255, -16 * 160) # clip for decorrelation

Here, clipped data would exclude the data using the flags from *identify_frequencies* function. Modilfy the last two parameters
to exclude the data properly.

There are several types of possible artifacts in the data. One of them is the sudden jump in the data. To remove such 
jumps, run the following command,

    jump_removed = bc.changepoint_remove(clipped, 1000, 17173)

Here, the value *1000* is the padding used to reallign the data, and the *17173*, is the index at which the jump occured.

---

The next step is very important for PCA method. We need to remove the median of the data. This process is called centerng
the data. To do so, run the command

    mean_sub = bc.mean_center(clipped)

Here we start the decorrelation analysis. First we do a simple PCA analysis by running

    res_pca = bc.pca_decor(mean_sub, n_components=3)

To see the result, plot the data using

    for i in range(10): 
    plt.plot(res_pca[i][:, 1], alpha = 0.5)

Now to use the result of the PCA analysed data for the next stage, we isolateand flatten the data using the command,

    flat = bc.flatten(res_pca, 0)

With this data available, we now apply the mask 

    mask = bc.get_sigmaclipped_mask(flat, avg_points=32, sigma=4, maxiters=4, axis=1)

The later stages of the chunkedPCA algorithm requires this mask to identify and execute the signal exclusion algorithm.
Now run,

    chunk_matrix = bc.make_chunk_matrix(mask, num_chunks=400)



