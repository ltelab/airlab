Initialization
==============
In matlab, navigate to the "airlab" directory
Run airpaths.m (add dir and subdirs to matlab path)


Image processing and calculation of particle desriptors
=======================================================
1) open PROCESS_OAP
2) set the correct path to campaigndir (where raw images are located), outputdir_mat, and outputdir_img
3) select the probe (2DS, HVPS or CPI)
4) select the raw image format ('image' for .png or .jpg; 'matfile' for .mat)
5) if your matlab version and computer allows, you can set process.parallel = true in order to run the code in parallel on your processor cores.
6) you can modify the other user parameters if desired but I would recommend to keep them as they are for best classification results.
7) run PROCESS_OAP


Snow habit classification
=========================
1) run predict_snow_habit(dir_data,probe,save_results,save_validation_files,truncated_threshold) with the following arguments:
	- dir_data = path to the processed images
	- probe = '2DS','HVPS' or 'CPI'
	- save_results = true
	- save_validation_files = true
	- threshold = empty or 0.3

save_validation_files is a new feature implemented to make visual check of the classification easier. If enabled, the program will save validation images displaying the cropped snowflake with its label and the associated probability. The images will be saved in dir_data/valid

