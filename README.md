# find-and-refine-affine
a ImageJ script that find affine transformation of landmarks, then refine the affine with elastix.
# input (see example folder)
* 5x_image: the czi file of the 5x image.
* 20x_image: the czi file of the 20x image.
* landmarks: the csv file of the landmarks exported from BigWarp.
* landmarks_channel: a txt file that specifies which channel of the image that the landmarks are based on.
* elastix script: the elastix.py file.
# output
* affine_lse.csv : saved in landmarks file's folder. The affine transformation matrix.
* affine_elastix.csv : saved in landmarks file's folder. The refined affine transformation matrix.
# python dependencies
itk, itk-elastix, numpy
# components
* main.groovy : the main script to run in ImageJ.
* jama-1.0.3.jar : a java package that supports matrix operations.
* elastix.py : the elastix script that was called by main.groovy to perform transformation refinement with elastix.
* temp : a folder that holds the temp files generated by the script.
