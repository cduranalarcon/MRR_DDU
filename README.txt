"MRR_processing.ipynb" is a jupyter notebook files to process and display data.

Requirements:

- Numpy
- pylab (similar to matplotlib, it is possible to call matplotlib instead, and modified the code accordingly)
- netCDF4 (for displaying) try with "pip install netCDF4" if it is not installed in cmd or look for an alternative whl files in https://www.lfd.uci.edu/~gohlke/pythonlibs/#netcdf4
- datetime
- copy
- warnings
- Last version of IMProToo is already downloaded in "...MRR\lib\IMProToo-master". No need for further installation. 

Output folders are created automatically (Plots, MK_processed, temp).

The used must define the name of the input directory and the first and last date of the dataset.

You can copy the "RawSpectra" file from the MRR directly into "...MRR\Data\Station_Name" 

There is no problem with gap of data.

Any question can be addressed to claudio.duran@univ-grenoble-alpes.fr 
