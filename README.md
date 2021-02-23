# LCMS_process
Script to process the results from the LCMS milk experiment

# imput for the script:

The imput file for the following should be as a xlsx file with a sheet for each sample. The sample_ID must be used to name the excel sheet. Each data sheet should contain at least the following columns (if more columns are in the file, they will be ignored by the script):

- Chromatogram	

- RT [min]	

- Area

The list of sample_ID to process should be given through a file called "subset_list.xlsx" with one column 'ID'. An example file is provided in this repo.

