Code to create file with indexes used to map from ASCII IMS file onto FV3 grid. 

Must be run for each FV3 resolution  - need to change resolution (2 places) in create_fv3_mapping.f90. 

Input (hard-coded into create_fv3_mapping.f90): 
-ASCII IMS example input file
-IMS lat lon files (created by Youlong). Note that the converted data has changed lat from its original
north-south to south-north. Such a change leads to a consistentency with the original ims ascii
file (south-north).
-FV3 grid tiles. 

To compile on hera, use build.sh. 

requires following modules: 

module load intel
module load netcdf/4.7.0

To run, submit executable (in big-mem queue)

There is a testcase and example submission script for C48 on hera at: 
/scratch2/BMC/gsienkf/Clara.Draper/DA_test_cases/snow/IMSobsproc/get_index/

Mike Barlage, Clara Draper, Youlong Xia.
