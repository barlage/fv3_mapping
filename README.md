Code to create file with indexes used to map observation lat, lon  onto FV3 grid. Currently coded for SMAP and IMS.

Namelist options: 

 tile_dim - model resolution
 otype    - filename stub for orography files (something like "oro_C384.mx025") 
 orog_path - path to orography files
 tile_path - path to grid files
 obs_source - observation resolution (either "IMS4km" or "IMS24km" for IMS, or "SMAP9km"
 coord_path - path to file with observation lat lon files.

On hera: 
 smap_coord_path = "/scratch2/NCEPDEV/land/data/DA/soil_moisture/SMAP/fix_coords/"
(coordinate file is:  "NSIDC0772_LatLon_EASE2_M09km_v1.0.nc") 

 ims_coord_path =  "/scratch2/NCEPDEV/land/data/DA/snow_ice_cover/IMS/fix_coords/" 
(coordinate files are:  imslat_4km_8bytes.bin  imslon_4km_8bytes.bin for IMS4km)

-above IMS lat lon files are created by Youlong. Note that the converted data has changed lat from its original
north-south to south-north. Such a change leads to a consistentency with the original ims ascii
file (south-north).

To compile on hera: 
>source mods_bash 
>use build.sh. 

To submit:  
change account in submit_mapping to your account.
>sbatch submit_mapping.sh 

Mike Barlage, Clara Draper, Youlong Xia.
