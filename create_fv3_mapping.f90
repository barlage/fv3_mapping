program create_fv3_mapping

  use netcdf
  implicit none
  
  integer, parameter :: source_i_size         = 1024
  integer, parameter :: source_j_size         = 1024
  logical            :: source_read           = .true.   ! true if you will read from file, false then set below
  character*100      :: source_filename       = "/Users/barlage/work/data/ims/ims_latlon.nc"
  character*100      :: source_latname        = "ims_lat"
  character*100      :: source_lonname        = "ims_lon"
  logical            :: include_source_latlon = .true.
  logical            :: perturb_source_latlon = .true.   ! if lat/lon not found, then add a small value to nudge off boundary
  real, parameter    :: perturb_value         = 1.d-4    ! a small adjustment to lat/lon to find [radians]
  integer, parameter :: fv3_size = 96
  integer, parameter :: fv3_grid = fv3_size*2 + 1
  character*100      :: fv3_path = "/Users/barlage/work/data/C96.mx100_frac/"

  integer :: fv3_search_order(6) = (/3,1,2,5,6,4/)
  integer :: quick_search_pad = 1

  real   , dimension(source_i_size,source_j_size) :: source_lat, source_lon
  real   , dimension(fv3_grid,fv3_grid,6)         :: fv3_lat, fv3_lon

  integer, dimension(source_i_size,source_j_size) :: lookup_tile, lookup_i, lookup_j
  
  real, dimension(4) :: lat_vertex, lon_vertex
  
  integer :: itile, tile_index, tile_i_index, tile_j_index, source_i_index, source_j_index
  integer :: tile_save, tile_i_save, tile_j_save, pad_i_min, pad_i_max, pad_j_min, pad_j_max
  logical :: found, inside_a_polygon
  real    :: lat2find, lon2find
  integer :: ncid, dimid, varid, status   ! netcdf identifiers
  integer :: dim_id_i, dim_id_j           ! netcdf dimension identifiers
  character*100 :: filename
  real, parameter :: deg2rad = 3.1415926535897931/180.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read source lat/lon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(source_read) then

    status = nf90_open(source_filename, NF90_NOWRITE, ncid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, source_latname, varid)
    status = nf90_get_var(ncid, varid , source_lat)
  
    status = nf90_inq_varid(ncid, source_lonname, varid)
    status = nf90_get_var(ncid, varid , source_lon)

    status = nf90_close(ncid)

  else

    do source_j_index = 1, source_j_size
    do source_i_index = 1, source_i_size
      source_lat(source_i_index,source_j_index) = -90.d0 + 0.5d0*(source_j_index-1)
      source_lon(source_i_index,source_j_index) =   0.d0 + 0.5d0*(source_i_index-1)
    end do
    end do

  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read FV3 tile information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do itile = 1, 6

    if(fv3_size < 100) then
      write(filename,'(a1,i2,a10,i1,a3)') "C", fv3_size, "_grid.tile", itile, ".nc"
    elseif(fv3_size < 1000) then
      write(filename,'(a1,i3,a10,i1,a3)') "C", fv3_size, "_grid.tile", itile, ".nc"
    elseif(fv3_size < 10000) then
      write(filename,'(a1,i4,a10,i1,a3)') "C", fv3_size, "_grid.tile", itile, ".nc"
    else
      stop "unknown fv3 size"
    end if

    filename = trim(fv3_path)//trim(filename)

    status = nf90_open(filename, NF90_NOWRITE, ncid)
      if (status /= nf90_noerr) call handle_err(status)

    status = nf90_inq_varid(ncid, "x", varid)
    status = nf90_get_var(ncid, varid , fv3_lon(:,:,itile))
  
    status = nf90_inq_varid(ncid, "y", varid)
    status = nf90_get_var(ncid, varid , fv3_lat(:,:,itile))
  
    status = nf90_close(ncid)

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through the source points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  lookup_tile = -9999
  lookup_i = -9999
  lookup_j = -9999
  
  tile_save = fv3_search_order(1)
  tile_i_save = 1
  tile_j_save = 1
  
  source_i_loop : do source_i_index = 1, source_i_size
  source_j_loop : do source_j_index = 1, source_j_size
  
    found = .false.
    lat2find = source_lat(source_i_index, source_j_index)
    lon2find = source_lon(source_i_index, source_j_index)
    
    if(lat2find < -90. .or. lat2find > 90.  .or. &
       lon2find <   0. .or. lon2find > 360.) cycle source_j_loop     ! skip if out of projection
    
    lat2find = deg2rad * lat2find
    lon2find = deg2rad * lon2find
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! check around the last found tile/i/j
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    pad_i_min = max(tile_i_save-quick_search_pad,1)
    pad_i_max = min(tile_i_save+quick_search_pad,fv3_size)
    pad_j_min = max(tile_j_save-quick_search_pad,1)
    pad_j_max = min(tile_j_save+quick_search_pad,fv3_size)
    
    tile_index = tile_save
    
    do tile_i_index = pad_i_min, pad_i_max
    do tile_j_index = pad_j_min, pad_j_max
      
      lat_vertex(1) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
      lat_vertex(2) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
      lat_vertex(3) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
      lat_vertex(4) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
      
      lon_vertex(1) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
      lon_vertex(2) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
      lon_vertex(3) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
      lon_vertex(4) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
        
      lat_vertex = lat_vertex * deg2rad
      lon_vertex = lon_vertex * deg2rad
      
      found = inside_a_polygon(lon2find, lat2find, 4, lon_vertex, lat_vertex)
        
      if(found) then
        lookup_tile(source_i_index,source_j_index) = tile_index
        lookup_i   (source_i_index,source_j_index) = tile_i_index
        lookup_j   (source_i_index,source_j_index) = tile_j_index
        tile_save = tile_index
        tile_i_save = tile_i_index
        tile_j_save = tile_j_index
        cycle source_j_loop
      end if
        
    end do
    end do
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! not found so do a general check
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print*, "Did not find, doing general search"
    
    do itile = 1, 6

      tile_index = fv3_search_order(itile)
      
      do tile_i_index = 1, fv3_size
      do tile_j_index = 1, fv3_size
      
        lat_vertex(1) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
        lat_vertex(2) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
        lat_vertex(3) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
        lat_vertex(4) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
      
        lon_vertex(1) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
        lon_vertex(2) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
        lon_vertex(3) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
        lon_vertex(4) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
        
        lat_vertex = lat_vertex * deg2rad
        lon_vertex = lon_vertex * deg2rad
      
        found = inside_a_polygon(lon2find, lat2find, 4, lon_vertex, lat_vertex)
        
        if(found) then
          lookup_tile(source_i_index,source_j_index) = tile_index
          lookup_i   (source_i_index,source_j_index) = tile_i_index
          lookup_j   (source_i_index,source_j_index) = tile_j_index
          tile_save = tile_index
          tile_i_save = tile_i_index
          tile_j_save = tile_j_index
          cycle source_j_loop
        end if
        
      end do
      end do
      
    end do
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! not found so do a general check with a perturbation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print*, "Did not find, add perturbation"
    
    lat2find = lat2find + perturb_value
    lon2find = lon2find + perturb_value
    
    do itile = 1, 6

      tile_index = fv3_search_order(itile)
      
      do tile_i_index = 1, fv3_size
      do tile_j_index = 1, fv3_size
      
        lat_vertex(1) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
        lat_vertex(2) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
        lat_vertex(3) = fv3_lat((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
        lat_vertex(4) = fv3_lat((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
      
        lon_vertex(1) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 1,tile_index)  ! LL
        lon_vertex(2) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 1,tile_index)  ! LR
        lon_vertex(3) = fv3_lon((tile_i_index - 1) * 2 + 3,(tile_j_index - 1) * 2 + 3,tile_index)  ! UR
        lon_vertex(4) = fv3_lon((tile_i_index - 1) * 2 + 1,(tile_j_index - 1) * 2 + 3,tile_index)  ! UL
        
        lat_vertex = lat_vertex * deg2rad
        lon_vertex = lon_vertex * deg2rad
      
        found = inside_a_polygon(lon2find, lat2find, 4, lon_vertex, lat_vertex)
        
        if(found) then
          lookup_tile(source_i_index,source_j_index) = tile_index
          lookup_i   (source_i_index,source_j_index) = tile_i_index
          lookup_j   (source_i_index,source_j_index) = tile_j_index
          tile_save = tile_index
          tile_i_save = tile_i_index
          tile_j_save = tile_j_index
          cycle source_j_loop
        end if
        
      end do
      end do
      
    end do
    
    if(.not.found) then
      print*, "Did not find in cube sphere:", source_lat(source_i_index, source_j_index), ",", source_lon(source_i_index, source_j_index)
      stop
    end if

  end do source_j_loop
     if(mod(source_i_index,10) == 0) print *, "finished loop: ",source_i_index, " of ", source_i_size
  end do source_i_loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create the output filename and netcdf file (overwrite old)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(fv3_size < 100) then
    write(filename,'(a13,i2,a3)') "fv3_mapping_C", fv3_size, ".nc"
  elseif(fv3_size < 1000) then
    write(filename,'(a13,i3,a3)') "fv3_mapping_C", fv3_size, ".nc"
  elseif(fv3_size < 10000) then
    write(filename,'(a13,i4,a3)') "fv3_mapping_C", fv3_size, ".nc"
  else
    stop "unknown fv3 size"
  end if

  status = nf90_create(filename, NF90_CLOBBER, ncid)
    if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.

  status = nf90_def_dim(ncid, "idim"   , source_i_size , dim_id_i)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "jdim"   , source_j_size , dim_id_j)
    if (status /= nf90_noerr) call handle_err(status)
  
! Define variables in the file.

  status = nf90_def_var(ncid, "tile", NF90_INT, (/dim_id_i, dim_id_j/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "fv3 tile location")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "tile_i", NF90_INT, (/dim_id_i, dim_id_j/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "fv3 i location in tile")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "tile_j", NF90_INT, (/dim_id_i, dim_id_j/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "fv3 j location in tile")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999)
      if (status /= nf90_noerr) call handle_err(status)

 if(include_source_latlon) then

  status = nf90_def_var(ncid, "source_lat", NF90_FLOAT, (/dim_id_i, dim_id_j/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "source latitude")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999.)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "source_lon", NF90_FLOAT, (/dim_id_i, dim_id_j/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "source longitude")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999.)
      if (status /= nf90_noerr) call handle_err(status)
  
 end if

  status = nf90_enddef(ncid)

  status = nf90_inq_varid(ncid, "tile", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lookup_tile)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "tile_i", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lookup_i)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "tile_j", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lookup_j)
    if (status /= nf90_noerr) call handle_err(status)

 if(include_source_latlon) then

  status = nf90_inq_varid(ncid, "source_lat", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , source_lat)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "source_lon", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , source_lon)
    if (status /= nf90_noerr) call handle_err(status)
  
 end if

 status = nf90_close(ncid)

end program

  subroutine latlon2xyz(siz,lon, lat, x, y, z)
  implicit none
  integer, intent(in) :: siz
  real, intent(in) :: lon(siz), lat(siz)
  real, intent(out) :: x(siz), y(siz), z(siz)
  
  integer n

  do n = 1, siz
    x(n) = cos(lat(n))*cos(lon(n))
    y(n) = cos(lat(n))*sin(lon(n))
    z(n) = sin(lat(n))
  enddo

  end subroutine latlon2xyz

  FUNCTION spherical_angle(v1, v2, v3)
    implicit none
    real, parameter :: EPSLN30 = 1.e-30
    real, parameter :: PI=3.1415926535897931
    real v1(3), v2(3), v3(3)
    real  spherical_angle
 
    real px, py, pz, qx, qy, qz, ddd;
  
    ! vector product between v1 and v2 
    px = v1(2)*v2(3) - v1(3)*v2(2)
    py = v1(3)*v2(1) - v1(1)*v2(3)
    pz = v1(1)*v2(2) - v1(2)*v2(1)
    ! vector product between v1 and v3 
    qx = v1(2)*v3(3) - v1(3)*v3(2);
    qy = v1(3)*v3(1) - v1(1)*v3(3);
    qz = v1(1)*v3(2) - v1(2)*v3(1);

    ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz);
    if ( ddd <= 0.0 ) then
      spherical_angle = 0. 
    else 
      ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd);
      if( abs(ddd-1) < EPSLN30 ) ddd = 1;
      if( abs(ddd+1) < EPSLN30 ) ddd = -1;
      if ( ddd>1. .or. ddd<-1. ) then
    !FIX to correctly handle co-linear points (angle near pi or 0) */
    if (ddd < 0.) then
      spherical_angle = PI
    else
      spherical_angle = 0.
    endif
      else
    spherical_angle = acos( ddd )
      endif
    endif  

    return

  END FUNCTION spherical_angle

  FUNCTION inside_a_polygon(lon1, lat1, npts, lon2, lat2)
    implicit none
    real, parameter :: EPSLN10 = 1.e-10
    real, parameter :: EPSLN8 = 1.e-8
    real, parameter :: PI=3.1415926535897931
    real, parameter :: RANGE_CHECK_CRITERIA=0.05
    real :: anglesum, angle, spherical_angle
    integer i, ip1
    real lon1, lat1
    integer npts
    real lon2(npts), lat2(npts)
    real x2(npts), y2(npts), z2(npts)
    real lon1_1d(1), lat1_1d(1)
    real x1(1), y1(1), z1(1)
    real pnt0(3),pnt1(3),pnt2(3)
    logical inside_a_polygon
    real max_x2,min_x2,max_y2,min_y2,max_z2,min_z2
    !first convert to cartesian grid */
    call latlon2xyz(npts,lon2, lat2, x2, y2, z2);
    lon1_1d(1) = lon1
    lat1_1d(1) = lat1
    call latlon2xyz(1,lon1_1d, lat1_1d, x1, y1, z1);
    inside_a_polygon = .false.
    max_x2 = maxval(x2)
    if( x1(1) > max_x2+RANGE_CHECK_CRITERIA ) return
    min_x2 = minval(x2)
    if( x1(1)+RANGE_CHECK_CRITERIA < min_x2 ) return
    max_y2 = maxval(y2)
    if( y1(1) > max_y2+RANGE_CHECK_CRITERIA ) return
    min_y2 = minval(y2)
    if( y1(1)+RANGE_CHECK_CRITERIA < min_y2 ) return
    max_z2 = maxval(z2)
    if( z1(1) > max_z2+RANGE_CHECK_CRITERIA ) return
    min_z2 = minval(z2)
    if( z1(1)+RANGE_CHECK_CRITERIA < min_z2 ) return

    pnt0(1) = x1(1)
    pnt0(2) = y1(1)
    pnt0(3) = z1(1)
    
    anglesum = 0;
    do i = 1, npts
       if(abs(x1(1)-x2(i)) < EPSLN10 .and. &
          abs(y1(1)-y2(i)) < EPSLN10 .and. &
          abs(z1(1)-z2(i)) < EPSLN10 ) then ! same as the corner point
      inside_a_polygon = .true.
      return
       endif
       ip1 = i+1
       if(ip1>npts) ip1 = 1
       pnt1(1) = x2(i)
       pnt1(2) = y2(i)
       pnt1(3) = z2(i)
       pnt2(1) = x2(ip1)
       pnt2(2) = y2(ip1)
       pnt2(3) = z2(ip1)

       angle = spherical_angle(pnt0, pnt2, pnt1);
!       anglesum = anglesum + spherical_angle(pnt0, pnt2, pnt1);
       anglesum = anglesum + angle
    enddo

    if(abs(anglesum-2*PI) < EPSLN8) then
       inside_a_polygon = .true.
    else
       inside_a_polygon = .false.
    endif

    return
    
  end FUNCTION inside_a_polygon

  subroutine handle_err(status)
    use netcdf
    integer, intent ( in) :: status
 
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine handle_err
