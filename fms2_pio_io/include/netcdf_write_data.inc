!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @file
!> @brief Routines for writing variable values to netcdf files with different
!! dimension sizes for the @ref netcdf_write_data interface

!> @addtogroup netcdf_io_mod
!> @{

!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_0d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(1) :: c
  integer :: ndims
  integer, dimension(:), allocatable :: start
  integer, dimension(:), allocatable :: dimsizes
  character, dimension(:), allocatable :: charbuf
  integer :: i
  integer :: tlen
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message
  
  append_error_msg = "netcdf_write_data_0d: file:"//trim(fileobj%path)//" variable: "//trim(variable_name)
  !if (fileobj%is_root) then !--- pio write routines are to be called by all PEs and not just root PE. 
    c(:) = 1
    if (present(corner)) then
      c(1) = corner
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 1) then
        call error("unlimited dimension must be the slowest varying dimension: "//trim(append_error_msg))
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(variable_data)
      type is (integer(kind=i4_kind))
        if (present(unlim_dim_level)) then
          err = pio_put_var(fileobj%file_desc, varid, start=c(1:1), count=[1], ival=[variable_data])
        else
          err = pio_put_var(fileobj%file_desc, varid, variable_data)
        endif
      type is (integer(kind=i8_kind))
        call error(trim(fileobj%path)//": 64 bit integers are not supported with PIO.")
      type is (real(kind=r4_kind))
        if (present(unlim_dim_level)) then
          err = pio_put_var(fileobj%file_desc, varid, start=c(1:1), count=[1], ival=[variable_data])
        else
          err = pio_put_var(fileobj%file_desc, varid, variable_data)
        endif
      type is (real(kind=r8_kind))
        if (present(unlim_dim_level)) then
          err = pio_put_var(fileobj%file_desc, varid, start=c(1:1), count=[1], ival=[variable_data])
        else
          err = pio_put_var(fileobj%file_desc, varid, variable_data)
        endif
      type is (character(len=*))
        call error(trim(fileobj%path)//": character write not implemented yet (netcdf_write_data.inc)")
        ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
        if (ndims .ne. 1) then
          call error("currently only scalar and 1d character writes are supported: "//trim(append_error_msg))
        endif
        allocate(start(ndims))
        start(:) = 1
        allocate(dimsizes(ndims))
        call get_variable_size(fileobj, variable_name, dimsizes, .false.)
        call allocate_array(charbuf, dimsizes)
        charbuf(:) = ""
        tlen = len_trim(variable_data)
        if (tlen .gt. dimsizes(1)) then
          call error("character buffer is too big; decrease length: "//trim(append_error_msg))
        endif
        do i = 1, tlen
          charbuf(i) = variable_data(i:i)
        enddo
        err = nf90_put_var(fileobj%ncid, varid, charbuf, start=start, count=dimsizes)
        deallocate(charbuf)
        deallocate(dimsizes)
        deallocate(start)
      class default
        call error("Unsupported variable type: "//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)
  !endif
end subroutine netcdf_write_data_0d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_1d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, dimension(1), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(1), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer :: ndims
  integer, dimension(:), allocatable :: start
  integer, dimension(:), allocatable :: dimsizes
  character, dimension(:,:), allocatable :: charbuf
  character(len=1024) :: sbuf
  integer :: i
  integer :: j
  integer :: tlen
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message
  integer :: true_ndim

  append_error_msg = "netcdf_write_data_1d: file:"//trim(fileobj%path)//" variable: "//trim(variable_name)

  !if (fileobj%is_root) then !--- pio write routines are to be called by all PEs and not just root PE. 
    c(:) = 1
    if (present(corner)) then
      c(1:1) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:1) = edge_lengths(:)
    else
      e(1:1) = shape(variable_data)
    endif
    true_ndim = 1
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      ! FMS doesn't count unlim (time) dimension towards the total number of variables dimension.
      ! Thus, when the variable has an unlimited dimension, we call PIO appropriately, i.e., by
      ! incrementing the number ndims to account for the time dimension.
      true_ndim = 2
      if (unlim_dim_index .ne. 2) then
        call error("unlimited dimension must be the slowest varying dimension: "//trim(append_error_msg))
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(variable_data)
      type is (integer(kind=i4_kind))
        err = pio_put_var(fileobj%file_desc, varid, start=c(1:true_ndim), count=e(1:true_ndim), ival=variable_data)
      type is (integer(kind=i8_kind))
        call error(trim(fileobj%path)//": 64 bit integers are not supported with PIO.")
      type is (real(kind=r4_kind))
        err = pio_put_var(fileobj%file_desc, varid, start=c(1:true_ndim), count=e(1:true_ndim), ival=variable_data)
      type is (real(kind=r8_kind))
        err = pio_put_var(fileobj%file_desc, varid, start=c(1:true_ndim), count=e(1:true_ndim), ival=variable_data)
      type is (character(len=*))
        call error(trim(fileobj%path)//": character write not implemented yet (netcdf_write_data.inc)")
        ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
        if (ndims .ne. 2) then
          call error("currently only scalar and 1d character writes are supported: "//trim(append_error_msg))
        endif
        allocate(start(ndims))
        start(:) = 1
        allocate(dimsizes(ndims))
        call get_variable_size(fileobj, variable_name, dimsizes, .false.)
        call allocate_array(charbuf, dimsizes)
        charbuf(:,:) = ""
        tlen = len(variable_data(1))
        if (tlen .gt. dimsizes(1)) then
          call error("character buffer is too big; decrease length: "//trim(append_error_msg))
        endif
        if (size(variable_data) .ne. dimsizes(2)) then
          call error("incorrect size of variable_data array: "//trim(append_error_msg))
        endif
        do j = 1, dimsizes(2)
          call string_copy(sbuf, variable_data(j))
          do i = 1, tlen
            charbuf(i,j) = sbuf(i:i)
          enddo
        enddo
        err = nf90_put_var(fileobj%ncid, varid, charbuf, start=start, count=dimsizes)
        deallocate(charbuf)
        deallocate(dimsizes)
        deallocate(start)
      class default
        call error("Unsupported variable type: "//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)
  !endif
end subroutine netcdf_write_data_1d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_2d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, dimension(2), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(2), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer,dimension(3) :: c
  integer, dimension(3) :: e
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  print *, "ERRORline", 234, "write_data.inc"; call error("PIO version not implemented!!!")

  append_error_msg = "netcdf_write_data_2d: file:"//trim(fileobj%path)//" variable: "//trim(variable_name)

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1:2) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:2) = edge_lengths(:)
    else
      e(1:2) = shape(variable_data)
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 3) then
        call error("unlimited dimension must be the slowest varying dimension: "//trim(append_error_msg))
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(variable_data)
      type is (integer(kind=i4_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (integer(kind=i8_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=r4_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=r8_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c ,count=e)
      class default
        call error("Unsupported variable type: "//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)
  endif
end subroutine netcdf_write_data_2d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_3d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(3), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(3), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(4) :: c
  integer, dimension(4) :: e
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  print *, "ERRORline", 301, "write_data.inc"; call error("PIO version not implemented!!!")

  append_error_msg = "netcdf_write_data_3d: file:"//trim(fileobj%path)//" variable: "//trim(variable_name)

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1:3) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:3) = edge_lengths(:)
    else
      e(1:3) = shape(variable_data)
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 4) then
        call error("unlimited dimension must be the slowest varying dimension: "//trim(append_error_msg))
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(variable_data)
      type is (integer(kind=i4_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (integer(kind=i8_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=r4_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=r8_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      class default
        call error("Unsupported variable type: "//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)
  endif
end subroutine netcdf_write_data_3d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_4d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, dimension(4), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(4), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer,dimension(5) :: c
  integer, dimension(5) :: e
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message
  integer :: ndim

  unlim_dim_index = -1
  c(:) = 1
  if (present(corner)) then
    c(1:4) = corner(:)
  endif
  e(:) = 1
  if (present(edge_lengths)) then
    e(1:4) = edge_lengths(:)
  else
    e(1:4) = shape(variable_data)
  endif
  if (present(unlim_dim_level)) then
    unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                             broadcast=.false.)
   if  (unlim_dim_index /= -1) then
     c(unlim_dim_index) = unlim_dim_level
   endif
  endif

  ndim = get_variable_num_dimensions(fileobj, variable_name)

  if (unlim_dim_index == -1) then 
    ! The field has NO time dimension
    if (ndim == 1) then
      if (size(variable_data,2)/=1 .or. size(variable_data,3)/=1 .or. size(variable_data,4)/=1) then
        call error("PIO data size assertion violation in netcdf_write_data_4d")
      endif
      call netcdf_write_data_1d(fileobj, variable_name, variable_data(:,1,1,1), unlim_dim_level, c(1), e(1))
      return
    elseif (ndim == 2) then
      if (size(variable_data,3)/=1 .or. size(variable_data,4)/=1) then
        call error("PIO data size assertion violation in netcdf_write_data_4d")
      endif
      call netcdf_write_data_2d(fileobj, variable_name, variable_data(:,:,1,1), unlim_dim_level, c(1:2), e(1:2))
      return
    elseif (ndim == 3) then
        if (size(variable_data,4)/=1) then
        call error("PIO data size assertion violation in netcdf_write_data_4d")
      endif
      call netcdf_write_data_3d(fileobj, variable_name, variable_data(:,:,:,1), unlim_dim_level, c(1:3), e(1:3))
      return
    endif
  else
    ! The field has time dimension
    if (ndim == 1) then
      if (size(variable_data,1)/=1 .or. size(variable_data,2)/=1 .or. size(variable_data,3)/=1 .or. size(variable_data,4)/=1) then
        call error("PIO data size assertion violation in netcdf_write_data_4d")
      endif
      call netcdf_write_data_0d(fileobj, variable_name, variable_data(1,1,1,1), unlim_dim_level, c(1))
      return
    elseif (ndim == 2) then
      if (size(variable_data,2)/=1 .or. size(variable_data,3)/=1 .or. size(variable_data,4)/=1) then
        call error("PIO data size assertion violation in netcdf_write_data_4d")
      endif
      call netcdf_write_data_1d(fileobj, variable_name, variable_data(:,1,1,1), unlim_dim_level, c(1:2), e(1:2))
      return
    elseif (ndim == 3) then
        if (size(variable_data,3)/=1 .or. size(variable_data,4)/=1) then
        call error("PIO data size assertion violation in netcdf_write_data_4d")
      endif
      call netcdf_write_data_2d(fileobj, variable_name, variable_data(:,:,1,1), unlim_dim_level, c(1:3), e(1:3))
      return
    endif
  endif
end subroutine netcdf_write_data_4d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_5d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in), target :: variable_data !< Data that will be
                                                                      !! written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(5), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(5), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer,dimension(6) :: c
  integer, dimension(6) :: e
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  print *, "ERRORline", 434, "write_data.inc"; call error("PIO version not implemented!!!")

  append_error_msg = "netcdf_write_data_5d: file:"//trim(fileobj%path)//" variable: "//trim(variable_name)

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1:5) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:5) = edge_lengths(:)
    else
      e(1:5) = shape(variable_data)
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 6) then
        call error("unlimited dimension must be the slowest varying dimension: "//trim(append_error_msg))
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(variable_data)
      type is (integer(kind=i4_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (integer(kind=i8_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=r4_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=r8_kind))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      class default
        call error("Unsupported variable type: "//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)
  endif
end subroutine netcdf_write_data_5d
