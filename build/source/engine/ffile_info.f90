! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module ffile_info_module
USE nrtype
USE netcdf
implicit none
private
public::ffile_info
contains


 ! ************************************************************************************************
 ! public subroutine ffile_info: read information on model forcing files
 ! ************************************************************************************************
 subroutine ffile_info(nHRU,err,message)
 ! used to read metadata on the forcing data file
 USE ascii_util_module,only:file_open
 USE netcdf_util_module,only:nc_file_open    ! open netCDF file
 USE summaFileManager,only:SETNGS_PATH       ! path for metadata files
 USE summaFileManager,only:INPUT_PATH        ! path for forcing files
 USE summaFileManager,only:FORCING_FILELIST  ! list of model forcing files
 USE data_struc,only:time_meta,forc_meta     ! model forcing metadata
 USE data_struc,only:forcFileInfo,data_step  ! info on model forcing file
 USE var_lookup,only:iLookTYPE               ! named variables to index elements of the data vectors
 USE get_ixname_module,only:get_ixtime,get_ixforce  ! identify index of named variable
 USE ascii_util_module,only:get_vlines      ! get a vector of non-comment lines
 USE ascii_util_module,only:split_line      ! split a line into words
 implicit none
 ! define input & output
 integer(i4b),intent(in)              :: nHRU           ! number of hydrologic response units
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 !netcdf file i/o related
 integer(i4b)                         :: ncid            ! netcdf file id
 integer(i4b)                         :: mode            ! netCDF file open mode
 integer(i4b)                         :: varid           ! netcdf variable id
 integer(i4b)                         :: dimId           ! netcdf dimension id
 character(LEN=nf90_max_name)         :: varName        ! character array of netcdf variable name
 integer(i4b)                         :: iNC             ! index of a variable in netcdf file
 integer(i4b)                         :: nvar            ! number of variables in netcdf local attribute file
 !the rest

 character(LEN=1024),allocatable      :: dataLines(:)   ! vector of lines of information (non-comment lines)
 integer(i4b),parameter               :: imiss = -999   ! missing data
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b),parameter               :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b),parameter               :: maxLines=1000000  ! maximum lines in the file
 character(LEN=256)                   :: filenameDesc   ! name of file that describes the forcing datafile

 integer(i4b)                         :: ivar           ! index of model variable
 integer(i4b)                         :: iFile
 integer(i4b)                         :: nFile
 integer(i4b)                         :: file_nHRU      ! number of HRUs in current forcing file
 real(dp)                             :: dataStep_iFile ! data step for a given forcing data file
 ! Start procedure here
 err=0; message="ffile_info/"
 ! ------------------------------------------------------------------------------------------------------------------
 ! (1) read from the list of forcing files
 ! ------------------------------------------------------------------------------------------------------------------
  ! build filename for forcing file list
 infile = trim(SETNGS_PATH)//trim(FORCING_FILELIST)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get a list of character strings from non-comment lines
 call get_vlines(unt,dataLines,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
 nFile = size(dataLines)
 ! allocate space for forcing information
 if(associated(forcFileInfo)) deallocate(forcFileInfo)
 allocate(forcFileInfo(nFile), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for forcFileInfo'; return; endif
 do iFile=1,nFile
! split the line into "words" (expect one word: the file describing forcing data for that index)
  read(dataLines(iFile),*,iostat=err) filenameDesc
  if(err/=0)then; message=trim(message)//'problem reading a line of data from file ['//trim(infile)//']'; return; endif
  forcFileInfo(iFile)%filenmDesc = trim(filenameDesc)
 enddo  ! (looping through files)
 ! ------------------------------------------------------------------------------------------------------------------
 ! (2) pull descriptive information from netcdf forcing file and check number of HRUs in each forcing file matches nHRU
 ! ------------------------------------------------------------------------------------------------------------------
 ! check that the time metadata is already populated
 if(.not.associated(time_meta))then; err=30; message=trim(message)//"TimeMetadataNonexistent"; return; endif
 ! check that the forcing metadata is already populated
 if(.not.associated(forc_meta))then; err=30; message=trim(message)//"ForcingMetadataNonexistent"; return; endif
 ! read
 do iFile=1,nFile
  ! allocate space for the column indices
  if(associated(forcFileInfo(iFile)%time_ix)) deallocate(forcFileInfo(iFile)%time_ix)
  if(associated(forcFileInfo(iFile)%data_ix)) deallocate(forcFileInfo(iFile)%data_ix)
  allocate(forcFileInfo(iFile)%time_ix(size(time_meta)),&
           forcFileInfo(iFile)%data_ix(size(forc_meta)),stat=err)
  if(err/=0)then; err=40; message=trim(message)//"problemAllocateStructureElement"; return; endif
  ! initialize column indices to missing
  forcFileInfo(iFile)%time_ix(:) = imiss
  forcFileInfo(iFile)%data_ix(:) = imiss
  ! build filename bfor actual forcing file
  infile = trim(INPUT_PATH)//trim(forcFileInfo(iFile)%filenmDesc)
  ! open file
  mode=nf90_NoWrite
  call nc_file_open(trim(infile), mode, ncid, err, cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  !how many variables are there?
  err = nf90_inquire(ncid, nvariables=nVar)
  call netcdf_err(err,message); if (err/=0) return
  !set forcing file name and nvar attributes
  forcFileInfo(iFile)%filenmData = trim(forcFileInfo(iFile)%filenmDesc)
  forcFileInfo(iFile)%ncols = nVar
  !loop through variables in netcdf file, check to see if everything needed to run the model exists and data_step is correct
  do iNC=1,nvar
   !inqure about current variable name, type, number of dimensions
   err = nf90_inquire_variable(ncid,iNC,name=varName)
   if(err/=0)then; message=trim(message)//'problem inquiring variable: '//trim(varName); return; endif
   !inquire nhru dimension size
   err = nf90_inq_dimid(ncid,'nHRU',dimId)
   if(err/=0)then; message=trim(message)//'nHRU dimension not present'; return; endif
   err = nf90_inquire_dimension(ncid,dimId,len=file_nHRU)
   if(file_nHRU /= nHRU)then; err=20;write(message,'(a,i0,a,i0)') trim(message)//'File HRU dimension: ',file_nHRU ,'not equal to nHRU: ',nHRU; return; endif
print *,'current var',iNC,trim(varName)
   select case(trim(varName))
    case('data_step' )
     ! read data_step from netcdf file
     err = nf90_inq_varid(ncid, "data_step", varId)
     if(err/=0)then; message=trim(message)//'problem inquiring data_step'; return; endif
     err = nf90_get_var(ncid,varid,dataStep_iFile)
     if(err/=0)then; message=trim(message)//'problem reading data_step'; return; endif
     !check data_step
     if(iFile == 1)then
      data_step = dataStep_iFile
print *,'data step!',data_step
     else
      if(abs(dataStep_iFile - data_step) > epsilon(dataStep_iFile))then
       write(message,'(a,i0,a)') trim(message)//'data step for forcing file ',iFile,'differs from the datastep of the first forcing file'
       err=20; return
      endif
     endif
    !forcing variables
    case('time','pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum')
     ivar = get_ixforce(trim(varname))
     if(ivar < 0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return; endif
     if(ivar>size(forcFileInfo(iFile)%data_ix))then
      err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
     endif
     ! put netcdf file variable index in the structure
     err = nf90_inq_varid(ncid, trim(varName), forcFileInfo(iFile)%data_ix(ivar))
     if(err/=0)then; message=trim(message)//"problem inquiring forcing variable[var="//trim(varName)//"]"; return; endif
    !for lat,lon,hruId do nothing
    case('hruId','latitude','longitude')
    ! check that variables are what we expect
    case default
     message=trim(message)//'unknown variable ['//trim(varName)//'] in local attributes file'
     err=20; return
   end select

  enddo ! (end of netecdf file variable loop)
 end do ! (loop through lines in the file)


 ! identify the first HRU to use a given data file
! do iFile=1,nHRU
!  forcFileInfo(iHRU)%ixFirstHRU = 0
!  do jHRU=1,iHRU-1
!   if(trim(forcFileInfo(iHRU)%filenmData) == trim(forcFileInfo(jHRU)%filenmData))then
!    forcFileInfo(iHRU)%ixFirstHRU = jHRU  ! index of first HRU to share the same data
!   endif
!  end do
! end do

 end subroutine ffile_info





!  call file_open(trim(infile),unt,err,cmessage)
!  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! get to the start of the variable descriptions
!  do iline=1,maxLines
!   read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit    ! read line of data
!   if (temp(1:1)/='!') exit  ! assume first line not comment is format code
!  end do ! looping through file to find the format code
!  ! read in format string
!  read(temp,*)ffmt
!  ! loop through the lines in the file
!  do iline=1,maxLines
!   ! read a line of data and exit if an error code (character read, so only possible error is end of file)
!   read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit
!   ! check that the line is not a comment
!   if (temp(1:1)=='!')cycle
!   ! save data into a temporary variables
!   read(temp,trim(ffmt),iostat=err) varname, dLim, vardata
!   if (err/=0) then; err=30; message=trim(message)//"errorReadLine[file="//trim(infile)//"; line="//trim(temp)//"]"; return; endif
!   ! check the delimiter
!   if(dLim(1:1)/='|')then; err=30; message=trim(message)//"incorrectFormat"//trim(infile); return; endif
!   !print*, 'varname = ', trim(varname)
!   !print*, 'vardata = ', trim(vardata)
!   ! put data into data structure
!   select case(trim(varname))
!    case('filenmData'); read(vardata,*) forcFileInfo(iHRU)%filenmData
!    case('ncols'     ); read(vardata,*) forcFileInfo(iHRU)%ncols
!    ! process the data step
!    case('data_step' )
!     read(vardata,*) dataStep_iHRU
!     if(iHRU == 1)then
!      data_step = dataStep_iHRU
!     else
!      if(abs(dataStep_iHRU - data_step) > epsilon(dataStep_iHRU))then
!       write(message,'(a,i0,a)') trim(message)//'data step for HRU ',iHRU,'differs from the datastep of the first HRU'
!       err=20; return
!      endif
!     endif
!    ! ***** identify the index of the time data variable
!    case('iyyy','im','id','ih','imin')
!     ivar = get_ixtime(trim(varname))
!     if(ivar < 0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return; endif
!     if(ivar>size(forcFileInfo(iHRU)%time_ix))then
!      err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
!     endif
!     ! put column index in the structure
!     read(vardata,*) forcFileInfo(iHRU)%time_ix(ivar)
!    ! ***** identity index for the forcing data variable
!    case('pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum')
!     ivar = get_ixforce(trim(varname))
!     if(ivar < 0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return; endif
!     if(ivar>size(forcFileInfo(iHRU)%data_ix))then
!      err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
!     endif
!     ! put column index in the structure
!     read(vardata,*) forcFileInfo(iHRU)%data_ix(ivar)
!    ! ***** error check
!    case default
!     message=trim(message)//'variableNotFound[var='//trim(varname)//'; file='//trim(infile)//']'
!     err=20; return
!   endselect
!  enddo ! (loop through lines in the file)
!  ! close file unit
!  close(unt)
! end do  ! (looping through files describing each HRU)
! ! identify the first HRU to use a given data file
! do iHRU=1,nHRU
!  forcFileInfo(iHRU)%ixFirstHRU = 0
!  do jHRU=1,iHRU-1
!   if(trim(forcFileInfo(iHRU)%filenmData) == trim(forcFileInfo(jHRU)%filenmData))then
!    forcFileInfo(iHRU)%ixFirstHRU = jHRU  ! index of first HRU to share the same data
!   endif
!  end do
! end do
! end subroutine ffile_info


! **********************************************************************************************************
 ! private subroutine netcdf_err: error control
 ! **********************************************************************************************************
 subroutine netcdf_err(err,message)
 ! used to handle errors for NetCDF calls
 use netcdf
 implicit none
 ! declare dummies
 integer(i4b), intent(inout)   :: err
 character(*), intent(inout)   :: message
 ! start procedure here
 if (err/=nf90_noerr) then
  message=trim(message)//"["//trim(nf90_strerror(err))//"]"
  err=200
 endif
 end subroutine netcdf_err


end module ffile_info_module
