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

module read_force_module
implicit none
private
public::read_force
contains


 ! ************************************************************************************************
 ! public subroutine read_force: read in forcing data
 ! ************************************************************************************************
 subroutine read_force(istep,iGRU,iHRU,iFile,iRead,ncid,err,message)
 ! provide access to subroutines
 USE netcdf                                            ! netcdf capability
 USE netcdf_util_module,only:nc_file_open              ! open netcdf file
 USE nrtype                                            ! variable types, etc.
 USE summaFileManager,only:INPUT_PATH                  ! path of the forcing data file
 USE time_utils_module,only:extractTime                ! extract time info from units string
 USE time_utils_module,only:compJulday                 ! convert calendar date to julian day
 USE time_utils_module,only:compcalday                 ! convert julian day to calendar date
 USE multiconst,only:secprday                          ! number of seconds in a day
 USE data_struc,only:time_meta                         ! time metadata
 USE data_struc,only:forcFileInfo                      ! forcing file info
 USE data_struc,only:data_step                         ! length of the data step (s)
 USE data_struc,only:dJulianStart                      ! julian day of start time of simulation
 USE data_struc,only:refTime,refJulday                 ! reference time
 USE data_struc,only:fracJulDay                        ! fractional julian days since the start of year
 USE data_struc,only:yearLength                        ! number of days in the current year
 USE data_struc,only:time_meta                         ! metadata structures
 USE data_struc,only:time_data                         ! time information
 USE data_struc,only:forc_data                         ! forcing data
 USE data_struc,only:gru_struc                         ! gru-hru mapping structure
 USE var_lookup,only:iLookTIME,iLookFORCE              ! named variables to define structure elements
 USE get_ixname_module,only:get_ixforce                ! identify index of named variable
 USE multiconst,only:integerMissing                    ! integer missing value
 implicit none
 ! define dummy variables
 integer(i4b),intent(in)           :: istep            ! time index AFTER the start index
 integer(i4b),intent(in)           :: iHRU             ! index of hydrologic response unit
 integer(i4b),intent(in)           :: iGRU             ! index of grouped response unit
 integer(i4b),intent(inout)        :: iFile            ! index of current forcing file in forcing file list
 integer(i4b),intent(inout)        :: iRead            ! index of read position in time dimension in current netcdf file
 integer(i4b),intent(inout)        :: ncid             ! netcdf file identifier
 integer(i4b),intent(out)          :: err              ! error code
 character(*),intent(out)          :: message          ! error message
 ! define local variables
 !netcdf related
 integer(i4b)                      :: varId             ! variable identifier
 integer(i4b)                      :: dimId             ! dimension identifier
 integer(i4b)                      :: mode              ! netcdf file mode
 integer(i4b)                      :: dimLen            ! dimension length
 integer(i4b)                      :: attLen            ! attribute length
 character(len = nf90_max_name)    :: varName           ! dimenison name
 integer(i4b)                      :: ncStart(2)        ! start array for reading hru forcing
 !rest
 real(dp),parameter                :: amiss= -1.d+30   ! missing real
 character(len=256)                :: infile           ! filename
 character(len=256)                :: cmessage         ! error message for downwind routine
 character(len=256)                :: refTimeString    ! reference time string
 logical(lgt)                      :: xist             ! .TRUE. if the file exists
 integer(i4b),parameter            :: baseUnit=28      ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                      :: iline            ! loop through lines in the file
 integer(i4b)                      :: iNC              ! loop through variables in forcing file
 integer(i4b)                      :: iVar             ! index of forcing variable in forcing data vector
 integer(i4b)                      :: iFFile           ! forcing file counter
 integer(i4b)                      :: nFile            ! number of forcing files
 integer(i4b)                      :: hruId            ! unique hru id
 real(dp)                          :: dsec             ! double precision seconds (not used)
 real(dp)                          :: juldayFirst      ! julian day of the first time step in the data file
 real(dp)                          :: startJulDay      ! julian day at the start of the year
 real(dp)                          :: currentJulday    ! Julian day of current time step
 logical(lgt),parameter            :: checkTime=.false.  ! flag to check the time
 logical(lgt)                      :: foundTime=.false.  ! flag to note if we've found the first simulation timestep
 real(dp)                          :: dataStepFracDay  ! fraction of day of data step
 real(dp)                          :: dataJulDay       ! julian day of current forcing data step being read
 real(dp)                          :: tempVar          ! temporary floating point variable
 real(dp),allocatable              :: fileTime(:)      ! array of time from netcdf file
 real(dp),allocatable              :: diffTime(:)      ! array of time differences

 ! Start procedure here
 err=0; message="read_force/"
 !determine the julDay of current model step (istep) we need to read
 if(istep==1)then
  currentJulDay = dJulianStart
 else
  currentJulDay = dJulianStart + (data_step*real(iStep-1,dp))/secprday
 endif
 !fraction of day for data_step
 dataStepFracDay = data_step/secprday

! **********************************************************************************************
! ***** part 0: if file not open, then open file and find initial time in file
! **********************************************************************************************
 if(ncid==integerMissing)then
  ! define the reference time for the model simulation
  ! get attribute from time variable
  iFile=1
  ! define filename
  infile=trim(INPUT_PATH)//trim(forcFileInfo(iFile)%filenmData)
  ! check if the forcing file exists
  inquire(file=trim(infile),exist=xist)
  if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   err=10; return
  endif
  ! open forcing data file
  mode=nf90_NoWrite
  call nc_file_open(trim(infile),mode,ncid,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  err = nf90_inq_varid(ncid,'time',varId)
  if(err/=0)then; message=trim(message)//'time data not present'; return; endif
  err = nf90_inquire_attribute(ncid,varId,'units',len = attLen)
  if(err/=0)then; message=trim(message)//'time units attribute not present'; return; endif
  err = nf90_get_att(ncid,varid,'units',refTimeString)
  if(err/=0)then; message=trim(message)//'trouble reading time units'; return; endif
  
  call extractTime(refTimeString,                         & ! input  = units string for time data
                   refTime%var(iLookTIME%iyyy),           & ! output = year
                   refTime%var(iLookTIME%im),             & ! output = month
                   refTime%var(iLookTIME%id),             & ! output = day
                   refTime%var(iLookTIME%ih),             & ! output = hour
                   refTime%var(iLookTIME%imin),dsec,      & ! output = minute/second
                   err,cmessage)                            ! output = error code and error message
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! convert the reference time to days since the beginning of time
  call compjulday(refTime%var(iLookTIME%iyyy),            & ! input  = year
                  refTime%var(iLookTIME%im),              & ! input  = month
                  refTime%var(iLookTIME%id),              & ! input  = day
                  refTime%var(iLookTIME%ih),              & ! input  = hour
                  refTime%var(iLookTIME%imin),dsec,       & ! input  = minute/second
                  refJulday,err,cmessage)                   ! output = julian day (fraction of day) + error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! find first timestep in any of the forcing files
  iFFile=1
  nFile=size(forcFileInfo(:)%data_id(1))
print *,'number forcing files: ',nFile
  do while(.not.foundTime)
   ! how many time steps in current file?
   err = nf90_inq_dimid(ncid,'time',dimId)
   err = nf90_inquire_dimension(ncid,dimId,len=dimLen)
   if(err/=0)then; message=trim(message)//'trouble reading time variable size'; return; endif
   ! read time variable from current file
   if(allocated(fileTime))then; deallocate(fileTime); endif
   if(allocated(diffTime))then; deallocate(diffTime); endif
   allocate(fileTime(dimLen),diffTime(dimLen),stat=err)
   if(err/=0)then; message=trim(message)//'trouble allocating space for time'; return; endif
   err = nf90_get_var(ncid,varId,fileTime,start=(/1/),count=(/dimLen/))
   fileTime=fileTime+refJulday
   !find difference of fileTime from currentJulday
   diffTime=abs(fileTime-currentJulday)
print *,'time: ',currentJulday,refJulDay,minval(diffTime),size(fileTime)
!print *,fileTime
   if(any(diffTime < 1e-9))then
    iRead=minloc(diffTime,1)
    foundTime = .True.
    iFile=iFFile
   else
    ! increment
    iFFile=iFFile+1
    !if we've gone beyone the end of the forcing file list, stop and error
    if(iFFile>nFile)then; err=99; message=trim(message)//'first requested simulation timestep not in any forcing file'; return; endif
    ! create new file name
    ! define filename
    infile=trim(INPUT_PATH)//trim(forcFileInfo(iFFile)%filenmData)
    ! check if the forcing file exists
    inquire(file=trim(infile),exist=xist)
   if(.not.xist)then
    message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
    err=10; return
    endif
    ! open next forcing data file
    mode=nf90_NoWrite
    call nc_file_open(trim(infile),mode,ncid,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   endif
  enddo ! (end of search for model first time step in forcing files)
print *, 'found first timestep: ',iRead
  !convert to time vector
  call compcalday(fileTime(iRead),                   & ! input  = julian day
                 time_data%var(iLookTIME%iyyy),      & ! output = year
                 time_data%var(iLookTIME%im),        & ! output = month
                 time_data%var(iLookTIME%id),        & ! output = day
                 time_data%var(iLookTIME%ih),        & ! output = hour
                 time_data%var(iLookTIME%imin),dsec, & ! output = minute/second
                 err,cmessage)                         ! output = error control
print *,'Opening new forcing file: ',trim(infile)
print *,'first time vector',time_data%var(:)
   
 endif  ! if the file is not yet open

 ! **********************************************************************************************
 ! ***** part 1: if file open, check to see if we've reached the end of the file, if so close it, 
 ! *****         and open new file
 ! *****         Then read the data
 ! **********************************************************************************************
 if(ncid>0)then
  !check to see if we've passed end of netcdf file
  if(iRead>forcFileInfo(iFile)%nTimeSteps)then
   err = nf90_close(ncid)
   if(err/=0)then; message=trim(message)//'problem closing file ['//trim(infile)//']'; return; endif
   ncid = -999
   !increment iFile so we open next forcing file
   iFile = iFile+1
   ! define new filename
   infile=trim(INPUT_PATH)//trim(forcFileInfo(iFile)%filenmData)
   ! check if the forcing file exists
   inquire(file=trim(infile),exist=xist)
   if(.not.xist)then
    message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
    err=10; return
   endif
   ! open forcing data file
   mode=nf90_NoWrite
   call nc_file_open(trim(infile),mode,ncid,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
print *,'Opening new forcing file: ',trim(infile),'step ',iRead
   !reset iRead since we opened a new file
   iRead=1
  endif

  ! **********************************************************************************************
  ! ***** part 1b: read data
  ! **********************************************************************************************

  ! initialize time and forcing data structures
  time_data%var(:) = integerMissing
  forc_data%var(:) = amiss
  !read time data
  err = nf90_inq_varid(ncid,'time',varId)
  !read time at iRead location in netcdf file
  err = nf90_get_var(ncid,varId,dataJulDay,start=(/iRead/))
  if(err/=0)then; message=trim(message)//'trouble getting time value'; return; endif
  dataJulDay = dataJulDay + refJulday
! print *,'reading step',iRead,iHRU,currentJulday,abs(currentJulday - dataJulDay),epsilon(julDayFirst)
  if(abs(currentJulday - dataJulDay) > 1e-9)then
   write(message,'(a,i0,f18.8,a,f18.8,a)') trim(message)//'date for time step: ',iStep,dataJulDay,' differs from the expected date: ',currentJulDay,' in file: '//trim(infile)
   err=40; return
  endif
  !convert julian day to time vector
  call compcalday(dataJulDay,                         & ! input  = julian day
                  time_data%var(iLookTIME%iyyy),      & ! output = year
                  time_data%var(iLookTIME%im),        & ! output = month
                  time_data%var(iLookTIME%id),        & ! output = day
                  time_data%var(iLookTIME%ih),        & ! output = hour
                  time_data%var(iLookTIME%imin),dsec, & ! output = minute/second
                  err,cmessage)                         ! output = error control
  ! check to see if any of the time data is missing
  if(any(time_data%var(:)==integerMissing))then
   do iline=1,size(time_data%var)
    if(time_data%var(iline)==integerMissing)then; err=40; message=trim(message)//"variableMissing[var='"//trim(time_meta(iline)%varname)//"']"; return; endif
   enddo
  endif
  !setup count,start arrays
  ncStart = (/iHru,iRead/)
  !get hruId
  err = nf90_inq_varid(ncid,'hruId',varId)
  if(err/=0)then; message=trim(message)//'hruId not present'; return; endif
  err = nf90_get_var(ncid,varId,hruId,start=(/iHru/))
  if(err/=0)then; message=trim(message)//'Trouble reading current hruId'; return; endif
  !is HRU being read a match to gru_struc?
 !print *,'netcdf open',iRead,iHRU,gru_struc(iGRU)%hru(iHRU)%hru_id,currentJulday
  if(gru_struc(iGRU)%hru(iHRU)%hru_id /= hruId)then
   write(message,'(a,i0,i0,a,i0,a,a)') trim(message)//'hruId for iHRU: ',iHRU,hruId,'differs from the expected:',     &
                                                       gru_struc(iGRU)%hru(iHRU)%hru_id,'in file',trim(infile)
   write(message,'(a)') trim(message)//'order of hruId in forcing file needs to match order in zLocalAttributes.nc'
   err=40; return
  endif  
  !read data into forcing structure
  do iNC=1,forcFileInfo(iFile)%nVars
  !inqure about current variable name
   err = nf90_inquire_variable(ncid,iNC,name=varName)
   if(err/=0)then; message=trim(message)//'problem inquiring variable: '//trim(varName); return; endif
   select case(trim(varName))
    !forcing variables
    case('pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum')
     !get index of forcing variable in forcing data vector
     ivar = get_ixforce(trim(varname))
     !get forcing data
     err=nf90_get_var(ncid,forcFileInfo(iFile)%data_id(ivar),forc_data%var(ivar),start=ncStart)
     !for time, convert days since reference to seconds since reference
    case('time')
     ivar = get_ixforce(trim(varname))
     ! get time from netcdf    
     err=nf90_get_var(ncid,forcFileInfo(iFile)%data_id(ivar),tempVar,start=(/iRead/))
     forc_data%var(ivar) = tempVar*secprday
     !for lat,lon,hruId,data_step do nothing
!    case('hruId','latitude','longitude','data_step')
!    ! check that variables are what we expect
!    case default
!     message=trim(message)//'unknown variable ['//trim(varName)//'] in local attributes file'
!     err=20; return
   end select
  enddo

!print *,'netcdf open',iRead,iHRU,currentJulday
 endif  !end ncid open check


 ! compute the julian day at the start of the year
 call compjulday(time_data%var(iLookTIME%iyyy),          & ! input  = year
                 1, 1, 1, 1, 0._dp,                      & ! input  = month, day, hour, minute, second
                 startJulDay,err,cmessage)                 ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! compute the fractional julian day for the current time step
 call compjulday(time_data%var(iLookTIME%iyyy),           & ! input  = year
                 time_data%var(iLookTIME%im),             & ! input  = month
                 time_data%var(iLookTIME%id),             & ! input  = day
                 time_data%var(iLookTIME%ih),             & ! input  = hour
                 time_data%var(iLookTIME%imin),0._dp,     & ! input  = minute/second
                 currentJulday,err,cmessage)                ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! compute the time since the start of the year (in fractional days)
 fracJulday = currentJulday - startJulDay

 ! compute the number of days in the current year
 yearLength = 365
 if(mod(time_data%var(iLookTIME%iyyy),4) == 0)then
  yearLength = 366
  if(mod(time_data%var(iLookTIME%iyyy),100) == 0)then
   yearLength = 365
   if(mod(time_data%var(iLookTIME%iyyy),400) == 0)then
    yearLength = 366
   endif
  endif
 endif

 ! test
 if(checkTime)then
  write(*,'(i4,1x,4(i2,1x),f9.3,1x,i4)') time_data%var(iLookTIME%iyyy),           & ! year
                                         time_data%var(iLookTIME%im),             & ! month
                                         time_data%var(iLookTIME%id),             & ! day
                                         time_data%var(iLookTIME%ih),             & ! hour
                                         time_data%var(iLookTIME%imin),           & ! minute
                                         fracJulday,                              & ! fractional julian day for the current time step
                                         yearLength                                 ! number of days in the current year
  !pause ' checking time'
 endif

 end subroutine read_force

end module read_force_module