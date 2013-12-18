# 1 "io_open.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "io_open.F"
c*************************************************************************
c                            IO_OPEN.F
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c open files
c
c             Input:
c                 iu              ==>  unit number (integer scalar)
c                 fname           ==>  file name (character*80)
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c                 format          ==>  format string (character*80)
c             Output:
c                 ierr            ==>  output from iostat
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    3/3/94
c Last revision: 1/30/98

      subroutine io_open(iu,fname,fopenstat,format,ierr)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu
      character*(*) fname,fopenstat,format

c...  Outputs: 
      integer ierr

c----
c...  Executable code 

      if( (fopenstat(1:6).eq.'append') .or. 
     &     (fopenstat(1:6).eq.'APPEND') ) then
         open(unit=iu, file=fname, status='old',

     &        position='append',



     &        form=format,iostat=ierr)
         if(ierr.ne.0) then
            write(*,*) 'Warning:  Could not open ',fname,' with'

            write(*,*) '          position=append.'



            write(*,*) '          Will open as status=new'
            open(unit=iu, file=fname, status='new',
     &           form=format,iostat=ierr)
         endif
      else
         open(unit=iu, file=fname, status=fopenstat,
     &        form=format,iostat=ierr)

      endif


      return
      end   ! io_open
