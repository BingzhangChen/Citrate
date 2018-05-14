subroutine  Readcsv(filename,nrow,ncol,dat)
implicit none
character(LEN=20),intent(in)      :: filename
integer,          intent(out)     :: nrow, ncol
character(LEN=10),dimension(ncol) :: header
integer,          parameter       :: stdout=6, funit = 9
real, dimension(nrow,ncol), intent(out):: dat
integer    :: err, ix, IOS
logical    :: I_opened
character(LEN=100000) :: cff

   ! Inquire whether the unit has been opened or not
   INQUIRE (funit, OPENED=I_opened) 

   if (I_opened) then
      print *, 'The file unit already open!'
      close(funit)
   endif

   OPEN(unit=funit,file=filename,status='OLD',iostat=err)

   IF (err /= 0) THEN
     write(stdout,*) 'open ', TRIM(filename),' fails'
     stop
     close(funit)

   ELSE
     read(funit,*,iostat=err) header  !read the first row

     if (err .gt. 0) then 
       write(stdout,*) 'The error is: ', err, 'at Row: ',ix, ' for ',TRIM(filename)
       CLOSE(funit)
       stop
     elseif (err .lt. 0) then
       write(stdout,*) 'The error is: ', err, 'at Row: ',ix, ' for ',TRIM(filename)
       write(stdout,*) 'End of file reached!'
       CLOSE(funit)
       stop
     else
       ! Determine nrow and ncol in the file

       do while (IOS == 0)
          read(fh, '(A)', iostat=ios) buffer
          if (ios == 0) then
             line = line + 1
             
             ! Find the first instance of whitespace.
             ! Split label and data.
             pos = scan(buffer, '    ')
             label = buffer(1:pos)
             buffer =
                                     buffer(pos+1:)do ix = 1, nrow  
         read(funit,*,iostat=IO) cff
         if (err .gt. 0) then 
           write(stdout,*) 'The error is: ', err, 'at Row: ',ix, ' for ',TRIM(filename)
           CLOSE(funit)
           stop
         elseif (err .lt. 0) then
           write(stdout,*) 'The error is: ', err, 'at Row: ',ix, ' for ',TRIM(filename)  
           write(stdout,*) 'End of file reached!'
           CLOSE(funit)
           stop
         else
           dat(ix,:) = cff(:)
         endif
       enddo  
     endif
     dat(:,:) = 0.0
   ENDIF
   CLOSE(funit)
! End reading data
End subroutine Readcsv

