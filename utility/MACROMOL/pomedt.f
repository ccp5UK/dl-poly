      program pomedt

c***********************************************************************
c
c     dl_poly utility to reprocesses files produced by the dl_poly
c     utility pdb2edt, for files with nucleotide residues.
c     In particular residues entries 'POM' (phosphate groups) are
c     seperated from the rest of the residue.
c
c     author t. forester feb 1996
c     copyright daresbury laboratory 1996
c
c***********************************************************************

      parameter(mxsit = 1000)
      character*80 oneline,line(mxsit)
      character*40 input,output

      write(*,*) 'name of file to be processed ?'
      read(*,'(a40)') input
      open(10,file=input)

      write(*,'(/,a)') 'name of output file ?'
      read(*,'(a40)') output
      open(11,file=output)

      do imega=1,2**30

        read(10,'(a80)',end=100,err=100) oneline
        if(oneline(6:12).ne.'RESIDUE') then

          write(11,'(a80)') oneline

        else

          do i = 1,mxsit

            read(10,'(a80)',end=50,err=50) line(i)
            if(line(i)(6:12).eq.'RESIDUE') goto 50

          enddo
          write(*,*) ' unsafe exit from loop '
 50       nlines = i-1
          backspace(10)
c
c     search for phosphate residue entries

          lfnd = 0
          ip = 0
          do i = 1,nlines
            if(line(i)(14:14).eq.'P') then
              lfnd = lfnd+1
              ip = i
            endif
          enddo
c
c     if no phosphate : then write it all out ...

          if(lfnd.eq.0.or.nlines.le.8) then
            write(11,'(a80)') oneline
            do i = 1,nlines
              write(11,'(a80)') line(i)
            enddo
          elseif(lfnd.eq.1) then

            if(ip.le.nlines/2) then

              write(11,'(5x,a,//,5x,a,/)')
     x             'RESIDUE     =  POM','BOND ARRAY BEGINS'

              do k = ip,ip+2

                if((line(k)(14:14).eq.'P').or.
     x             (line(k)(14:14).eq.'O')) then

                  write(11,'(a80)') line(k)

                endif

              enddo

              write(11,*)

              write(11,'(a80)') oneline
              do i = 1,nlines

                if((i.lt.ip).or.(i.gt.ip+2))
     x               write(11,'(a80)') line(i)

              enddo

            else

              write(11,'(a80)') oneline
              do i = 1,nlines

                if((i.lt.ip).or.(i.gt.ip+2))
     x               write(11,'(a80)') line(i)

              enddo

              write(11,'(5x,a,//,5x,a,/)')
     x             'RESIDUE     =  POM','BOND ARRAY BEGINS'

              do k = ip,ip+2

                if((line(k)(14:14).eq.'P').or.
     x             (line(k)(14:14).eq.'O')) then

                  write(11,'(a80)') line(k)

                endif

              enddo

              write(11,*)

            endif

          else

            write(*,*) 'error - more than one POM group in residue!!'

            write(11,'(a80)') oneline
            do i = 1,nlines

              if((i.lt.ip).or.(i.gt.ip+2))
     x             write(11,'(a80)') line(i)

            enddo

              write(11,'(5x,a,//,5x,a,/)')
     x             'RESIDUE     =  POM','BOND ARRAY BEGINS'

            do k = ip,ip+2

              if((line(k)(14:14).eq.'P').or.
     x             (line(k)(14:14).eq.'O')) then

                write(11,'(a80)') line(k)

              endif

            enddo

            write(11,*)

          endif

        endif

      enddo
 100  close(11)
      close(10)
      end
