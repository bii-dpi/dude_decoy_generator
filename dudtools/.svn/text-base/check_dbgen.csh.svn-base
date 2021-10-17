#!/bin/tcsh

# might want to report all errors, even if db was made?

foreach h ( * )
  if ( -d $h ) then
  foreach i ( $h/* )
    if ( -d $i ) then
    foreach j ( $i/* )
      if ( -d $j ) then
      foreach k ( $j/* )
        if ( -d $k ) then
        set nonomatch
        foreach l ( $k/TEMP* )
          if ( -d $l ) then
            if ( ! -e $l/db.db.bz2 ) then
              echo $l missing db.db.bz2
              cat $l/sge.log | grep -v "^TEMP" | grep -v "^Run"
            else if ( `ls -la $l/db.db.bz2 | awk '{print $5}'` < 100 ) then
              echo $l db.db.bz2 too small
              cat $l/sge.log | grep -v "^TEMP" | grep -v "^Run"
            else if  ( `grep -ic 'error\|fail' $l/sge.log` > 0 ) then
              echo $l failures
              cat $l/sge.log | grep -v "^TEMP" | grep -v "^Run"
            endif
          endif
        end
        unset nonomatch
        endif
      end
      endif
    end
    endif
  end
  endif
end

