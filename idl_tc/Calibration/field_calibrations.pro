pro field_calibrations,field,field_lastscans,calfield,calfile_indicies,session,field_num

; This is another directory of which fields were scaned each night and what 
; the corresponding calibration scans are.

if session eq 01 then begin
  if field_num eq 0 then begin
    field = 'wigglez22hr'
    field_lastscans = '161-168'
    calfield = '3c348'
    calfile_indicies = [0]
  endif
  if field_num eq 2 then begin
    field = 'wigglez1hr'
    field_lastscans = '201-208'
    calfield = '3c48'
    calfile_indicies = [1]
  endif
endif

    
