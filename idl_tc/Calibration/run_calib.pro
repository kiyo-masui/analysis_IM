pro run_calib

; Driver script for docalib_svdfil.pro (the calibration version).
; This script runs calibration_scans.pro to find out about the
; scan numbers of and calibrators of each session
; and runs the appropriate docalib_svdfill.pro, such that it acctually
; finds the data file.


; These used for looping.  Keep them up to date and do not change this to
; (for example) shorten the loop.  Instead change it below.
sessions = ['00', '01', '02', '03', '04', '05', '06', '08', '10', '11', '12', $
               '13', '14', '15' ]
calfiles = [   3,    2,    2,    2,    4,    2,    2,    2,    2,    2,    3, $
                  2,    1,    3 ]
nsessions = n_elements(sessions)

; steps is an array of 3 flags.  Setting a flag to 0 will skip teh
; corresponding analysis step [docalib_svdfill, cal_Tcal, cal_tcal_matrix]
steps = [0,0,0,1]

; now only worry about the ones I want to process NOW.
whichsessions = [ 0 ]
; whichsessions = indgen(nsessions)

sessions = sessions(whichsessions)
calfiles = calfiles(whichsessions)
nsessions = n_elements(sessions)

for ii=0, nsessions-1 do begin
  session = sessions[ii]
  nfiles = calfiles[ii]
  for jj=0,nfiles-1 do begin
    ; now get the field name and the scan number from the 'library'
    calibration_scans,field,scans,session,jj
    ; now we start processing the data one file at a time.
    if steps[0] then begin
      docalib_svdfill,field=field,ffdir=scans,dname=session
    endif
  endfor

  for jj=0,nfields-1 do begin
    if steps[1] then begin
      cal_Tcal,field=field,ffdir=scans,dname=session
    endif
    if steps[2] then begin
      cal_tcal_matrix, field=field,ffdir=scans,dname=session
    end  
endfor

end
