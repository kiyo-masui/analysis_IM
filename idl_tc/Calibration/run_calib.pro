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
fields =   [  -1,    2]
nsessions = n_elements(sessions)

; steps is an array of 3 flags.  Setting a flag to 0 will skip teh
; corresponding analysis step [docalib_svdfill, cal_Tcal, cal_tcal_matrix]
steps = [1,1,1,1]

; now only worry about the ones I want to process NOW.
whichsessions = [ 1 ]
; whichsessions = indgen(nsessions)

sessions = sessions(whichsessions)
calfiles = calfiles(whichsessions)
nsessions = n_elements(sessions)

for ii=0, nsessions-1 do begin
  session = sessions[ii]
  nfiles = calfiles[ii]
  nfields = fields[ii]
  for jj=0,nfiles-1 do begin
    ; now get the field name and the scan number from the 'library'
    calibration_scans,calfield,scans, session,jj
    ; now we start processing the data one file at a time.
    if steps[0] then begin
      docalib_svdfill,field=calfield,ffdir=scans,dname=session
    endif
  endfor

  for jj=0,nfields-1 do begin
    ; Get the filed and which calibration scans to process from the library
    field_calibrations,field,field_lastscans,calfield,calfile_indicies, session,jj
    n_cals = n_elements(calfield_indicies)
    ffdir = []
    for kk=1, n_cals-1 do begin
      calibration_scans,calfield,scans, session,calfield_indicies[jj]
      ffdir = [ffdir,scans]
    endfor
    if steps[1] then begin
      cal_Tcal,calfield=calfield,ffdir=ffdir,dname=session
    endif
    if steps[2] then begin
      cal_tcal_matrix, field=calfield,ffdir=ffdir,dname=session
    end  
endfor

end
