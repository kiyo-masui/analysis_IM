pro calibration_scans, field, scans, session, filenum

; This is a library of which calibrator fields are pointed to and what the scan
; number are for each calibration.  'field' and 'scans' are outputs and are
; strings.  'session' and 'filenum' are inputs.  'seesion' is a string (like '00')
; and 'filenum' is an integer asking for which calibraiontion file in that
; session. 0 labels the first file.

if session eq '00' then begin
  if filenum eq 0 then begin
    field = '3c348'
    scans = ['5-8']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['161-164']
  endif
  if filenum eq 2 then begin
    field = '3c48'
    scans = ['165-168']
  endif
endif
 
if session eq '01' then begin
  if filenum eq 0 then begin
    field = '3c348'
    scans = ['5-8']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['209-212']
  endif
endif

if session eq '02' then begin
  if filenum eq 0 then begin
    field = '3c348'
    scans = ['5-8']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['161-164']
  endif
endif

if session eq '03' then begin
  if filenum eq 0 then begin
    field = '3c286'
    scans = ['13-16']
  endif
  if filenum eq 1 then begin
    field = '3c286'
    scans = ['184-187']
  endif
endif

if session eq '04' then begin
  if filenum eq 0 then begin
    field = '3c348'
    scans = ['5-8']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['297-300']
  endif
  if filenum eq 2 then begin
    field = '3c48'
    scans = ['304-307']
  endif
  if filenum eq 3 then begin
    field = '3c48'
    scans = ['356-359']
  endif
endif

if session eq '05' then begin
  if filenum eq 0 then begin
    field = '3c286'
    scans = ['5-8']
  endif
  if filenum eq 1 then begin
    field = '3c286'
    scans = ['161-164']
  endif
endif



if session eq '06' then begin
  if filenum eq 0 then begin
    field = '3c348'
    scans = ['112-115']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['222-225']
  endif
endif

if session eq '08' then begin
  if filenum eq 0 then begin
    field = '3c48'
    scans = ['127-130']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['59-62']
  endif
endif

if session eq '10' then begin
  if filenum eq 0 then begin
    field = '3c286'
    scans = ['53-56']
  endif
  if filenum eq 1 then begin
    field = '3c286'
    scans = ['137-140']
  endif
endif

if session eq '11' then begin
  if filenum eq 0 then begin
    field = '3c348'
    scans = ['7-10']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['163-166']
  endif
endif

if session eq '12' then begin
  if filenum eq 0 then begin
    field = '3c48'
    scans = ['5-8']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['133-136']
  endif
  if filenum eq 2 then begin
    field = '3c48'
    scans = ['193-196']
  endif
endif

if session eq '13' then begin
  if filenum eq 0 then begin
    field = '3c48'
    scans = ['11-14']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['119-122']
  endif
endif

if session eq '14' then begin
  if filenum eq 0 then begin
    field = '3c348'
    scans = ['9-12']
  endif
endif

if session eq '15' then begin
  if filenum eq 0 then begin
    field = '3c348'
    scans = ['18-21']
  endif
  if filenum eq 1 then begin
    field = '3c48'
    scans = ['174-177']
  endif
  if filenum eq 2 then begin
    field = '3c48'
    scans = ['250-253']
  endif
endif

end
