pro calibration_svd,Tsky_off,Tsky_on,noise,flagtot_off,flagtot_on,$
;flagtot_off,flagtot_on,average_freq,avedata_freq,avedata_time,data,
ra,dec,freq,pol,parallatic,fout=fout,thisif=thisif,scan=scan,recn=recn,tave=tave,threshould=threshould,fn=fn,daisy=daisy,period=period,galactic=galactic,drift=drift,sigmacut=sigmacut,mkplot=mkplot
;
; sdfits convention:
; https://wikio.nrao.edu/bin/view/Main/SdfitsDetails
; 
; polarization convention
;
; RR -1  LL -2  RL -3  LR -4
; XX -5  YY -6  XY -7  YX -8
;
; currently, data convention is  -5, -7, -8, -6, noise source low on
; XX ON, XX OFF, XY ON, XY OFF, YX ON, YX OFF, YY ON, YY OFF, etc.  

; recn=[rec0,rec1]:  record number to start-finish
;fn='/home/sdfits/AGBT08A_081_00/AGBT08A_081_00.raw.acs.fits'
;fn='/home/scratch/tchang/AGBT06C_048_08.raw.acs.fits'

; plots
if keyword_set(mkplot) then set_plot,'ps'

Npol=4
;

; read data from the passed file name.  Read only specified records if supplied.
if not keyword_set(recn) then begin
    data=mrdfits(fn,1,header)
    rec0=0
    rec1=n_elements(data.scan)
endif else begin
    rec0=recn(0)
    rec1=recn(1)
    data=mrdfits(fn,1,header,range=[rec0,rec1])
endelse

; parameters to decide on the length of each timedump
;
; timedump = Npol x NIF x Noisecal
;          = 4 x 8 x 2

; ignore the overlapping frequency channels for now, just treat them
; as 2048x8

; figure out the IF
if not keyword_set(thisif) then thisif=0
if not keyword_set(threshould) then threshould=0.02

; The rough frequency of each IF spectrometer window
obsfreq=round(data.CRVAL1/1e6)
restfreq=data.RESTFREQ

; just get a sorted array of spectrometer IF frequencies
dif=obsfreq(uniq(obsfreq, sort(obsfreq)))
IFnum=n_elements(dif) ; Total number of IFs.

; check to make sure a sane IF number was passed
if thisif gt (IFnum-1) then begin
    print,'wrong IF setting'
    print, 'setting IF to',IFnum-1
    thisif=IFnum-1
endif
ifn=strtrim(string(thisif),2)
print,'object name:',data(0).object
print,'number of IF:',IFnum
print,'thisif:',thisif

;ind=where(data.crval1 eq dif(thisif) and data.object eq
;data(rec0).object and data.scan eq scan,cind) 

; The records corresponding to the this IF.  There are 2 x NPol of them.
; without checking scan number
ind=where(obsfreq eq dif(thisif) and data.object eq data(0).object,cind)

print,'observing freq:',dif(thisif)
print,'rest freq:',data(0).RESTFREQ

scannum=cind ; most appropriately described as nrec/nif.
; = n times x npol x 2 (cal on off) x 2


print,'scan number:',scannum
print, 'ending records of first IF:',ind(cind-1)

ntime=scannum
;ntime=rec1-rec0+1
;if (ntime gt scannum) then ntime=scannum
print,'number of scans being process:',ntime
print,'should set rec1-rec0+1 to be:',scannum

; find out tave by setting it to the timestamps of 1 scan
; No idea what this is or why it's an average time .... -km
; It looks more like the number of times....  -km
if not keyword_set(tave) then begin
    ind=where(obsfreq eq dif(thisif) and data.object eq data(0).object and data.scan lt (data(0).scan+6),cind2)
    tave=cind2/(Npol*2)
    print,'tave = ',tave
endif


; now deal with one If at a time, Throw away other IFs.
data=data(ind)


; calculate the "off" signal by averaging over time

Nif=1
nrec=Npol*Nif*2 ; number of records per time bin
nave=(ntime/nrec)/tave
nfreq=n_elements(data(0).data)
freq=((dindgen(nfreq)-(data(0).CRPIX1-1))*data(0).cdelt1/1e6)+double(dif(thisif))
Tsky=dblarr(nave*tave,nrec,nfreq) ; array to store Temperature
Televation=dblarr(nave*tave) ; three arrays for coordinates.
Tazimuth=dblarr(nave*tave)
LST=dblarr(nave*tave)

print,'nrec,ntime,nave,nfreq,ntime/nrec',nrec,ntime,nave,nfreq,ntime/nrec

; check NAN's.  Array for flags.  0 means bad entry.
flagtot=make_array(nave*tave,nrec,nfreq,/double,value=1.)

pol=intarr(nave*tave,nrec) ; some sort of polarization flag that tells you
; what kind of polarization your looking at. -km

; This loop stores data in Tsky then sets any NAN records to 0
for i=0,nave*tave-1 do begin
   ; store all the information related to pointing
   Televation(i)=data(i*nrec).ELEVATIO
   Tazimuth(i)=data(i*nrec).AZIMUTH
   LST(i)=data(i*nrec).LST
   for j=0,nrec-1 do begin
       ; store the sky temperature
       Tsky(i,j,*)=data(i*nrec+j).data
       if keyword_set(whatever) then begin  ; don't do any of this.
       if (i eq 26 and j eq 0 ) then begin
          help,Tsky(i,j,*)
          print,'i,j',i,j
          ind=where(~finite(Tsky),cind)
          print,'NAN:',cind
          if cind lt 2048 then plot,Tsky(i,j,*)
       endif
       if (i eq 28 and j eq 0 ) then begin
          help,Tsky(i,j,*)
          print,'i,j',i,j
          ind=where(~finite(Tsky),cind)
          print,'NAN:',cind
          if cind lt 2048 then plot,Tsky(i,j,*)
       endif
       if (i eq 43 and j eq 1 ) then begin
          help,Tsky(i,j,*)
          print,'i,j',i,j
          ind=where(~finite(Tsky),cind)
          print,'NAN:',cind
          if cind lt 2048 then plot,Tsky(i,j,*)
       endif
       if (i eq 45 and j eq 1 ) then begin
          help,Tsky(i,j,*)
          print,'i,j',i,j
          ind=where(~finite(Tsky),cind)
          print,'NAN:',cind
          if cind lt 2048 then plot,Tsky(i,j,*)
          endif
       if (i eq 43 and j eq 0 ) then begin
          help,Tsky(i,j,*)
          print,'i,j',i,j
          ind=where(~finite(Tsky),cind)
          print,'NAN:',cind
          if cind lt 2048 then plot,Tsky(i,j,*)
       endif 
       if (i eq 45 and j eq 0 ) then begin
          help,Tsky(i,j,*)
          print,'i,j',i,j
          ind=where(~finite(Tsky),cind)
          print,'NAN:',cind
          if cind lt 2048 then plot,Tsky(i,j,*)
       endif
    endif

       ; These lines check for NANs and sets those entries to 0
       ; Also flags them so they don't get used later.
       ind=where(~finite(Tsky(i,j,*)),cind)
       ; add this back
       if (cind gt 0) then begin
           flagtot(i,j,ind)=dblarr(cind)       
           Tsky(i,j,ind)=dblarr(cind)
       endif
       ; Fill an array with the polarization keys (-5,-7,-6,-8) for each time 
       ; bin and each sub record.  
      pol(i,j)=data(i*nrec+j).crval4
   endfor
endfor

;if keyword_set(whatever) then begin



; These two for loops could be easily vectorized.  Don't bother with the if
; statement in the first loop because you zero it in the second loop anyway.

; do hanning smoothing
; basically a low pass filter in frequency space.
; This probably reduces discretization noise. -km
Tsky2=Tsky
for i=0,nave*tave-1 do begin
    for j=0,nrec-1 do begin
        for k=1,nfreq-2 do begin
                                ;if ( finite((Tsky2(i,j,k-1))) and
                                ;finite((Tsky2(i,j,k))) and
                                ;finite((Tsky2(i,j,k+1))) ) then $
            if ( flagtot(i,j,k-1) and flagtot(i,j,k) and flagtot(i,j,k+1)) then $
                Tsky(i,j,k)=Tsky2(i,j,k-1)*0.25d0 + $
                  Tsky2(i,j,k)*0.5d0 + $
                  Tsky2(i,j,k+1)*0.25d
        endfor
    endfor
endfor

; get rid of points that couldn't be smoothed due to adjacent bad
; entries.  This loop could be merged with the last one (with an 'else').
for i=0,nave*tave-1 do begin
    for j=0,nrec-1 do begin
        for k=1,nfreq-2 do begin
            if ( flagtot(i,j,k) eq 0) then begin
               Tsky(i,j,(k-1))=0.
               Tsky(i,j,(k+1))=0.
               ; Should we not also plag these entries? -km
            endif
        endfor
    endfor
endfor

; do hanning smoothing
;data2=data
;for i=0,nave*tave-1 do begin
;    for j=0,nrec-1 do begin
;        for k=1,nfreq-2 do begin
;            data(i*nrec+j).data(k)=data2(i*nrec+j).data(k-1)*0.25d0 + $
;              data2(i*nrec+j).data(k)*0.5d0 + $
;              data2(i*nrec+j).data(k+1)*0.25d0
;        endfor
;    endfor
;endfor


; calculate the noise source and subtract it from the noise-on data
; Tsig_off

;allowcate memory
noise_time=dblarr(nave,nrec/2,nfreq)
noise_bpass=noise_time
noise=dblarr(nave*tave,nrec/2,nfreq) ; stores noise cal temperature
noiseflag=noise ; stores a flag if either the Cal on or cal off are bad
noisetot=noise ; stores the sum of cal on and cal off.  Note sure what for. -km
;
; do median in time
for j=0,nrec-1,2 do begin
   ;subtract cal off from cal on to get the cal temp
   noise(*,j/2,*)=Tsky(*,j,*)-Tsky(*,j+1,*)
   noisetot(*,j/2,*)=Tsky(*,j,*)+Tsky(*,j+1,*)
   ; both can on and cal off need to be 1 or the noise bin is flagged
   noiseflag(*,j/2,*)=flagtot(*,j,*)*flagtot(*,j+1,*)
   for k=0,nfreq-1 do begin     
      medd=(Tsky(*,j,k)-Tsky(*,j+1,k))
      medd2=Tsky(*,j+1,k)
      ind=where(flagtot(*,j,k)*flagtot(*,j+1,k) eq 1,cind)
      ind=where(flagtot(*,j+1,k) eq 1,cind)
      if cind gt 1 then begin
         noise_time(0,j/2,k)=median(medd(ind),/double)
         noise_bpass(0,j/2,k)=median(medd2(ind),/double)
         for i=0,tave*nave-1 do begin
            ;noise(i,j/2,k)=noise(i,j/2,k)/noise_bpass(0,j/2,k)
            ;noise(i,j/2,k)=noise(i,j/2,k)/noise_time(0,j/2,k)-1.
             ;noise(i,j/2,k)=noise(i,j/2,k)-noise_time(0,j/2,k)
         endfor
      endif else noiseflag(*,j/2,k)=0
    endfor
endfor

; take time median of the calibrated noise
noisecal=dblarr(nrec/2,nfreq)
for j=0,nrec/2-1 do begin
   for k=0,nfreq-1 do begin
      ind=where(noiseflag(*,j,k) eq 1,cind)
      if cind gt 1 then begin
         medd=noise(*,j,k)
         noisecal(j,k)=median(medd(ind),/double)
         ;noise(*,j,k)=noise(*,j,k)-noisecal(j,k)
      endif
   endfor
endfor


; subtract the noise source from noise-on data
for i=0,nave*tave-1 do begin
   for j=0,nrec-1,2 do begin
      for k=0,nfreq-1 do begin
         Tsky(i,j,k)=Tsky(i,j,k)-noise_time(i/tave,j/2,k)
      endfor
   endfor
endfor

if keyword_set(svdnoise) then begin
;
noisesvd=dblarr(nave*tave,nrec/2,nfreq)
noisesvd1mode=noisesvd
noisesvd2mode=noisesvd
noisesvd5mode=noisesvd
noisesvd10mode=noisesvd
;   
; do a svd to the noise source
noise=noise*noiseflag
noisetot=noisetot*noiseflag
;
;do a sigmacut
du=50e6/double(nfreq)
dt=0.5
sigma=sqrt(2.)/sqrt(du*dt)*5.
;print,'sigma',min(sigma),max(sigma),sqrt(variance(sigma))
ind=where(abs(noise) gt sigma,cind)
if cind gt 0 then noiseflag(ind)=0.
noise=noise*noiseflag

for j=0,nrec/2-1 do begin
   noise2d=reform(noise(*,j,*))
   la_svd,noise2d, s,u,v,/double,status=status
   print,'SVD status:',status

   w=dblarr(n_elements(s))
   w2=w
   w5=w
   w10=w
   w(0)=s(0)
   w2(0:1)=s(0:1)
   w5(0:4)=s(0:4)
   w10(0:9)=s(0:9)

   noisesvd(*,j,*) = u ## diag_matrix(s) ## transpose(v)
   noisesvd1mode(*,j,*) =u ## diag_matrix(w) ## transpose(v)
   noisesvd2mode(*,j,*) = u ## diag_matrix(w2) ## transpose(v)
   noisesvd5mode(*,j,*) =u ## diag_matrix(w5) ## transpose(v)
   noisesvd10mode(*,j,*) = u ## diag_matrix(w10) ## transpose(v)
endfor

noiseres = noise-noisesvd
noiseres1mode = noise-noisesvd1mode
noiseres2mode = noise-noisesvd2mode
noiseres5mode = noise-noisesvd5mode
noiseres10mode = noise-noisesvd10mode

device,filename='Tnoise.'+strtrim(string(thisif),2)+'.ps'
plot,noise(100,0,*)
plot,noise(150,1,*)
plot,noise(200,2,*)
plot,noise(250,3,*)
plt_image,transpose(reform(noise(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noisesvd(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noiseres(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noisesvd1mode(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noiseres1mode(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noisesvd2mode(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noiseres2mode(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noisesvd5mode(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noiseres5mode(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noisesvd10mode(*,0,*))),/scalable,/colbar
plt_image,transpose(reform(noiseres10mode(*,0,*))),/scalable,/colbar

device,/close

endif

if keyword_set(whatever) then begin

; flagging the RFI's
;  -- use noise_cal off
;flagtot=make_array(nave*tave,nrec,nfreq,/double,value=1.)
flagnumon=0L
flagnumoff=0L
dnu=50e6/double(nfreq)
dt=0.5d
sigma=2./(dnu*dt) ;where did the 2 come from?
for i=0,nave*tave-1 do begin
   ; noise-off
    r= (Tsky(i,3,*)^2+Tsky(i,5,*)^2.) / $
      (Tsky(i,1,*)*Tsky(i,7,*))
;    r=(data(i*nrec+1*2+1).data^2.+data(i*nrec+2*2+1).data^2.) / $
;      (data(i*nrec+0*2+1).data*data(i*nrec+3*2+1).data)
;    plot,r
;    read,test
;    ind=where(r gt threshould,cind)
    ind=where(r gt sigma*(5.^2),cind)
    if (cind gt 0) then begin
       for j=1,nrec-1,2 do begin
          flagtot(i,j,ind)=dblarr(cind)
       endfor
       flagnumoff=flagnumoff+cind
    endif
    ;
    ; noise-on
    r= (Tsky(i,2,*)^2+Tsky(i,4,*)^2.) / $
      (Tsky(i,0,*)*Tsky(i,6,*))
    ind=where(r gt sigma*(5.^2),cind)
    if (cind gt 0) then begin
        for j=0,nrec-1,2 do begin
            flagtot(i,j,ind)=dblarr(cind)
        endfor
        flagnumon=flagnumon+cind
     endif
endfor

print,'flagged points off',flagnumoff
print,'flagged percentage off',double(flagnumoff)/double(n_elements(Tsky(*,0,*)))
print,'flagged points on',flagnumon
print,'flagged percentage on',double(flagnumon)/double(n_elements(Tsky(*,0,*)))

;print, 'total points',nfreq*nave*tave


; Tsig_off
avedata_time=dblarr(nave,nrec,nfreq)
;
; do median in time
for j=0,nrec-1 do begin
    for k=0,nfreq-1 do begin
        for i=0,nave-1 do begin
            medd=Tsky(i*tave:(i+1)*tave-1,j,k)
            flagpart=flagtot(i*tave:(i+1)*tave-1,j,k)
            ind=where(flagpart eq 1,cind)
            if (cind gt 1) then begin
                med=median(medd(ind),/double)
                avedata_time(i,j,k)=med
            endif 
            if (cind eq 1) then avedata_time(i,j,k)=medd(ind)
        endfor
    endfor
endfor



;; what is this?   why is it here?  10/19/2008
;;subtract off the time median from the polarization data
;for j=2,5 do begin
;    for k=0,nfreq-1 do begin
;       Tsky(*,j,k)=Tsky(*,j,k)-avedata_time(0,j,k)
;     endfor
; endfor

if keyword_set(mkplot) then begin

device,filename='Tpol-noise-med'+ifn+'.ps'
plt_image, transpose(reform((Tsky(*,2,*)>(-0.1)<(0.1))*flagtot(*,2,*))),/colbar,/scalable
device,/close

endif



; do T=Tsky/Tref, where Tref=avedata_time
;ind=where(avedata_time eq 0.,cind)
;print,'zero elements of Tref:',cind
;if (cind eq 0) then begin
;   for i=0,nave*tave-1 do begin
;      Tsky(i,*,*)=Tsky(i,*,*)/avedata_time(i/tave,*,*)-1.
;   endfor
;endif else begin
for i=0,nave*tave-1 do begin
   for j=0,nrec-1 do begin
;      avetimeslice=avedata_time(i/tave,j,*)
      if j eq 0 or j eq 1 or j eq 6 or j eq 7 then $
         avetimeslice=avedata_time(i/tave,j,*) else $
            if j eq 2 or j eq 4 then $
               avetimeslice=sqrt(avedata_time(i/tave,0,*)*avedata_time(i/tave,6,*)) else $
                  avetimeslice=sqrt(avedata_time(i/tave,1,*)*avedata_time(i/tave,7,*))
      for k=0,nfreq-1 do begin
         if (avetimeslice(k) eq 0. or ~finite(avetimeslice(k))) then begin
            Tsky(i,j,k)=0.
            flagtot(i,j,k)=0 
         endif else begin
            Tsky(i,j,k)=Tsky(i,j,k)/avetimeslice(k)
         endelse
      endfor
   endfor
endfor
;
for j=0,1 do begin
   Tsky(*,j,*)=Tsky(*,j,*)-1.
endfor
for j=6,7 do begin
   Tsky(*,j,*)=Tsky(*,j,*)-1.
endfor

;endelse

if keyword_set(mkplot) then begin

device,filename='Tsky-noise'+ifn+'.ps'
plt_image, transpose(reform(Tsky(*,0,*)*flagtot(*,0,*))),/colbar,/scalable, xtitle='!6time',ytitle='freq'
device,/close

device,filename='Tsky_nonoise'+ifn+'.ps'
plt_image, transpose(reform(Tsky(*,1,*)*flagtot(*,1,*))),/colbar,/scalable, xtitle='!6time',ytitle='freq'
device,/close

endif


if keyword_set(sigmacut) then begin

; do median in central frequency
flagout=0L
Tshould=0.03
fstart=floor(nfreq/5.)
fend=ceil(nfreq/5.*4.)
avedata_freq=dblarr(nave*tave,nrec)
for i=0,nave*tave-1 do begin
   for j=0,nrec-1 do begin
      flagpart=flagtot(i,j,fstart:fend-1)
      Tskypart=Tsky(i,j,fstart:fend-1)
      ind=where(flagpart eq 1 and abs((Tskypart)-1.) lt (Tshould),cind)
      if (cind gt 1) then $
         med=median(Tskypart(ind),/double) $
      else if (cind eq 1) then $
         med=Tskypart(ind) 
      if (cind eq 0) then begin
         med=1.
         flagtot(i,j,*)=dblarr(nfreq)
         flagout=flagout+1L
      endif
      Tsky(i,j,*)=Tsky(i,j,*)/med-1.
      avedata_freq(i,j)=med
   endfor
endfor
print,'all freq flagged:',flagout

endif

endif



;
;
;

; calculate ra,dec
ra=dblarr(nave*tave)
dec=ra
if (keyword_set(daisy) or keyword_set(drift)) then begin
   GBTLAT=38.433119
   GBTLON=-79.839833
   for i=0,nave*tave-1 do begin
      t1=data(i*nrec).date_obs
      year=fix(strmid(t1,0,4))
      month=fix(strmid(t1,5,2))
      day=fix(strmid(t1,8,2))
      hour=float(strmid(t1,11,2))
      minute=float(strmid(t1,14,2))
      sec=float(strmid(t1,17,5))
      time=[year,month,day,hour,minute,sec]
      az=data(i*nrec).CRVAL2
      el=data(i*nrec).CRVAL3
      juldate,time,jd
      HOR2EQ, el, az, (jd+2400000.d), ra1, dec1, LAT=GBTLAT, LON=GBTLON, $
              REFRACT_= 0
      ra(i)=ra1
      dec(i)=dec1
   endfor
endif else begin
    for i=0,nave*tave-1 do begin
        ra(i)=data(i*nrec).crval2
        dec(i)=data(i*nrec).crval3
    endfor
endelse
if keyword_set(galactic) then begin
    gl=dblarr(nave*tave)
    gb=gl
    for i=0,nave*tave-1 do begin
        gl(i)=data(i*nrec).crval2
        gb(i)=data(i*nrec).crval3
    endfor
    euler,ra,dec,gl,gb,2
endif 



; calculate the parallatic angle -- according to Myers 
;                                   http://www.aoc.nrao.edu/~smyers/Synth2004/MyersPolarization04print.pdf
; Paralax of the atmosphere changes polarization angle, so this is important.
GBTLAT=38.433119/!radeg
hourangle=(LST/3600.d0*15.-ra)/!radeg
parallatic=atan(cos(GBTLAT)*sin(hourangle), sin(GBTLAT)*cos(dec/!radeg)-cos(GBTLAT)*sin(dec/!radeg)*cos(hourangle))


; subtract off periodic struture in time for Daisy scan (corresponding
; to same position on the sky)
;if keyword_set(daisy) then begin
;   if not keyword_set(period) then period=20
;   index=findgen(nave*tave)
;   for i=0,period-1 do begin
;      ind=where(index mod period eq i,cind)
;      for j=0,nrec-1 do begin
;         for k=0,nfreq-1 do begin
;            tmp=Tsky(ind,j,k)
;            ftmp=flagtot(ind,j,k)
;            ind2=where(ftmp eq 1,cind)
 ;           tmpt=0.
;            if cind gt 0 then $
;               tmpt=mean(tmp(ind2))
;            Tsky(ind,j,k)=Tsky(ind,j,k)-tmpt
;         endfor
;      endfor
;   endfor
   ;
   ; subtract off the residual time fluctuation;
;   Tsky2=Tsky
;   for i=0,nave*tave-1 do begin
;      ind=indgen(period)-period/2+i
;      ind2=where(ind ge 0 and ind lt (nave*tave),cind)
;      for j=0,nrec-1 do begin
;         for k=0,nfreq-1 do begin
;            tmp=Tsky2(ind(ind2),j,k) ;/double(cind)
;            ftmp=flagtot(ind(ind2),j,k)
;            ind3=where(ftmp eq 1,cind3)
;            tmpt=0.
;            if cind3 gt 0 then $
;               tmpt=mean(tmp(ind3))
;            Tsky(i,j,k)=Tsky(i,j,k)-tmpt
;         endfor
;      endfor
;   endfor
;endif

if keyword_set(whatever) then begin

; for drfit scans
; taking the median of time to be the bandpass and divide out the
; bandpass variation, then we want to subtract off the receiver noise,
; which can be estimated by taking the average of the scan (in time),
; but need to flag out the RFI's first by flagging out 5-sigma
; points, and in order to calculate the sigma (and compared it to
; theoretical expectations) one needs to subtract off the foregrounds
; first - by taking the median in some frequency range
;
; subtract off the foreground
; calculate the median in some freq range and subtract it off
; not counting the flagged points
dnfreq=8
dfreq=nfreq/dnfreq
Tshould=0.02
average_freq=dblarr(nave*tave,nrec,dnfreq)
;average_freq=dblarr(nave*tave,nrec,nfreq)
for i=0,nave*tave-1 do begin
   for j=0,nrec-1 do begin
      for k=0,dnfreq-1 do begin
         medd=Tsky(i,j,k*dfreq:(k+1)*dfreq-1)
         flagpart=flagtot(i,j,k*dfreq:(k+1)*dfreq-1)
 ;        ind=where(flagpart eq 1,cind)
         ind=where(flagpart eq 1 and abs((medd)) lt (Tshould),cind)
         if (cind gt 1) then begin
            med=median(medd(ind),/double)
            average_freq(i,j,k)=med
;            for ll=0,dfreq-1 do begin
;;               average_freq(i,j,k*dfreq+ll)=med
;               Tsky(i,j,k*dfreq+ll)=Tsky(i,j,k*dfreq+ll)-med
;            endfor
         endif
      endfor
   endfor
endfor


if keyword_set(whatever) then begin

if keyword_set(mkplot) then begin

device,filename='average_freq'+ifn+'.ps'
plt_image, (reform(average_freq(0:128,1,*))),/colbar,/scalable
device,/close

device,filename='Tsky-fg'+ifn+'.ps'
plt_image, transpose(reform((Tsky(*,1,*))*flagtot(*,1,*))),/colbar,/scalable
device,/close

device,filename='Tpol-fg'+ifn+'.ps'
plt_image, transpose(reform((Tsky(*,3,*))*flagtot(*,3,*))),/colbar,/scalable
device,/close

endif


; sigmacut
; sigma=1./sqrt(dnu dt)
; dt=0.5 sec, dnu=50MHz/nfreq
;
;; ===== no sigma-cut is applied now =========
; == applying sigma-cut for drift scans ========
;if keyword_set(sigmacut) then begin

dnu=50e6/double(nfreq)
dt=0.5d
sigma=1./sqrt(dnu*dt)
sigma_th=5.
ind=where(abs(Tsky) gt (sigma_th*sigma),cind)
if (cind gt 0) then flagtot(ind)=dblarr(cind)
;ind_unpol=[0,1,6,7]
;indt=where(finite(flagtot(*,ind_unpol,*)))
;ind=where(abs(Tsky(indt) gt (3.*sigma),cind)
;if (cind gt 0) then flagtot(indt(ind))=dblarr(cind)
ind_unpol=where(abs(Tsky(*,1,*)) gt (sigma_th*sigma),cind_ll)
ind_unpol=where(abs(Tsky(*,7,*)) gt (sigma_th*sigma),cind_rr)
ind_unpol=where(abs(Tsky(*,3,*)) gt (sigma_th*sigma),cind_lr)
ind_unpol=where(abs(Tsky(*,5,*)) gt (sigma_th*sigma),cind_rl)
ind_unpol=where(abs(Tsky(*,0,*)) gt (sigma_th*sigma),cindon_ll)
ind_unpol=where(abs(Tsky(*,6,*)) gt (sigma_th*sigma),cindon_rr)
ind_unpol=where(abs(Tsky(*,2,*)) gt (sigma_th*sigma),cindon_lr)
ind_unpol=where(abs(Tsky(*,4,*)) gt (sigma_th*sigma),cindon_rl)
;if (cind gt 0) then flagtot(ind)=dblarr(cind)
print,'** sigmacut **', cind, double(cind)/n_elements(Tsky)
print,'sigmacut ll', cind_ll, double(cind_ll)/n_elements(Tsky(*,0,*))
print,'sigmacut rr', cind_rr, double(cind_rr)/n_elements(Tsky(*,0,*))
print,'sigmacut lr', cind_lr, double(cind_lr)/n_elements(Tsky(*,0,*))
print,'sigmacut rl', cind_rl, double(cind_rl)/n_elements(Tsky(*,0,*))
print,'sigmacut ll noise-on', cindon_ll, double(cindon_ll)/n_elements(Tsky(*,0,*))
print,'sigmacut rr noise-on', cindon_rr, double(cindon_rr)/n_elements(Tsky(*,0,*))
print,'sigmacut lr noise-on', cindon_lr, double(cindon_lr)/n_elements(Tsky(*,0,*))
print,'sigmacut rl noise-on', cindon_rl, double(cindon_rl)/n_elements(Tsky(*,0,*))

;;endif

if keyword_set(mkplot) then begin

device,filename='Tsky_sigmacut'+ifn+'.ps'
plt_image, transpose(reform(Tsky(*,1,*)*flagtot(*,1,*))),/colbar,/scalable
device,/close

device,filename='Tpol_sigmacut'+ifn+'.ps'
plt_image, transpose(reform(Tsky(*,3,*)*flagtot(*,3,*))),/colbar,/scalable
device,/close

endif

; add back the mean foreground
for i=0,nave*tave-1 do begin
   for j=0,nrec-1 do begin
      for k=0,dnfreq-1 do begin
;         medd=Tsky(i,j,k*dfreq:(k+1)*dfreq-1)
;         flagpart=flagtot(i,j,k*dfreq:(k+1)*dfreq-1)
;         ind=where(flagpart eq 1,cind)
;         if (cind gt 1) then begin
            for ll=0,dfreq-1 do begin
               Tsky(i,j,k*dfreq+ll)=Tsky(i,j,k*dfreq+ll)+average_freq(i,j,k)
           endfor
;       endif
   endfor
endfor
endfor



; now calculate the average in time
average_time=dblarr(nrec,nfreq)
for j=0,nrec-1 do begin
    for k=0,nfreq-1 do begin
       medd=Tsky(*,j,k)
       flagpart=flagtot(*,j,k)
       ind=where(flagpart eq 1,cind)
       if (cind gt 0) then begin
          med=mean(medd(ind),/double)
          average_time(j,k)=med
          Tsky(*,j,k)=Tsky(*,j,k)-med
       endif
    endfor
endfor

if keyword_set(whatever) then begin
; now calculate the average in time                                                         
average_noise=dblarr(nrec/2,nfreq)
for j=0,nrec/2-1 do begin
    for k=0,nfreq-1 do begin
       average_noise(j,k)=mean(noise(*,j,k),/double)
       noise(*,j,k)=noise(*,j,k)-average_noise(j,k)
    endfor
endfor

endif

endif


if keyword_set(mkplot) then begin

device,filename='noise_final'+ifn+'.ps'
plt_image, transpose(reform(noise(*,0,*))),/colbar,/scalable, xtitle='!6time'
device,/close

device,filename='Tsky-noise_final'+ifn+'.ps'
plt_image, transpose(reform(Tsky(*,0,*)*flagtot(*,0,*))),/colbar,/scalable, xtitle='!6time',ytitle='freq'
device,/close

device,filename='Tsky_nonoise_final'+ifn+'.ps'
plt_image, transpose(reform(Tsky(*,1,*)*flagtot(*,1,*))),/colbar,/scalable, xtitle='!6time',ytitle='freq'
device,/close

endif

;read,test

; apply flagging and rescale to temperature [K]
;Tsky=Tsky*flagtot*40.d0

if keyword_set(sigmacut) then begin
; restore the orignial value
for i=0,nave*tave-1 do begin
   for j=0,nrec-1 do begin
      Tsky(i,j,*)=(Tsky(i,j,*)+1.d)*avedata_freq(i,j)
   endfor
endfor
endif

;for i=0,nave*tave-1 do begin
;   for j=0,nrec-1 do begin
;      for k=0,nfreq-1 do begin
;         Tsky(i,*,*)=Tsky(i,j,k)*avedata_time(i/tave,j,k)
;      endfor
;   endfor
;endfor
;
; apply flagging
;Tsky=Tsky<(0.05)>(-0.05)
;Tsky=Tsky<100.
Tsky=(Tsky*flagtot*40.d0)
noise=noise<(1.)>(-1.)
;Tsky=average_freq

endif



; seperate noise on and noise off
Tsky=(Tsky*flagtot)
noise=noise<(1.)>(-1.)
ind=indgen(nrec)
ind_on=where(ind mod 2 eq 0)
ind_off=where(ind mod 2 eq 1)
Tsky_off=Tsky(*,ind(ind_off),*)
Tsky_on=Tsky(*,ind(ind_on),*)
flagtot_off=flagtot(*,ind(ind_off),*)
flagtot_on=flagtot(*,ind(ind_on),*)


fdir='/home/scratch/tchang/'
if keyword_set(fout) then begin
   openw,1,fdir+fout+'.t'
   printf,1,Tsky
   close,1
   ;
   openw,1,fdir+fout+'.ra'
   printf,1,ra
   close,1
   ;
   openw,1,fdir+fout+'.dec'
   printf,1,dec
   close,1
   ;
   openw,1,fdir+fout+'.freq'
   printf,1,freq
   close,1
   ;
   openw,1,fdir+fout+'.pol'
   printf,1,pol
   close,1
   ;
   openw,1,fdir+fout+'.params'
   printf,1,'nt, nrec, nfreq'
   printf,1, (nave*tave), (nrec/2), nfreq
   close,1
   ;
endif

end


