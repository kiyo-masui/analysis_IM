pro cal_Tcal,dname=dname,field=field,calfield=calfield,ffdir=ffdir

; read in Tcal/Tsys data files

if not keyword_set(field) then field='3c218'  ; for the radio/optical cube setup
if not keyword_set(dname) then dname='04'

outdirroot = '/mnt/raid-project/gmrt/kiyo/wiggleZ/calibration/'

if not keyword_set(ffdir) then begin
  if dname eq '00' then begin
     ffdir=['5-8'];,'34-37']
  endif
  ;
  if dname eq '02' then begin
     ffdir=['5-8']
  endif
  ;
  if dname eq '03' then begin
     ffdir=['5-8']
  endif
  ;
  if dname eq '04' then begin
     ffdir=['5-8']
  endif
endif
;

nscan=4
npol=4
nwin=8
nfreq=2048
nfile=n_elements(ffdir)
t3ctcal=dblarr(nfile,nwin,npol/2,nfreq,nscan/2)
freqtot=dblarr(nwin,nfreq)
tsystcal=dblarr(nfile,nwin,npol/2,nfreq,nscan/2)

for ifile=0,nfile-1 do begin
   

   fdir = outdirroot+dname+'/'+calfield+'_'+ffdir(ifile)+'/'
   ;fdir='/cita/d/raid-project/gmrt/tchang/GBT/zCOSMOS/Calibration/'+dname+'/3c218_'+ffdir(ifile)+'/'
   print, 'Running cal_Tcal on directory: ' + fdir

   for iwin=0,nwin-1 do begin

      ; take 1/scan0-1/scan1=T_3c48/Tcal,
      ; take scan average of (0,1) and (2,3)
      tcaltsys=dblarr(nscan,npol,nfreq)
      freq=dblarr(nfreq)
      openr,1,fdir+'Freq.scan0'+strtrim(string(iwin),2)+'.txt'
      readf,1,freq
      close,1
      freqtot(iwin,*)=freq
      tmp=dblarr(npol,nfreq)
      ;
      for scan=0,nscan-1 do begin
                                ; scans 0,2 are position on = Tcal/(T_3c48+Tsys),
                                ; scans 1,3 are position off = Tcal/Tsys
                                ; ipol=0,3 are LL, RR
                                ; ipol=1,2 are LR, RL
         openr,1,fdir+'TcalTsys.scan'+strtrim(string(scan),2)+strtrim(string(iwin),2)+'.txt'
         readf,1,tmp
         close,1
         tcaltsys(scan,*,*)=tmp
         ;
      endfor
      ;
      ; now calculate the ratio (Tsys+T_3c)/Tcal - Tsys/Tcal
      for iscan=0,nscan/2-1 do begin
         for ipol=0,3,3 do begin
            ;print,tcaltsys(iscan*2,ipol,0:100)
            ;read,test
            ;print,tcaltsys(iscan*2+1,ipol,0:100)
            ;read,test
            ;print,1./tcaltsys(iscan*2,ipol,0:100)-1./tcaltsys(iscan*2+1,ipol,0:100)
            ;read,test
            t3ctcal(ifile,iwin,ipol<1,*,iscan)=1./tcaltsys(iscan*2,ipol,*)-1./tcaltsys(iscan*2+1,ipol,*)
            tsystcal(ifile,iwin,ipol<1,*,iscan)=1./tcaltsys(iscan*2+1,ipol,*)
         endfor
      endfor
      
   endfor

endfor

fdir = outdirroot+dname+'/'+field+'/'
spawn 'mkdir -p '+fdir

t3ctcal=t3ctcal<100.>(-100.)

set_plot,'ps'
device,filename=fdir+'tsystcal.ps'

f=reform(freqtot(0,*))
t=reform(tsystcal(0,0,0,*,0))
plot,f,t,xtitle='!6freq [MHz]',ytitle='Tsys/Tcal [K]'
f=reform(freqtot(0,*))
t=reform(tsystcal(0,0,1,*,0))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
;
f=reform(freqtot(4,*))
t=reform(tsystcal(0,4,0,*,0))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(4,*))
t=reform(tsystcal(0,4,1,*,0))
oplot,f,t,linestyle=2

device,/close

; This data from Baars, Genzel, Pauliny-Toth and Witzel, 1997
; Might want to eventually test the consitancy of these, as temporial
; variation on time scales of decades is expected. -km
if (field eq '3c48') then begin
  const = 2.345
  alpha = 0.071
  delta = -.138
endif
if (field eq '3c348') then begin
  const = 4.963
  alpha = -1.052
  delta = 0
endif
if (field eq '3c218') then begin
  const = 4.497
  alpha = -0.910
  delta = 0
endif
if (field eq '3c286') then begin ; This one is polarized
  const = 1.480
  alpha = 0.292
  delta = -0.124
endif
;
; calculate tcal
;f3c=(freqtot/1000.)^(-0.85584) * 21.37  ;Jy  3c48
;f3c=(freqtot/1400.)^(-0.79949634) * 40.8
;const       3.9912530
;alpha     -0.72159880
;const= 3.9912530
;alpha= -0.72159880
; for 3c218??  May 14, 2010


f3c=10^(const+alpha*alog10(freqtot)+delta*alog10(freqtot)*alog10(freqtot))
t3c=f3c*2  ;gain=2 K/Jy.  Is this not frequency dependant? -km.


tcal=1./t3ctcal
for ifile=0,nfile-1 do begin
   for iwin=0,nwin-1 do begin
      for ipol=0,1 do begin
         for iscan=0,1 do begin
            tcal(ifile,iwin,ipol,*,iscan)=tcal(ifile,iwin,ipol,*,iscan)*t3c(iwin,*)
         endfor
      endfor
   endfor
endfor
;
;average nscan and nfile
tcal2=tcal
tcal=dblarr(2,nwin,nfreq)
tsys=dblarr(2,nwin,nfreq)
for ipol=0,1 do begin
   for iwin=0,nwin-1 do begin
      for ifreq=0,nfreq-1 do begin
         tcal(ipol,iwin,ifreq)=mean(tcal2(*,iwin,ipol,ifreq,*))
         tsys(ipol,iwin,ifreq)=mean(tsystcal(*,iwin,ipol,ifreq,*))
      endfor
   endfor
endfor
tsys=tsys*tcal

; write Tcal
;fdir='/cita/scratch/cottontail/tchang/GBT/Processed/'+dname+'/3c48_'+ffdir(ifile)+'/'
openw,1,fdir+'Tcal.txt'
printf,1,tcal
close,1

openw,1,fdir+'Tsys.txt'
printf,1,tsys
close,1

set_plot,'ps'
device,filename=fdir+'tcal.ps'

tcal=tcal<7.>(0.)

f=reform(freqtot(0,*))
t=reform(tcal(0,0,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(0,*))
t=reform(tcal(1,0,*))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
;
f=reform(freqtot(1,*))
t=reform(tcal(0,1,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(1,*))
t=reform(tcal(1,1,*))
oplot,f,t,linestyle=2
;
f=reform(freqtot(2,*))
t=reform(tcal(0,2,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(2,*))
t=reform(tcal(1,2,*))
oplot,f,t,linestyle=2
;
f=reform(freqtot(3,*))
t=reform(tcal(0,3,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(3,*))
t=reform(tcal(1,3,*))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
;
f=reform(freqtot(4,*))
t=reform(tcal(0,4,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(4,*))
t=reform(tcal(1,4,*))
oplot,f,t,linestyle=2
;
f=reform(freqtot(5,*))
t=reform(tcal(0,5,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(5,*))
t=reform(tcal(1,5,*))
oplot,f,t,linestyle=2
;
f=reform(freqtot(6,*))
t=reform(tcal(0,6,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(6,*))
t=reform(tcal(1,6,*))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
;
f=reform(freqtot(7,*))
t=reform(tcal(0,7,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
f=reform(freqtot(7,*))
t=reform(tcal(1,7,*))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_cal [K]'
;
device,/close

;
; plot Tsys
device,filename=fdir+'tsys.ps'

tsys=tsys<70.>(10.)

f=reform(freqtot(0,*))
t=reform(tsys(0,0,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
f=reform(freqtot(0,*))
t=reform(tsys(1,0,*))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
;
f=reform(freqtot(1,*))
t=reform(tsys(0,1,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
f=reform(freqtot(1,*))
t=reform(tsys(1,1,*))
oplot,f,t,linestyle=2
;
f=reform(freqtot(2,*))
t=reform(tsys(0,2,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
f=reform(freqtot(2,*))
t=reform(tsys(1,2,*))
oplot,f,t,linestyle=2
;
f=reform(freqtot(3,*))
t=reform(tsys(0,3,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
f=reform(freqtot(3,*))
t=reform(tsys(1,3,*))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
;
f=reform(freqtot(4,*))
t=reform(tsys(0,4,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
f=reform(freqtot(4,*))
t=reform(tsys(1,4,*))
oplot,f,t,linestyle=2
;
f=reform(freqtot(5,*))
t=reform(tsys(0,5,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
f=reform(freqtot(5,*))
t=reform(tsys(1,5,*))
oplot,f,t,linestyle=2
;
f=reform(freqtot(6,*))
t=reform(tsys(0,6,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
f=reform(freqtot(6,*))
t=reform(tsys(1,6,*))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
;
f=reform(freqtot(7,*))
t=reform(tsys(0,7,*))
plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
f=reform(freqtot(7,*))
t=reform(tsys(1,7,*))
oplot,f,t,linestyle=2
;plot,f,t,xtitle='!6freq [MHz]',ytitle='T_sys [K]'
;
device,/close

;endfor

end


; 1.0 GHz 21.37 JY
; 0.9 GHz 23.39 JY
; 0.8 GHz 25.87 JY
; 0.7 GHz 29.00 JY
; 0.6 GHz 33.09 JY
;
;
; flux = (nu(MHz)/1000)^(-0.85584) * 21.37  Jy  for 3c48
;

; 3c218
;freq_flux=[[1.40E+09, 4.08E+01],$
;           [9.60E+08, 6.54E+01],$
;           [7.50E+08, 8.36E+01],$
;           [6.35E+08, 9.71E+01],$
;           [4.68E+08, 1.15E+02],$
;           [4.08E+08, 1.32E+02],$
;;           [3.65E+08, 7.29E+01],$
;           [1.60E+08, 2.46E+02]]



