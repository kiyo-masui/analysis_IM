pro cal_tcal_matrix,dname=dname,ffdir=ffdir,field=field

; calculate T_cal(nu) as a function of the binned (x,z) matrix
; -- basically take T_cal(nu) and convert it into T_cal(z) 

if not keyword_set(field) then field='3c218'  ; for the radio/optical cube setup
if not keyword_set(dname) then dname='03'

if not keyword_set(ffdir) then begin
  if dname eq '00' then begin
     ffdir=['5-8'];,'34-37']
  endif

  if dname eq '02' then begin
     ffdir=['5-8']
  endif

  if dname eq '03' then begin
     ffdir=['5-8']
  endif

  if dname eq '04' then begin
     ffdir=['5-8']
  endif
endif

outdirroot = '/mnt/raid-project/gmrt/kiyo/wiggleZ/calibration/'

;
; parameters
nscan=4
npol=4
nwin=8
nfreq=2048

fcaldir = ffdir
;read in Tcal
ffcaldir=outdirroot+dname+'/'+field+'_'+fcaldir(0)+'/'
tcal=dblarr(2,nwin,nfreq)
openr,1,ffcaldir+'Tcal.txt'
readf,1,tcal
close,1

; tsys
tsys=dblarr(2,nwin,nfreq)
openr,1,ffcaldir+'Tsys.txt'
readf,1,tsys
close,1

; calculate the corresponding frequencies
nfile=1
nkeep=nfreq*3/5
freqtot=dblarr(nwin*nkeep)
tcaltot=dblarr(2,nwin*nkeep)
tsystot=dblarr(2,nwin*nkeep)

for ifile=0,nfile-1 do begin
   fdir=outdirroot+dname+'/'+field+'_'+fcaldir(ifile)+'/'
   for iwin=0,nwin-1 do begin
      ; take 1/scan0-1/scan1=T_3c48/Tcal,
      ; take scan average of (0,1) and (2,3)
      tcaltsys=dblarr(nscan,npol,nfreq)
      freq=dblarr(nfreq)
      openr,1,fdir+'Freq.scan0'+strtrim(string(iwin),2)+'.txt'
      readf,1,freq
      close,1
      ; arrange in one big long frequency array
      freqtot((iwin*nkeep):(iwin+1)*nkeep-1)=freq(nfreq/5+1:nfreq*4/5-1)
      tcaltot(*,(iwin*nkeep):(iwin+1)*nkeep-1)=tcal(*,iwin,nfreq/5+1:nfreq*4/5-1)
      tsystot(*,(iwin*nkeep):(iwin+1)*nkeep-1)=tsys(*,iwin,nfreq/5+1:nfreq*4/5-1)
   endfor
endfor

; change to redshift domain
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz
Da_z,(freq0/freqtot-1.),distslice
distbin=(floor((distslice-d0)/dd)) ;d0 is at bin0
;
tcalz=dblarr(2,nz)
tsysz=dblarr(2,nz)
weightz=dblarr(nz)
freqbin=dblarr(nz)
for i=0,nz-1 do begin
   ind=where(distbin eq i,cind)
   if cind gt 0 then begin
      tcaltmp=tcaltot(0,ind)
      tcalz(0,i)=median(tcaltmp,/double)
      tcaltmp=tcaltot(1,ind)
      tcalz(1,i)=median(tcaltmp,/double)
      tsystmp=tsystot(0,ind)
      tsysz(0,i)=median(tsystmp,/double)
      tsystmp=tsystot(1,ind)
      tsysz(1,i)=median(tsystmp,/double)
      weightz(i)=double(cind)
      freqbin(i)=median(freqtot(ind),/double)
   endif
;      tcalz(0,i)=total(tcaltot(0,ind))/double(cind)
;      tcalz(0,i)=total(tcaltot(1,ind))/double(cind)
endfor

;print,'freqbin:',freqbin

spawn, 'mkdir -p '+dname

openw,1,dname+'/tcal_z.txt'
printf,1,tcalz
close,1
;
openw,1,dname+'/tsys_z.txt'
printf,1,tsysz
close,1
;
openw,1,dname+'/tcal_weight.txt'
printf,1,weightz
close,1
;
openw,1,dname+'/tcal_freq.txt'
printf,1,freqbin
close,1
;
set_plot,'ps'
device,filename=dname+'/tcal_z.ps'
ind=where(freqbin ne 0.)
plot,freqbin(ind),reform(tcalz(0,ind)),xtitle='!6freq [MHz]',ytitle='T!dcal!n [K]'
oplot,freqbin(ind),reform(tcalz(1,ind)),linestyle=2
device,/close

set_plot,'ps'
device,filename=dname+'/tsys_z.ps'
ind=where(freqbin ne 0.)
plot,freqbin(ind),reform(tsysz(0,ind)),xtitle='!6freq [MHz]',ytitle='T!dsys!n [K]'
oplot,freqbin(ind),reform(tsysz(1,ind)),linestyle=2
device,/close
end

;useful freqs:        6480
;         555         172
;       838.31055       680.03418
;tsys_median:       27.533531       30.636456
;mean tsys_median       29.084993
;tsys_mean       29.010866       31.705120
;tsys_median:       25.845503       28.719227



; 1.0 GHz 21.37 JY
; 0.9 GHz 23.39 JY
; 0.8 GHz 25.87 JY
; 0.7 GHz 29.00 JY
; 0.6 GHz 33.09 JY
;
;
; flux = (nu(MHz)/1000)^(-0.85584) * 21.37  Jy
;

