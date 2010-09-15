PRO combine_maps,corrv,corre,radio,weight,optical,field=field,fin=fin,radiocut=radiocut,fdir=fdir,noise=noise,sim=sim,daisy=daisy,outfile=outfile,dname=dname,svdfg=svdfg,mode=mode,dosigcut=dosigcut
seed=1
if not keyword_set(radiocut) then radiocut=0.1
if not keyword_set(fdir) then fdir=''
if not keyword_set(field) then field='f3'
;if not keyword_set(dname) then dname='sept07'
if not keyword_set(mode) then mode=''


if field eq 'f3' then begin
   ;
   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
   print,'thisdir: ',thisdir
   fdir=['/f3_53-58/','/f3_45-50/','/f3_27-29/']
   subdir=['/aug29/f3_05-10/', '/aug29/f3_15-20/',  '/aug29/f3_21-26/',  '/aug29/f3_27-32/',  $
           '/aug29/f3_37-42/',  '/aug29/f3_47-52/',  '/aug29/f3_53-58/', $
           '/aug30/f3_9-14/', '/aug30/f3_15-20/',  '/aug30/f3_21-26/',  '/aug30/f3_27-32/', $ 
           '/aug30/f3_33-38/',  '/aug30/f3_39-44/',  '/aug30/f3_45-50/', $ 
           '/sept07/f3_09-14/', '/sept07/f3_15-20/', $
           '/sept07/f3_21-26/', '/sept07/f3_27-29/']
   if not keyword_set(dname) then dname=['aug29','aug30','sept07']
   if n_elements(dname) eq 1 then begin
      if dname eq 'aug29' then begin
         fdir='/f3_53-58/'
         subdir=subdir(0:6)
      endif
      if dname eq 'aug30' then begin
         fdir='/f3_45-50/'
         subdir=subdir(7:13)
      endif
      if dname eq 'sept07' then begin
         fdir='/f3_27-29/'
         subdir=subdir(14:17)
      endif
   endif
   ffdir=thisdir+dname+fdir
   submapdir=thisdir+subdir
   nx=18
   ny=4
   ;
endif


if field eq 'f4' then begin
   ;
   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
   print,'thisdir: ',thisdir
   fdir=['/f4_69-74/','/f4_67-72/','/f4_33-36/']
   subdir=['/aug29/f4_63-68/', '/aug29/f4_69-74/','/aug30/f4_55-60/', $
           '/aug30/f4_61-66/', '/aug30/f4_67-72/',$
           '/sept06/f4_09-14/',  '/sept06/f4_15-20/',  '/sept06/f4_21-26/', $
           '/sept06/f4_27-32/',  '/sept06/f4_33-36/']
   if not keyword_set(dname) then dname=['aug29','aug30','sept06']
   if n_elements(dname) eq 1 then begin
      if dname eq 'aug29' then begin
         fdir='/f4_69-74/'
         subdir=subdir(0:1)
      endif
      if dname eq 'aug30' then begin
         fdir='/f4_67-72/'
         subdir=subdir(2:4)
      endif
      if dname eq 'sept06' then begin
         fdir='/f4_33-36/'
         subdir=subdir(5:9)
      endif
   endif
   ffdir=thisdir+dname+fdir
   submapdir=thisdir+subdir
   nx=12
   ny=4
   ;
endif
set_plot,'ps'
device,/color, bits_per_pixel=8

deepdir='/cita/scratch/cottontail/tchang/GBT/Deep2/'
;deepdir='/Users/tzu/Projects/GBT/Deep2/'
; field center
if (field eq 'f1') then begin
 ;  f0=[214.250,52.50]
   f0=[96.47,59.75]
   fbin=[2,14]
   fn=deepdir+'deep2_f1.dat'
   readcol,fn,ratmp,dectmp,oz,magi,format='d,d,d,d'
   euler,ratmp,dectmp,ora,odec,1
endif 
if (field eq 'f2') then begin
   f0=[252.43229,34.938242]     ;field center coordinate
   fsize=[1.43867,0.500210]    ;apparent size on the sky
   fbin=[12,4]
   fnopt=deepdir+'deep2_f2.dat'
   readcol,fnopt,ora,odec,oz,magi,format='d,d,d,d'
endif
if (field eq 'f3') then begin
   f0=[352.49077,0.13736937]
   fsize=[2.24112, 0.538381]
;   fbin=[18,5]
   fbin=[18,4]
   fnopt=deepdir+'deep2_f3.dat'
   readcol,fnopt,ora,odec,oz,magi,format='d,d,d,d'
endif
if (field eq 'f4') then begin
   f0=[37.255199,0.58632926]
   fbin=[12,4]
   fnopt=deepdir+'deep2_f4.dat'
   readcol,fnopt,ora,odec,oz,magi,format='d,d,d,d'
endif


; cross correlate GBT with the deep2 fields

; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
freq0=1420.405751 ;MHz
beam_sigma=beam/2.35482
dd=2.
nz=1120/dd
npol=4
nday=n_elements(dname)
gainxtot=dblarr(nday,npol,nx*ny)
gainztot=dblarr(nday,npol,nz)
Da_z,oz,dist
distbin=(floor((dist-d0)/dd)) ;d0 is at bin0

for iday=0,nday-1 do begin   
   ; read in gain_x and gain_z
   gainx=dblarr(npol,nx*ny)
   gainz=dblarr(npol,nz)
   ;
   gdir='~/projects/GBT/pros/'+dname(iday)+'/'+field+'/'
   openr,1,gdir+'deep2_gain_x_b4svd.txt'
   readf,1,gainx
   close,1
   ;
   openr,1,gdir+'deep2_gain_z_b4svd.txt'
   readf,1,gainz
   close,1
   ;
   gainxtot(iday,*,*)=gainx
   gainztot(iday,*,*)=gainz
   ;
endfor
;
;
; read in sub radio maps to do bootstrap error estimation
nsubmap=n_elements(submapdir)*npol
subrmap=dblarr(nsubmap,nx,ny,nz)
subwmap=dblarr(nsubmap,nx,ny,nz)
subweightmap=subwmap
if field eq 'f3' and n_elements(dname) eq 3 then dayindx=[0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2]
if field eq 'f3' and n_elements(dname) eq 1 then begin
   if dname eq 'aug29' then dayindx=[0,0,0,0,0,0,0]
   if dname eq 'aug30' then dayindx=[0,0,0,0,0,0,0]
   if dname eq 'sept07' then dayindx=[0,0,0,0]
endif
if field eq 'f4' and n_elements(dname) eq 3 then dayindx=[0,0,1,1,1,2,2,2,2,2]
if field eq 'f4' and n_elements(dname) eq 1 then begin
   if dname eq 'aug29' then dayindx=[0,0]
   if dname eq 'aug30' then dayindx=[0,0,0]
   if dname eq 'sept06' then dayindx=[0,0,0,0,0]
endif
;
;
nzz=370
radio=dblarr(nx,ny,nzz)
weight=dblarr(nx,ny,nzz)
;
for imap=0,nsubmap/npol-1 do begin
   ;
   rtmp=dblarr(npol,nx,ny,nz)
   wtmp=rtmp
   ;
   openr,1,submapdir(imap)+field+'.radio'+mode+'_weighted_aftsvd'
   readf,1,rtmp
   close,1
      ;
   openr,1,submapdir(imap)+field+'.weight'+mode+'_aftsvd'
   readf,1,wtmp
   close,1
      ;
      ;
   if field eq 'f3' then begin
      ;
      rtmp(*,*,3,*)=0
      wtmp(*,*,3,*)=0
      ;
      rtmp(*,13:nx-1,*,*)=0
      wtmp(*,13:nx-1,*,*)=0
      ;
   endif
   if field eq 'f4' then begin
      ;
      rtmp(*,*,0,*)=0
      wtmp(*,*,0,*)=0
      ;
   endif
   ;
   rtmp=reform(rtmp,npol,nx*ny,nz)
   wtmp=reform(wtmp,npol,nx*ny,nz)
   for ipol=0,npol-1 do begin
      for ix=0,nx*ny-1 do begin
         for iz=0,nz-1 do begin
            if (gainxtot(dayindx(imap),ipol,ix)*gainztot(dayindx(imap),ipol,iz) ne 0.) then begin
               rtmp(ipol,ix,iz)=rtmp(ipol,ix,iz)/gainxtot(dayindx(imap),ipol,ix)/gainztot(dayindx(imap),ipol,iz) 
            endif else begin
               wtmp(ipol,ix,iz)=0.
            endelse
         endfor
      endfor
   endfor
   ;
   rtmp=reform(rtmp,npol,nx,ny,nz)
   wtmp=reform(wtmp,npol,nx,ny,nz)
   ;
   rtmp2=rtmp(*,*,*,2+555/3:560-4)
   wtmp2=wtmp(*,*,*,2+555/3:560-4)
   for ipol=0,npol-1 do begin
      radio=radio+reform(rtmp2(ipol,*,*,*)*wtmp2(ipol,*,*,*))
      weight=weight+reform(wtmp2(ipol,*,*,*))
   endfor
   ;
endfor
;
nz=nzz
ind=where(weight gt 0,cind)
radio(ind)=radio(ind)/weight(ind)
;
openw,1,submapdir(imap-1)+field+'.radio'+mode+'_weighted_aftsvd_combined'
printf,1,radio
close,1
;
openw,1,submapdir(imap-1)+field+'.weight'+mode+'_aftsvd_combined'
printf,1,weight
close,1

return
;
; use last day's gain to set thermal threshould to kelvin
nfreq=2048
dnu=50e6/double(nfreq)
dt=0.5d
thermal=(dnu*dt)*(weight)
thermal=reform(thermal,nx*ny,nz)
gainx=dblarr(nx*ny)
gainz=dblarr(nz)
for ix=0,nx*ny-1 do begin
   gainx(ix)=mean(gainxtot(*,*,ix))
endfor
for iz=0,nz-1 do begin
   gainz(iz)=mean(gainztot(*,*,2+555/3+iz))
endfor
for ix=0,nx*ny-1 do begin
   for iz=0,nz-1 do begin
      if (gainx(ix)*gainz(iz) ne 0.) then begin
         thermal(ix,iz)=thermal(ix,iz)/gainx(ix)/gainz(iz) 
      endif
   endfor
endfor
;
neweight=dblarr(nx*ny,nz)
weight2d=reform(weight,nx*ny,nz)
radio2d=reform(radio,nx*ny,nz)
for iz=0,nz-1 do begin
   ind=where(weight2d(*,iz) gt 0,cind)
   if cind gt 1 then begin
      radiotmp=radio2d(*,iz)
      neweight(*,iz)=(1./variance(radiotmp(ind)))<thermal(*,iz)
   endif
endfor
neweight=reform(neweight,nx,ny,nz)
weight=neweight


;
nerror=100
;correlation=dblarr(nerror+1,nz)
;correlation_weight=dblarr(nerror+1,nz)
correlation2=dblarr(nerror+1)
correlation_weight2=correlation2
opticalreal=dblarr(nx,ny,nz)
optmaskreal=dblarr(nx,ny,nz)
lumweight=dblarr(nx,ny,nz)
nlag=61
lagarr=findgen(nlag)-30.
;correlation_lag=dblarr(nerror+1,nz,nlag)
correlation_lag=dblarr(nerror+1,nlag)

booterr_lag=dblarr(nerror+1,nlag)
boottot_lag=dblarr(nerror+1,nlag)

optical=dblarr(nx,ny,nz)
opticalerr=dblarr(nx,ny,nz)

; read in all sub maps to calculate the errors
; 
; make the optical map
for ierr=0,nerror do begin

   ; use magi as optical weight
   magi=magi-max(magi)
   lum=10^(-magi/2.5)
   
   if ierr eq 0 then begin
      ; bin the optical survey
      for i=0,n_elements(ora)-1 do begin
         ; figure out spatial bin
         x=(floor( (ora(i)-f0(0))*cos(odec(i)/!radeg) /dra)+nx/2) ;>0<(fbin(0)-1)
         y=(floor( (odec(i)-f0(1))/ddec) +ny/2)                   ;>0<(fbin(1)-1)
         if ierr eq 0 then zi=distbin(i) else begin
            testi=round(randomu(seed)*double(n_elements(ora)-1))
            zi=distbin(testi)
         endelse
         ;
         ; fill in the nearby cube slice
         if (zi ge 0 and zi lt nz) then begin
            for ii=-1, 1 do begin
               for jj=-1, 1 do begin
                  xi=x+ii
                  yi=y+jj
                  if xi ge 0 and xi lt nx and yi ge 0 and yi lt ny then begin
                     ycor=(yi-ny/2)*ddec+f0(1)
                     xcor=(xi-nx/2)*dra/cos(ycor/!radeg)+f0(0)
                     bdist= ( ((ora(i)-xcor)*cos(odec(i)/!radeg))^2. + $
                              (odec(i)-ycor)^2. )
                     optical(xi,yi,zi)=optical(xi,yi,zi)+exp(-bdist/beam_sigma^2./2.) ;*lum(i)
                     ;lumweight(xi,yi,zi)=lumweight(xi,yi,zi)+lum(i)
                  endif
               endfor
            endfor
         endif
         ;
      endfor
   endif
   ;
   ;
   ; convolve with the optical redshift space correlation, only in the
   ; redshift space
   ; kernel= exp(-sqrt(2)*|v|/sigma12) 
   ;
   if keyword_set(kernal) then begin
      sigma12=6./dd             ;Mpc
      temp=findgen(101)-50.
      kernel=exp(-sqrt(2.)*abs(temp)/sigma12)
      for ii=0,nx-1 do begin
         for jj=0,ny-1 do begin
            opttmp=reform(optical(ii,jj,*))
            result=convol(opttmp,kernel)
            optical(ii,jj,*)=result
         endfor
      endfor
   endif
   ;
   ;
   if ierr eq 0 then begin
      ;first calculate freq average
      optmask=dblarr(nx,ny,nz)
      optmask2d=dblarr(nx,ny)
      for i=0,nx-1 do begin
         for j=0,ny-1 do begin
            med=total(reform(optical(i,j,*)))
            if med ne 0 then med=1
            optmask2d(i,j)=med
            for k=0,nz-1 do begin
               optmask(i,j,k)=med
            endfor
         endfor
      endfor
      ;
      opticaltho=optical
      nnz=2
      noptz=nnz
      thobar=dblarr(noptz)
      doptz=nz/noptz
      index=indgen(nz)
      for i=0,noptz-1 do begin
         ind=where((index/doptz) eq i,cind)
         if cind ne doptz then print,'problem in optical rebinning!'
         ;
         ; ignore ny=3 dec strips here!!
         thobar(i)=mean(optical(*,*,ind))
         if field eq 'f3' then thobar(i)=mean(optical(*,0:ny-1-1,ind))
         if field eq 'f4' then thobar(i)=mean(optical(*,1:ny-1,ind))
         if ierr eq 0 and thobar(i) eq 0. then begin
            print,'cind',cind
            print,'thobar eq 0',thobar
            print,'thobar eq 0'
            ;read,test
            thobar(i)=1.
         endif
         opticaltho(*,*,ind)=opticaltho(*,*,ind)/thobar(i)-1.
      endfor
      optical=opticaltho
   endif

   if ierr eq 0 then optmaskreal=optmask

   if ierr eq 0 then begin
      openw,1,ffdir(nday-1)+field+'.opticaltho'
      printf,1,optical
      close,1
   endif


   ; subtract off spatial average                                                  
   opt1d=dblarr(nz)
   for i=0,nz-1 do begin
   ;
   ; ignore ny=3 dec strip here!!
      optslice=reform(optical(*,*,i))
      if field eq 'f3' then optslice=reform(optical(*,0:ny-1-1,i))
      if field eq 'f4' then optslice=reform(optical(*,1:ny-1,i))
      opt1d(i)=total(optslice*optmask2d)/total(optmask2d)
      optical(*,*,i)=optical(*,*,i)-opt1d(i)
   endfor

   if ierr eq 0 then opticalreal=optical

   if ierr eq 0 then begin
      openw,1,ffdir(nday-1)+field+'.opticalthobar'
      printf,1,optical
      close,1
      openw,1,ffdir(nday-1)+field+'.opticalthobar_weight'
      printf,1,weight
      close,1
   endif

   radioerr=dblarr(nx,ny,nz)
   weighterr=dblarr(nx,ny,nz)

   ;noptz=37  ; seems to be a mistake here, noptz should be 10 for a
   ;74Mpc strip for correlation.  06/16/09
   noptz=10
   dnoptz=(nz/noptz)
   npick=nz-dnoptz-1
   for ii=0,noptz-1 do begin
      ; randomly pick a number between 0 to nz-dnoptz-1
      nstart=round(randomu(seed)*double(npick))
      radioerr(*,*,ii*dnoptz:((ii+1)*dnoptz-1))=radio(*,*,nstart:(nstart+dnoptz-1))
      weighterr(*,*,ii*dnoptz:((ii+1)*dnoptz-1))=weight(*,*,nstart:(nstart+dnoptz-1))
   endfor

   ;make a random shuffle of the radio cube
   radioran=dblarr(nx,ny,nz)
   weightran=dblarr(nx,ny,nz)
   for iz=0,nz-1 do begin
      testi=round(randomu(seed)*double(nz-1))
      radioran(*,*,iz)=radio(*,*,testi)
      weightran(*,*,iz)=weight(*,*,testi)
   endfor

   ;if ierr eq 0 then opticalerr=opticalreal
   opticalerr=opticalreal
   opterrmask=optmaskreal
   if ierr eq 0 then radioerr=radio
   if ierr eq 0 then radioran=radio
   if ierr eq 0 then weighterr=weight
   if ierr eq 0 then weightran=weight
 
; NOTE to self:
; opitcalreal is the real optical catalog in 3D
; opticalerr is the piece-wise shuffled real optical calatog. if
; ierr=0 opticalerr=opticalreal 
; optical is the random, individually shuffled real optical
; catalog. if ierr=0, optical=opticalreal

   ; cross-correlation that includes potential radio fg correlation
   corre=radioerr*opticalreal

   ; bootstrap correlation
   correboots=radioran*opticalreal

   ; with radio fg correlation, all z
   correlation2(ierr)=total(corre*weighterr*opterrmask)/total(weighterr*opterrmask)

   ; bootstrap, all z
   correlation_weight2(ierr)=total(correboots*weightran*optmask)/total(weightran*optmask)
   ;

   ; weight by redshift
;   for ii=0,nz-1 do begin
;      ind=where(weighterr(*,*,ii)*opterrmask(*,*,ii) gt 0,cind1)
;      ind=where(weightran(*,*,ii)*optmaskreal(*,*,ii) gt 0,cind2)
;      if cind1 gt 0 then $
;         correlation(ierr,ii)=total(corre(*,*,ii)*weighterr(*,*,ii)*opterrmask(*,*,ii)) $
;                              /total(weighterr(*,*,ii)*opterrmask(*,*,ii))
;      if cind2 gt 0 then $
;         correlation_weight(ierr,ii)=total(correboots(*,*,ii)*weightran(*,*,ii)*optmaskreal(*,*,ii;)) $
;                                     /total(weightran(*,*,ii)*optmaskreal(*,*,ii))
;   endfor



   if ierr eq 0 then begin
      corrold=total(corre*weight*optmaskreal)/total(weight*optmaskreal)
      print,'original correlation:',corrold
   endif

   ; calculate correlation fn at different lag
   for ii=0,nlag-1 do begin
      radiotmp=dblarr(nx,ny,nz)
      weighttmp=dblarr(nx,ny,nz)
      lag=long(lagarr(ii))
      if lag ge 0 then begin
         radiotmp(*,*,0:nz-1-lag)=radioerr(*,*,lag:nz-1)
         weighttmp(*,*,0:nz-1-lag)=weighterr(*,*,lag:nz-1)
      endif
      if lag lt 0 then begin
         radiotmp(*,*,-lag:nz-1)=radioerr(*,*,0:nz-1+lag)
         weighttmp(*,*,-lag:nz-1)=weighterr(*,*,0:nz-1+lag)
      endif
;      for jj=0,nz-1 do begin
      ind=where(weighttmp*opterrmask gt 0,cind)
      if cind gt 0 then begin
         correlation_lag(ierr,ii)=total(radiotmp*opticalerr*weighttmp*opterrmask) $
                                  /total(weighttmp*opterrmask)
      endif
   endfor
   ;
   ;
   ; bootstrap all the submaps
      ; calculate correlation fn at different lag
   booterr=dblarr(nsubmap,nlag)
   rboottot=dblarr(nx,ny,nz)
   wboottot=rboottot
   for imap=0,nsubmap-1 do begin
      if ierr eq 0 then num=imap else num=round(randomu(seed)*double(nsubmap-1))
      rboottot=rboottot+reform(subrmap(num,*,*,*)*subweightmap(num,*,*,*))
      wboottot=wboottot+reform(subweightmap(num,*,*,*))
      rboot=reform(subrmap(num,*,*,*))
      wboot=reform(subwmap(num,*,*,*))
      ;
      for ii=0,nlag-1 do begin
         radiotmp=dblarr(nx,ny,nz)
         weighttmp=dblarr(nx,ny,nz)
         lag=long(lagarr(ii))
         if lag ge 0 then begin
            radiotmp(*,*,0:nz-1-lag)=rboot(*,*,lag:nz-1)
            weighttmp(*,*,0:nz-1-lag)=wboot(*,*,lag:nz-1)
;            radiotmp(*,*,nz-lag:nz-1)=rboot(*,*,0:lag-1)
;            weighttmp(*,*,nz-lag:nz-1)=wboot(*,*,0:lag-1)
         endif
         if lag lt 0 then begin
            radiotmp(*,*,-lag:nz-1)=rboot(*,*,0:nz-1+lag)
            weighttmp(*,*,-lag:nz-1)=wboot(*,*,0:nz-1+lag)
         endif
         booterr(imap,ii)=total(radiotmp*opticalerr*weighttmp*opterrmask) $
                          /total(weighttmp*opterrmask)
      endfor
   endfor
   for ii=0,nlag-1 do begin
      booterr_lag(ierr,ii)=mean(booterr(*,ii))
   endfor
   ;
   if ierr eq 0 then print,'submap correlations',booterr(*,0)
   ;
   if ierr eq 0 then print,'submap correlation',booterr_lag(0,0)
   ;
                                ; calculate correlation fn at
                                ; different lag for the bootsrapped
                                ; maps
   ind=where(wboottot gt 0,cind)
   rboottot(ind)=rboottot(ind)/wboottot(ind)
   ;
   thermal=(dnu*dt)*(wboottot)
   thermal=reform(thermal,nx*ny,nz)
   gainx=dblarr(nx*ny)
   gainz=dblarr(nz)
   for ix=0,nx*ny-1 do begin
      gainx(ix)=mean(gainxtot(*,*,ix))
   endfor
   for iz=0,nz-1 do begin
      gainz(iz)=mean(gainztot(*,*,iz))
   endfor
   for ix=0,nx*ny-1 do begin
      for iz=0,nz-1 do begin
         if (gainx(ix)*gainz(iz) ne 0.) then begin
            thermal(ix,iz)=thermal(ix,iz)/gainx(ix)/gainz(iz) 
         endif
      endfor
   endfor
   ;
   thermal=reform(thermal,nx,ny,nz)
   for iz=0,nz-1 do begin
      ind=where(wboottot(*,*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=rboottot(*,*,iz)
         wboottot(*,*,iz)=(1./variance(radiotmp(ind)))<thermal(*,*,iz)
      endif
   endfor
   for ii=0,nlag-1 do begin
      radiotmp=dblarr(nx,ny,nz)
      weighttmp=dblarr(nx,ny,nz)
      lag=long(lagarr(ii))
      if lag ge 0 then begin
         radiotmp(*,*,0:nz-1-lag)=rboottot(*,*,lag:nz-1)
         weighttmp(*,*,0:nz-1-lag)=wboottot(*,*,lag:nz-1)
      endif
      if lag lt 0 then begin
         radiotmp(*,*,-lag:nz-1)=rboottot(*,*,0:nz-1+lag)
         weighttmp(*,*,-lag:nz-1)=wboottot(*,*,0:nz-1+lag)
      endif
;      for jj=0,nz-1 do begin
      ind=where(weighttmp*opterrmask gt 0,cind)
      if cind gt 0 then begin
         boottot_lag(ierr,ii)=total(radiotmp*opticalerr*weighttmp*opterrmask) $
                                  /total(weighttmp*opterrmask)
      endif
   endfor
   ;
   if ierr eq 0 then print,'submaptot correlation',boottot_lag(0,0)
   ;
endfor




; cooadd the correlation
;correrr=dblarr(nz)
;correrrslice=dblarr(nz)
;for ii=0,nz-1 do begin
;   weightslice=weighterr(*,*,ii)
;   ind=where(weightslice gt 0,cind)
;   if cind gt 0 then begin
;      correrrslice(ii)=variance(correlation(1:nerror,ii))
;   endif
;   weightslice=weightran(*,*,ii)
;   ind=where(weightslice gt 0,cind)
;   print,'weightslice gt 0',cind
;   if cind gt 0 then begin
;      correrr(ii)=(variance(correlation_weight(1:nerror,ii)))
;      print,'correrr(ii)',ii,correrr(ii)
;   endif
;endfor
;
;
;ind=where(correrr ne 0 and finite(correrr),cind)
;print,'correrr:cind',cind
;ind3=where(correrrslice ne 0 and finite(correrrslice),cind)
;print,'correrrslice:cind',cind
;
;ind2=where(1./correrr(ind) lt 1e9)
;ind4=where(1./correrrslice(ind3) lt 1e9)
;
;xcorr0=total(correlation_weight(0,ind(ind2))/correrr(ind(ind2)))/total(1./correrr(ind(ind2)))
;xcorr=total(correlation(0,ind3(ind4))/correrrslice(ind3(ind4)))/total(1./correrrslice(ind3(ind4)))
;corr=dblarr(nerror)
;corr2=dblarr(nerror)
;corr3=corr
;for i=1,nerror do begin
;   corr(i-1)=total(correlation(i,ind3(ind4))/correrrslice(ind3(ind4)))/total(1./correrrslice(ind3(ind4)))
;   corr2(i-1)=total(correlation_weight(i,ind(ind2))/correrr(ind(ind2)))/total(1./correrr(ind(ind2)))
;endfor
;

corrlagvalue=dblarr(nlag)
corrlagerr=dblarr(nlag)
corrlagerr_mean=dblarr(nlag)
corrcovar=dblarr(nlag,nlag)
booterr_mean=dblarr(nlag)
booterr_err=dblarr(nlag)
boottot_mean=dblarr(nlag)
boottot_err=dblarr(nlag)

;xcorrlag=dblarr(nerror+1,nlag)
;for kk=0,nlag-1 do begin
;   correrrlag=dblarr(nz)
;   for ii=0,nz-1 do begin
;      ind=where(weight(*,*,ii) gt 0,cind)
;      if cind gt 0 then begin
;         correrrlag(ii)=variance(correlation_lag(1:nerror,ii,kk))
;      endif
;   endfor
;   ind=where(1./correrrlag lt 1e9)
;   for jj=0,nerror do begin
;      xcorrlag(jj,kk)=total(correlation_lag(jj,ind,kk)/correrrlag(ind))/total(1./correrrlag(ind))
;   endfor
;   corrlagvalue(kk)=xcorrlag(0,kk)
;   corrlagerr(kk)=sqrt(total(xcorrlag(1:nerror,kk)^2.)/double(nerror))
;   corrlagerr_mean(kk)=mean(xcorrlag(1:nerror,kk))
   ;
;endfor
;corrlagvalue=corrlagvalue-corrlagerr_mean
;

for kk=0,nlag-1 do begin
   ;
   corrlagvalue(kk)=correlation_lag(0,kk)
   corrlagerr_mean(kk)=mean(correlation_lag(1:nerror,kk))
   corrlagerr(kk)=sqrt(total((correlation_lag(1:nerror,kk)-corrlagerr_mean(kk))^2.))/double(nerror)
   booterr_mean(kk)=mean(booterr_lag(1:nerror,kk))
   booterr_err(kk)=sqrt(total((booterr_lag(1:nerror,kk)-booterr_lag(0,kk))^2.))/double(nerror)
   ;
   boottot_mean(kk)=mean(boottot_lag(1:nerror,kk))
   boottot_err(kk)=sqrt(total((boottot_lag(1:nerror,kk)-boottot_lag(0,kk))^2.))/double(nerror)
   ;
endfor
corrlagvalue=corrlagvalue-corrlagerr_mean
corrlagvalue_boot=corrlagvalue-booterr_mean
;
;
;calculate the noise covariance matrix
for jj=0,nlag-1 do begin
   for kk=0,nlag-1 do begin
      corrcovar(jj,kk)=total(correlation_lag(1:nerror,jj)*correlation_lag(1:nerror,kk))/(double(nerror))
   endfor
endfor

xerr=booterr_err(0)
xerr_mean=booterr_mean(0)
xerr3=boottot_err(0)
xerr3_mean=boottot_mean(0)
;xerr=sqrt(variance(corr))
;xerr_mean=mean(corr)
xerr2=sqrt(variance(correlation2(1:nerror)))
xerr2_mean=mean(correlation2(1:nerror))
;xerr3=sqrt(variance(corr2))
;xerr3_mean=mean(corr2)
xerr4=sqrt(variance(correlation_weight2(1:nerror)))
xerr4_mean=mean(correlation_weight2(1:nerror))
;
;
;print,'redshift-weighted correlation with fg-weight',xcorr,xcorr-xerr_mean,xerr_mean

;corrv=xcorr-xerr_mean
;corre=xerr

;print,'redshift-weighted correlation',xcorr0,xcorr0-xerr3_mean,xerr3_mean
;print,'error of correlation with fg-correlation and z-weight',xerr
;print,'error of correlation with bootstrap and z-weight',xerr3
print,'original correlation:',corrold,corrold-xerr2_mean,corrold-xerr4_mean
print,'mean of randomized correlation',xerr2_mean,xerr4_mean,xerr_mean,xerr3_mean
print,'error of correlation with fg-correlation and no z-weight',xerr2
print,'error of correlation with bootstrap and no z-weight',xerr4
print,'bootstrap error',xerr
print,'bootstraptot error',xerr3
print,'rms HI after sigmacut:',sqrt(total(radio^2.*weight^2.)/total(weight^2.))
;print,'3-sigma flag:',cind2,double(cind2)/double(nx*ny*double(nz))
;
print,''
print,'correlation at different lag:',corrlagvalue
print,'correlatioin error at different lag:',corrlagerr
;
;print,'covariance',corrcovar
;print,'diagonal:'
;for ii=0,nlag-1 do print,sqrt(corrcovar(ii,ii))
openw,1,'corrlag.radio.dat' ,/append
printf,1,corrlagvalue
printf,1,corrlagerr
close,1

openw,1,'corrcovar.radio.dat' ,/append
printf,1,corrcovar
close,1
;printf,5,'redshift-weighted correlation',xcorr
;printf,5,'error of correlation',xerr
;printf,5,'error of correlation with normal weight',xerr2
;printf,5,'original correlation:',corrold
;printf,5,'error of correlation',xerr3
;printf,5,'error of correlation with normal weight',xerr4

;close,5

if keyword_set(outfile) then begin

openw,10,'correlation_orig.dat',/append
printf,10,corrold
close,10

openw,10,'correlation_error.dat',/append
printf,10,xerr2
close,10

endif


END

      
      
