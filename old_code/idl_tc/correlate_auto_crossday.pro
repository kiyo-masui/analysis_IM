PRO correlate_auto_crossday,auto,meanauto,err,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,radio,weight,optical,field=field,fin=fin,radiocut=radiocut,fdir=fdir,noise=noise,sim=sim,daisy=daisy,outfile=outfile,dname=dname,svdfg=svdfg,mode=mode,combine=combine,dosigcut=dosigcut,svdmode=svdmode

seed=1
if not keyword_set(radiocut) then radiocut=0.1
if not keyword_set(fdir) then fdir=''
if not keyword_set(field) then field='1hr'
;if not keyword_set(dname) then dname='sept07'
if not keyword_set(mode) then mode=''
if not keyword_set(svdmode) then svdmode=''
if not keyword_set(half) then half='_half'

thisdir='/cita/d/raid-cita/tchang/wiggleZ/processed/'
if field eq '1hr' then begin
   fdir=['/june14/97-104/','/june14/153-160/']
   dname=['00','00']
endif
if field eq '22hr' then begin
   fdir=['/june14/25-32/','/june14/41-48/']
   dname=['00','00']
endif
ffdir=thisdir+fdir

set_plot,'ps'
device,/color, bits_per_pixel=8

; cross correlate GBT with the deep2 fields

; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz
beam_sigma=beam/2.35482

if not keyword_set(fbin) then fbin=[40,20]
npol=4
nx=fbin(0)
ny=fbin(1)
nday=n_elements(ffdir)
;
radionday=dblarr(nday,nx*ny,nz)
weightday=dblarr(nday,nx*ny,nz)
;
radioaddtot=dblarr(nx*ny,nz)
weightaddtot=dblarr(nx*ny,nz)
gainxtot=dblarr(nday,npol,nx*ny)
gainztot=dblarr(nday,npol,nz)
;
;
;
; read in data
for iday=0,nday-1 do begin
   ;
   niday=strtrim(string(iday),2)
   radio=dblarr(npol,nx,ny,nz)
   weight=radio
   radioadd=dblarr(nx*ny,nz)
   weightadd=dblarr(nx*ny,nz)
   ;
   openr,1,ffdir(iday)+field+'.radiotot'+mode+'_weighted_aftsvd'+half
   readf,1,radio
   close,1
   ;
   openr,1,ffdir(iday)+field+'.weighttot'+mode+'_aftsvd'+half
   readf,1,weight
   close,1
   ;
   ;
   ; apply calibratoin
   ; reform the arrays
   radio=reform(radio,npol,nx*ny,nz)
   weight=reform(weight,npol,nx*ny,nz)
   ;
   help,weight
   help,radio
   ;
   ; read in gain_x and gain_z
   gainx=dblarr(npol,nx*ny)
   gainz=dblarr(npol,nz)
   ;
   gdir='~/projects/GBT/pros/wiggle_z/Calibration/'+dname(iday)+'/'
   openr,1,gdir+field+'.wigglez_gain_x_b4svd.txt'
   readf,1,gainx
   close,1
   ;
   openr,1,gdir+field+'.wigglez_gain_z_b4svd.txt'
   readf,1,gainz
   close,1
   ;
   ;
   radio=reform(radio,npol,nx,ny,nz)
   weight=reform(weight,npol,nx,ny,nz)
   ;
   ;
   radioall=radio
   weightall=weight
   neweight=dblarr(npol,nx,ny,nz)
   nfreq=2048
   dnu=50e6/double(nfreq)
   dt=0.5d
   thermal=(dnu*dt)*(weight)
   for ipol=0,npol-1 do begin
      for iz=0,nz-1 do begin
         ind=where(weight(ipol,*,*,iz) gt 0,cind)
         if cind gt 1 then begin
            radiotmp=radio(ipol,*,*,iz)
            neweight(ipol,*,*,iz)=(1./variance(radiotmp(ind)));<thermal(ipol,*,*,iz)
         endif
      endfor
   endfor
   neweightall=neweight
   ;
   ;
   ;
   for ipol=0,Npol-1 do begin
      ;
      radio=reform(radioall(ipol,*,*,*))
      weight=reform(weightall(ipol,*,*,*))
      neweight=reform(neweightall(ipol,*,*,*))
      print,'max(radio),min(radio),mean(radio),variance(radio)'
      print,max(radio),min(radio),mean(radio),variance(radio)
      ind=where(weight gt 0,cind)
      ;
      ; do a svd on the spatial-redsfhit space
      ;
      if keyword_set(svdfg) then begin
         ;
         radiosvd=dblarr(nx*ny,nz)
         radiores=radiosvd
         radiosvd1mode=dblarr(nx*ny,nz)
         radiores1mode=radiosvd
         radiosvd2mode=dblarr(nx*ny,nz)
         radiores2mode=radiosvd
         radio2dtot=reform(radio,nx*ny,nz)
         weight2dtot=reform(weight,nx*ny,nz)
         neweight2dtot=reform(neweight,nx*ny,nz)
         ;
         ;do weighting
         radio2dtot=radio2dtot*neweight2dtot
         ;
         la_svd,radio2dtot, s,u,v,/double,status=status
         print,'SVD status:',status
         ;
         help,s
         help,u
         help,v
         ;
         w=dblarr(n_elements(s))
         help,w
         w2=w
         w(0)=s(0)
         svdm=long(svdmode)
         w2(0:(svdm-1))=s(0:(svdm-1))
         ;
         radiosvd=u ## diag_matrix(s) ## transpose(v)
         radiosvd1mode=u ## diag_matrix(w) ## transpose(v)
         radiosvd2mode=u ## diag_matrix(w2) ## transpose(v)
         ;
         radiores=radio2dtot-radiosvd
         radiores1mode=radio2dtot-radiosvd1mode
         radiores2mode=radio2dtot-radiosvd2mode
         ;
         ;
         ;undo weighting
         radiores1mode=radiores1mode/neweight2dtot
         radiores2mode=radiores2mode/neweight2dtot
         radio2dtot=radio2dtot/neweight2dtot
         ;
         ; do another 3-sigma cut
         sigcut=3.
         nfreq=2048
         dnu=50e6/double(nfreq)
         dt=0.5d
         sigmal=1./sqrt(dnu*dt)
         ;
         ind=where(weight2dtot gt 0,cind)
         sigmacut=sigcut*sigmal/sqrt(weight2dtot(ind))
         ;
         ind2=where(abs(radiores2mode(ind)) lt sigcut*sigmal/sqrt(weight2dtot(ind)),cind2)
         print,'nx,ny,nz',nx,ny,nz,double(nx*ny*double(nz))
         print,'3-sigma flag:',cind2,double(cind2)/double(nx*ny*double(nz))
         ;
         radio=dblarr(nx*ny,nz)
         weight=radio
         radio(ind(ind2))=radiores2mode(ind(ind2))
         weight(ind(ind2))=weight2dtot(ind(ind2))

         ;
      endif else begin
         ;
         radio=reform(radio,nx*ny,nz)
         weight=reform(weight,nx*ny,nz)
         ;
      endelse
      ;
      if keyword_set(dosigcut) then begin
      ; do another 3-sigma cut
      sigcut=3.
      nfreq=2048
      dnu=50e6/double(nfreq)
      dt=0.5d
      sigmal=1./sqrt(dnu*dt)
      ;
      ind=where(weight gt 0,cind)
      sigmacut=sigcut*sigmal/sqrt(weight(ind))
      ;
      ind2=where(abs(radio(ind)) gt sigcut*sigmal/sqrt(weight(ind)),cind2)
      print,'nx,ny,nz',nx,ny,nz,double(nx*ny*double(nz))
      print,'3-sigma flag:',cind2,double(cind2)/double(nx*ny*double(nz))
      ;
      if cind2 gt 0 then begin
         radio(ind(ind2))=0.
         weight(ind(ind2))=0
      endif
      ;
   endif
      ;
      help,radio
      help,weight
      help,gainx
      help,gainz
      ; set back to kelvin
      for ix=0,nx*ny-1 do begin
         for iz=0,nz-1 do begin
            if (gainx(ipol,ix)*gainz(ipol,iz) ne 0.) then begin
               radio(ix,iz)=radio(ix,iz)/gainx(ipol,ix)/gainz(ipol,iz)
            endif else begin
               weight(ix,iz)=0.
            endelse
         endfor
      endfor
      ;
      ;
      radioadd=radioadd+radio*weight
      weightadd=weightadd+weight
      radioaddtot=radioaddtot+radio*weight
      weightaddtot=weightaddtot+weight
      ;
   endfor
   ;
   ind=where(weightadd gt 0,cind)
   radioadd(ind)=radioadd(ind)/weightadd(ind)
   radionday(iday,*,*)=radioadd
   weightday(iday,*,*)=weightadd
   ;
endfor
;
help,gainx
help,gainz
;
radioday=reform(radionday,nday,nx,ny,nz);   Tsys=40K - change the unit to Kelvin
weightday=reform(weightday,nday,nx,ny,nz)

ind=where(weightaddtot gt 0,cind)
radioaddtot(ind)=radioaddtot(ind)/weightaddtot(ind)
;
;
;
; use last day's gain to set thermal threshould to kelvin
thermal=(dnu*dt)*(weightday)
thermal=reform(thermal,nday,nx*ny,nz)
for iday=0,nday-1 do begin
   gainx=dblarr(nx*ny)
   gainz=dblarr(nz)
   for ix=0,nx*ny-1 do begin
      gainx(ix)=mean(gainxtot(iday,*,ix))
   endfor
   for iz=0,nz-1 do begin
      gainz(iz)=mean(gainztot(iday,*,iz))
   endfor
   for ix=0,nx*ny-1 do begin
      for iz=0,nz-1 do begin
         if (gainx(ix)*gainz(iz) ne 0.) then begin
            thermal(iday,ix,iz)=thermal(iday,ix,iz)/gainx(ix)/gainz(iz) 
         endif
      endfor
   endfor
endfor
;
neweight=dblarr(nday,nx*ny,nz)
weightx=dblarr(nday,nx*ny,nz)
for iday=0,nday-1 do begin
   weight2d=reform(weightday(iday,*,*,*))
   weight2d=reform(weight2d,nx*ny,nz)
;   radio2d=reform(radioday(iday,*,*,*))
;   radio2d=reform(radio2d,nx*ny,nz)
   radio2d=radioaddtot
   for iz=0,nz-1 do begin
      ind=where(weight2d(*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=radio2d(*,iz)
         neweight(iday,*,iz)=(1./variance(radiotmp(ind)));<thermal(iday,*,iz)
      endif
   endfor
endfor
;
;
;
weightday=reform(weightday,nday,nx*ny,nz)
for iday=0,nday-1 do begin
   for ix=0,nx*ny-1 do begin
;      weightx(iday,ix,*)=total(weightday(iday,ix,*))
       weightx(iday,ix,*)=total(weightaddtot(ix,*))
   endfor
endfor
neweight=neweight*weightx
neweight=reform(neweight,nday,nx,ny,nz)
;weight=neweight
;
;
auto=dblarr(nday*(nday-1)/2)
count=0
for iday=0,nday-1 do begin
   r1=reform(radioday(iday,*,*,*))
   w1=reform(neweight(iday,*,*,*))
   for jday=iday+1,nday-1 do begin
      r2=reform(radioday(jday,*,*,*))
      w2=reform(neweight(jday,*,*,*))
      auto(count)=total(r1*r2*w1*w2)/total(w1*w2)
      count=count+1
   endfor
endfor

radio=reform(radioday(0,*,*,*)*neweight(0,*,*,*))
radio=reform(radio,nx*ny,nz)
; make a plot
device,filename=field+'autosvd.'+mode+svdmode+'.ps'
plt_image, transpose(radio),/scalable,/colbar
device,/close


if keyword_set(error) then begin
;
;
;
nerror=1000
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



; read in all sub maps to calculate the errors
; 
; make the optical map
for ierr=0,nerror do begin

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
   opticalerr=radio
   opterrmask=weight
   opticalreal=radio
   optmask=weight
   optmaskreal=weight
;   opticalerr=opticalreal
;   opterrmask=optmaskreal
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
         radiotmp(*,*,0:nz-1-lag)=radioran(*,*,lag:nz-1)
         weighttmp(*,*,0:nz-1-lag)=weightran(*,*,lag:nz-1)
      endif
      if lag lt 0 then begin
         radiotmp(*,*,-lag:nz-1)=radioran(*,*,0:nz-1+lag)
         weighttmp(*,*,-lag:nz-1)=weightran(*,*,0:nz-1+lag)
      endif
;      for jj=0,nz-1 do begin
      ind=where(weighttmp*opterrmask gt 0,cind)
      if cind gt 0 then begin
         correlation_lag(ierr,ii)=total(radiotmp*opticalreal*weighttmp*opterrmask) $
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
;   if keyword_set(boot) then begin
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
;endif
   for ii=0,nlag-1 do begin
      booterr_lag(ierr,ii)=mean(booterr(*,ii))
   endfor
   ;
   if ierr eq 0 then print,'submap correlations',booterr(*,0)
   ;
   if ierr eq 0 then print,'submap correlation',booterr_lag(0,30)
   ;
                                ; calculate correlation fn at
                                ; different lag for the bootsrapped
                                ; maps
   ind=where(wboottot gt 0,cind)
   rboottot(ind)=rboottot(ind)/wboottot(ind)
   ;
   if ierr eq 0 then begin
      ;
      openw,1,ffdir(nday-1)+field+'.radio'+mode+'_weighted_aftsvd_combined2'
      printf,1,rboottot
      close,1
      ;
      openw,1,ffdir(nday-1)+field+'.weight'+mode+'_aftsvd_combined2'
      printf,1,wboottot
      close,1
      ;
   endif
   ;
   thermal=(dnu*dt)*(wboottot)
   thermal=reform(thermal,nx*ny,nz)
;   gainx=dblarr(nx*ny)
;   gainz=dblarr(nz)
;   for ix=0,nx*ny-1 do begin
;      gainx(ix)=mean(gainxtot(*,*,ix))
;   endfor
;   for iz=0,nz-1 do begin
;      gainz(iz)=mean(gainztot(*,*,iz))
;   endfor
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
   if ierr eq 0 then print,'submaptot correlation',boottot_lag(0,30)
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
;corrlagvalue=corrlagvalue-corrlagerr_mean
;corrlagvalue_boot=corrlagvalue-booterr_mean
;
;
;calculate the noise covariance matrix
for jj=0,nlag-1 do begin
   for kk=0,nlag-1 do begin
      corrcovar(jj,kk)=total(correlation_lag(1:nerror,jj)*correlation_lag(1:nerror,kk))/(double(nerror))
   endfor
endfor

xerr=booterr_err(30)
xerr_mean=booterr_mean(30)
xerr3=boottot_err(30)
xerr3_mean=boottot_mean(30)
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

corrv1=corrold
corrv2=booterr_lag(0,30)
corre1=xerr4
corre2=xerr2
corre3=xerr3
corrm1=xerr4_mean
corrm2=xerr2_mean
;
;print,'redshift-weighted correlation',xcorr0,xcorr0-xerr3_mean,xerr3_mean
;print,'error of correlation with fg-correlation and z-weight',xerr
;print,'error of correlation with bootstrap and z-weight',xerr3
print,'original correlation:',corrold,corrold-xerr2_mean,corrold-xerr4_mean
print,'original correlation:', correlation2(0), correlation_weight2(0),corrlagvalue(30)
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
openw,1,'corrlag.radio.'+field+mode+'.dat' ;,/append
printf,1,corrlagvalue
close,1

openw,1,'correrr.radio.'+field+mode+'.dat' ;,/append
printf,1,corrlagerr
close,1

openw,1,'corrcovar.radio.'+field+mode+'.dat'; ,/append
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

endif


print,'auto-correlation [K]',sqrt(auto)
if nday gt 2 then begin
   err=sqrt(variance(auto)/double(nday-1))
   err=sqrt(err)
   print,'err',err
endif
print,'mean',mean(sqrt(auto))




END

      
      
