PRO correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,radio,weight,optical,field=field,fin=fin,radiocut=radiocut,fdir=fdir,noise=noise,sim=sim,daisy=daisy,outfile=outfile,dname=dname,svdfg=svdfg,mode=mode,combine=combine,svdmode=svdmode,test=test,mask=mask,dosigcut=dosigcut,sigcut5=sigcut5

seed=1
if not keyword_set(radiocut) then radiocut=0.1
if not keyword_set(fdir) then fdir=''
if not keyword_set(field) then field='f3'
;if not keyword_set(dname) then dname='sept07'
if not keyword_set(mode) then mode=''
if mode eq '15mode' then mode=''
if not keyword_set(signal) then signal='100muK'
if not keyword_set(svdmode) then svdmode=''

if field eq 'f3' then begin
   ;
;   if keyword_set(sim) then
;   thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/'
   if keyword_set(sigcut5) then begin
      if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'_5sigcut/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed_5sigcut/'
   endif else begin
      if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
   endelse
   ;
   print,'thisdir: ',thisdir
   fdir=['/f3_53-58/','/f3_45-50/','/f3_27-29/']
   subdir=['/aug29/f3_05-10/', '/aug29/f3_15-20/',  '/aug29/f3_21-26/',  '/aug29/f3_27-32/',  $
           '/aug29/f3_37-42/',  '/aug29/f3_47-52/',  '/aug29/f3_53-58/', $
           '/aug30/f3_15-20/',  '/aug30/f3_21-26/',  '/aug30/f3_27-32/', $ 
           '/aug30/f3_33-38/',  '/aug30/f3_39-44/',  '/aug30/f3_45-50/', $ 
           '/aug30/f3_9-14/', '/sept07/f3_09-14/', '/sept07/f3_15-20/', $
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
   ;
endif


if field eq 'f4' then begin
   ;
   if keyword_set(sigcut5) then begin
      if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'_5sigcut/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed_5sigcut/'
   endif else begin
      if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
   endelse
   ;if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
;   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/'
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
   ;
endif
set_plot,'ps'
device,/color, bits_per_pixel=8

; cross correlate GBT with the deep2 fields

; GBT parameters
beam=15./60.d0 ; beam size
;dra=beam/2.d
;ddec=beam/2.d
dra=3./60.d0
ddec=3./60.d0
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz
beam_sigma=beam/2.35482


deepdir='/cita/h/home-1/tchang/projects/GBT/Deep2/'
;deepdir='/cita/scratch/cottontail/tchang/GBT/Deep2/'
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
;   f0=[352.49077,0.13736937]
   f0=[352.35077,0.13736937]
   fsize=[2.24112, 0.538381]
   fbin=[40,6]
;   fbin=[18,5]
;   fbin=[18,4]
;   f0=[352.491,0.165758]
;   fbin=round([18,4]*1.5)
   fnopt=deepdir+'deep2_f3.dat'
   readcol,fnopt,ora,odec,oz,magi,format='d,d,d,d'
endif
if (field eq 'f4') then begin
   f0=[37.255199,0.58632926]
   fbin=[30,6]
;   fbin=[12,4]
;   f0=[37.2807,0.615825]
;   fbin=round([12,4]*1.5)
   fnopt=deepdir+'deep2_f4.dat'
   readcol,fnopt,ora,odec,oz,magi,format='d,d,d,d'
endif


npol=4
nx=fbin(0)
ny=fbin(1)
nday=n_elements(ffdir)
;
newnz=370
radioadd=dblarr(nx*ny,newnz)
weightadd=dblarr(nx*ny,newnz)
Da_z,oz,dist
distbin=(floor((dist-d0)/dd)) ;d0 is at bin0
gainxtot=dblarr(nday,npol,nx*ny)
gainztot=dblarr(nday,npol,370)
;
if not keyword_set(test) then begin

; read in data
for iday=0,nday-1 do begin
   ;
   nz=560
   niday=strtrim(string(iday),2)
   radio=dblarr(npol,nx,ny,nz)
   weight=radio
   ;
   openr,1,ffdir(iday)+field+'.radiotot'+mode+'_weighted_aftsvd'
   readf,1,radio
   close,1
   ;
   openr,1,ffdir(iday)+field+'.weighttot'+mode+'_aftsvd'
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
   device,filename='data'+niday+'.ps'
   plt_image,transpose(reform(radio(0,*,*))),/scalable,/colbar
   plt_image,transpose(reform(weight(0,*,*))),/scalable,/colbar
   ;
   ;
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
   ;
   radio=reform(radio,npol,nx,ny,nz)
   weight=reform(weight,npol,nx,ny,nz)
   ;
   ;
   ; reform the gain too
   gainz=gainz(*,2:nz-4)
   nzz=555
   gainz=gainz(*,nzz/3:nzz-1)
   ;
   gainxtot(iday,*,*)=gainx
   gainztot(iday,*,*)=gainz
   ;
;   if field eq 'f3' then begin
      ;
;      radio(*,*,3,*)=0
;      weight(*,*,3,*)=0
      ;
; get rid of the ix=13 for every dec strip
;nx=13
;      radio(*,13:nx-1,*,*)=0
;      weight(*,13:nx-1,*,*)=0
;      radio(*,0:3,*,*)=0
;      weight(*,0:3,*,*)=0
      ;
;   endif
   ;
   ;
;   if field eq 'f4' then begin
      ;
;      radio(*,*,0,*)=0
;      weight(*,*,0,*)=0
      ;
;   endif
   ;
   help,radio
   help,weight
   ;
   radiotmp=reform(radio,npol,nx*ny,nz)
   weightmp=reform(weight,npol,nx*ny,nz)
   plt_image,transpose(reform(radiotmp(0,*,*))),/scalable,/colbar
   plt_image,transpose(reform(weightmp(0,*,*))),/scalable,/colbar
   ;
   device,/close
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
   device,filename='radio'+niday+'.ps'
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
      nz=560
      zstart=2
      radio=radio(*,*,zstart:nz-4)
      weight=weight(*,*,zstart:nz-4)
      neweight=neweight(*,*,zstart:nz-4)
      nz=555
      nnz=3
      if (field ne 'f1' and field ne '3c286') then begin
         help,radio
         help,weight
         radio=radio(*,*,nz/3:nz-1)
;        ;help,radio
         weight=weight(*,*,nz/3:nz-1)
         neweight=neweight(*,*,nz/3:nz-1)
         if iday eq 0 and ipol eq 0 then begin
            inddist=where(distbin ge (zstart+nz/3),count)
            print,'sources with z>0.7',count
            distbin=distbin(inddist)
            ora=ora(inddist)
            odec=odec(inddist)
            distbin=distbin-nz/3-zstart
         endif
         nz=(nz/3)*2
         nnz=2
      endif else begin
         if iday eq 0 and ipol eq 0 then begin
            inddist=where(distbin ge zstart,count)
            print,'sources with z>zstart',count
            distbin=distbin(inddist)
            ora=ora(inddist)
            odec=odec(inddist)
            distbin=distbin-zstart
         endif
      endelse
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
         plt_image,transpose(radio2dtot),/scalable,/colbar
         ;
         radio2dtot=radio2dtot*neweight2dtot
         plt_image,transpose(radio2dtot),/scalable,/colbar
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
         plt_image,transpose(radiores1mode),/scalable,/colbar
         plt_image,transpose(radiores2mode),/scalable,/colbar
         ;
         ;undo weighting
         radiores1mode=radiores1mode/neweight2dtot
         radiores2mode=radiores2mode/neweight2dtot
         radio2dtot=radio2dtot/neweight2dtot
         ;
         ; do another 3-sigma cut
         sigcut=3.
;         if mode eq '0mode' then sigcut=8.
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
         plt_image,transpose(radio),/scalable,/colbar
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
;      if mode eq '0mode' then sigcut=8.
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
      ;
   endfor
   ;
endfor
;
ind=where(weightadd gt 0,cind)
radioadd(ind)=radioadd(ind)/weightadd(ind)
; 
; calculate the weight(x), sum over redshifts
weightx=dblarr(nx*ny,nz)
for ix=0,nx*ny-1 do begin
   weightx(ix,*)=total(weightadd(ix,*))
endfor
;
;
help,radioadd
help,weightadd
help,gainx
help,gainz
;
radio=reform(radioadd,nx,ny,nz);   Tsys=40K - change the unit to Kelvin
weight=reform(weightadd,nx,ny,nz)
weightx=reform(weightx,nx,ny,nz)
;
;
; write outputs
openw,1,ffdir(nday-1)+field+'.radio.processed.combined'+mode+'svd'+svdmode
printf,1,radio
close,1
;
openw,1,ffdir(nday-1)+field+'.weight.processed.combined'+mode+'svd'+svdmode
printf,1,weight
close,1
;
if keyword_set(combine) then begin
   ;
   radio=dblarr(nx,ny,nz)
   weight=dblarr(nx,ny,nz)
   ;
   openr,1,ffdir(nday-1)+field+'.radio'+mode+'_weighted_aftsvd_combined2'
   readf,1,radio
   close,1
   ;
   openr,1,ffdir(nday-1)+field+'.weight'+mode+'_aftsvd_combined2'
   readf,1,weight
   close,1
   ;
   print,'Hola2!'
   read,test
   ;
endif
;
;
; use last day's gain to set thermal threshould to kelvin
thermal=(dnu*dt)*(weight)
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
neweight=dblarr(nx*ny,nz)
weight2d=reform(weight,nx*ny,nz)
radio2d=reform(radio,nx*ny,nz)
for iz=0,nz-1 do begin
   ind=where(weight2d(*,iz) gt 0,cind)
   if cind gt 1 then begin
      radiotmp=radio2d(*,iz)
      neweight(*,iz)=(1./variance(radiotmp(ind)));<thermal(*,iz)
   endif
endfor
;
;
;plt_image,transpose(radio2d),/scalable,/colbar
;plt_image,transpose(radio2d*weight2d),/scalable,/colbar
;plt_image,transpose(radio2d*neweight),/scalable,/colbar
;
device,/close
;
device,filename=field+'.skymap'+mode+'.ps'
plt_image,transpose(radio2d),/scalable,/colbar
plt_image,transpose(radio2d*weight2d),/scalable,/colbar
plt_image,transpose(radio2d*neweight),/scalable,/colbar
plt_image,transpose(radio2d*neweight*reform(weightx,nx*ny,nz)),/scalable,/colbar
device,/close
;
neweight=reform(neweight,nx,ny,nz)

endif else begin
   
   nz=370
   radio=dblarr(nx,ny,nz)
   weight=radio
   neweight=radio
   help,radio
   ;
   openr,1,'~/projects/GBT/pros/new_weight/'+field+'.modelmap.'+mode+'.dat'
   readf,1,radio
   close,1
   ;
   openr,1,'~/projects/GBT/pros/new_weight/'+field+'.diffnewweight.'+mode+'.dat'
   readf,1,neweight
   close,1
   ;
   openr,1,'~/projects/GBT/pros/new_weight/'+field+'.diffweight.'+mode+'.dat'
   readf,1,weight
   close,1
   ;
   if keyword_set(mask) then begin
      mask=dblarr(nx,ny,nz)
      mask2=make_array(nx,ny,nz,value=1.)
      openr,1,'~/projects/GBT/pros/new_weight/'+field+'.diffmask.'+mode+'.dat'
      readf,1,mask
      close,1
      ;
;      ind=where(abs(mask) ge 20.,cind)
      ind=where(mask ge 100. or mask lt -25.,cind)
      print,'masked points:',cind
      if cind gt 0 then mask2(ind)=0.
      radio=radio*mask2
   endif
   ; calculate the weight(x), sum over redshifts
   weightx=dblarr(nx*ny,nz)
   weightadd=reform(weight,nx*ny,nz)
   for ix=0,nx*ny-1 do begin
      weightx(ix,*)=total(weightadd(ix,*))
   endfor
   weightx=reform(weightx,nx,ny,nz)

endelse


weightorig=weight
weight=neweight
;
;
if keyword_set(bootstrap) then begin

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
for imap=0,nsubmap/npol-1 do begin
   ;
   nzz=560
   rtmp=dblarr(npol,nx,ny,nzz)
   wtmp=rtmp
   nz=370
   rtmp2=dblarr(npol,nx,ny,nz)
   wtmp2=rtmp
   openr,1,submapdir(imap)+field+'.radio'+mode+'_weighted_aftsvd'
   readf,1,rtmp
   close,1
      ;
   openr,1,submapdir(imap)+field+'.weight'+mode+'_aftsvd'
   readf,1,wtmp
   close,1
      ;
   rtmp2=rtmp(*,*,*,(2+555/3):(2+555-1))
   wtmp2=wtmp(*,*,*,(2+555/3):(2+555-1))
      ;
;   if field eq 'f3' then begin
      ;
;      rtmp2(*,*,3,*)=0
;      wtmp2(*,*,3,*)=0
      ;
;      rtmp2(*,13:nx-1,*,*)=0
;      wtmp2(*,13:nx-1,*,*)=0
      ;
;   endif
;   if field eq 'f4' then begin
      ;
;      rtmp2(*,*,0,*)=0
;      wtmp2(*,*,0,*)=0
      ;
;   endif
;   rtmp2(*,0:3,*,*)=0
;   wtmp2(*,0:3,*,*)=0
   ;
   if keyword_set(dosigcut) then begin
   ; do a sigmacut
   sigcut=3.
   nfreq=2048
   dnu=50e6/double(nfreq)
   dt=0.5d
   sigmal=1./sqrt(dnu*dt)
   ;
   ind=where(wtmp2 gt 0,cind)
   sigmacut=sigcut*sigmal/sqrt(wtmp2(ind))
   ;
   ind2=where(abs(rtmp2(ind)) gt sigcut*sigmal/sqrt(wtmp2(ind)),cind2)
   print,'nx,ny,nz',nx,ny,nz,double(nx*ny*double(nz))
   print,'3-sigma flag:',cind2,double(cind2)/double(nx*ny*double(nz)*npol)
      ;
   if cind2 gt 0 then begin
      rtmp2(ind(ind2))=0.
      wtmp2(ind(ind2))=0
   endif
      ;
endif
   ;
   wtmp2=reform(wtmp2,npol,nx*ny,nz)
   rtmp2=reform(rtmp2,npol,nx*ny,nz)
   ;
   for ipol=0,npol-1 do begin
      thermal=dblarr(nx*ny,nz)
      for ix=0,nx*ny-1 do begin
         for iz=0,nz-1 do begin
            if (gainxtot(dayindx(imap),ipol,ix)*gainztot(dayindx(imap),ipol,iz) ne 0.) then begin
               thermal(ix,iz)=(dnu*dt)*(wtmp2(ipol,ix,iz))/gainxtot(dayindx(imap),ipol,ix)/gainztot(dayindx(imap),ipol,iz) 
               rtmp2(ipol,ix,iz)=rtmp2(ipol,ix,iz)/gainxtot(dayindx(imap),ipol,ix)/gainztot(dayindx(imap),ipol,iz) 
            endif else begin
               wtmp2(ipol,ix,iz)=0.
            endelse
         endfor
      endfor
      thermal=reform(thermal,nx,ny,nz)
      for iz=0,nz-1 do begin
         ind=where(wtmp2(ipol,*,iz) gt 0,cind)
         if cind gt 1 then begin
            radiotmp=rtmp2(ipol,*,iz)
            subwmap(imap*npol+ipol,*,*,iz)=(1./variance(radiotmp(ind)));<thermal(*,*,iz)
         endif
      endfor
   endfor
   ;
   rtmp2=reform(rtmp2,npol,nx,ny,nz)
   subrmap((imap*npol):(npol*(imap+1)-1),*,*,*)=rtmp2
;   wtmp2=reform(submap,npol,nx,ny,nz)
;   subweightmap((imap*npol):(npol*(imap+1)-1),*,*,*)=wtmp2
   wtmp2=reform(wtmp2,npol,nx,ny,nz)
   subweightmap((imap*npol):(npol*(imap+1)-1),*,*,*)=wtmp2
   ;
endfor
;
endif



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

optical=dblarr(nx,ny,nz)
opticalerr=dblarr(nx,ny,nz)

; read in all sub maps to calculate the errors
; 
; make the optical map

for ierr=0,nerror do begin

if keyword_set(mk_opt) then begin
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
            for ii=-3, 3 do begin
               for jj=-3, 3 do begin
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
;         if field eq 'f3' then thobar(i)=mean(optical(*,0:ny-1-1,ind))
;         if field eq 'f4' then thobar(i)=mean(optical(*,1:ny-1,ind))
         if field eq 'f3' then thobar(i)=mean(optical(*,*,ind))
         if field eq 'f4' then thobar(i)=mean(optical(*,*,ind))
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
;      if field eq 'f3' then optslice=reform(optical(*,0:ny-1-1,i))
;      if field eq 'f4' then optslice=reform(optical(*,1:ny-1,i))
      optslice=reform(optical(*,*,i))
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

endif else begin
   
   openr,1,'~/projects/GBT/pros/optical/'+field+'.opticaldensity.dat'
   readf,1,optical
   close,1

endelse


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
opticalreal=optical
optmaskreal=optmask

radioerr=dblarr(nx,ny,nz)
weighterr=dblarr(nx,ny,nz)
weightorigerr=dblarr(nx,ny,nz)

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
      weightorigerr(*,*,ii*dnoptz:((ii+1)*dnoptz-1))=weightorig(*,*,nstart:(nstart+dnoptz-1))
   endfor

   ;make a random shuffle of the radio cube
   radioran=dblarr(nx,ny,nz)
   weightran=dblarr(nx,ny,nz)
   weightorigran=dblarr(nx,ny,nz)
   for iz=0,nz-1 do begin
      testi=round(randomu(seed)*double(nz-1))
      radioran(*,*,iz)=radio(*,*,testi)
      weightran(*,*,iz)=weight(*,*,testi)
      weightorigran(*,*,iz)=weightorig(*,*,testi)
   endfor

   ;if ierr eq 0 then opticalerr=opticalreal
   opticalerr=opticalreal
   opterrmask=optmaskreal
   if ierr eq 0 then radioerr=radio
   if ierr eq 0 then radioran=radio
   if ierr eq 0 then weighterr=weight
   if ierr eq 0 then weightran=weight
   if ierr eq 0 then weightorigerr=weightorig
   if ierr eq 0 then weightorigran=weightorig
 
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
   correlation2(ierr)=total(corre*weighterr*opterrmask*weightx)/total(weighterr*opterrmask*weightx)

   ; bootstrap, all z
   correlation_weight2(ierr)=total(correboots*weightran*optmask*weightx)/total(weightran*optmask*weightx)
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
      corrold=total(corre*weight*optmaskreal*weightx)/total(weight*optmaskreal*weightx)
      print,'original correlation:',corrold
   endif

   ; calculate correlation fn at different lag
   for ii=0,nlag-1 do begin
      radiotmp=dblarr(nx,ny,nz)
      weighttmp=dblarr(nx,ny,nz)
      weightorigtmp=dblarr(nx,ny,nz)
      lag=long(lagarr(ii))
      if lag ge 0 then begin
         radiotmp(*,*,0:nz-1-lag)=radioran(*,*,lag:nz-1)
         weighttmp(*,*,0:nz-1-lag)=weightran(*,*,lag:nz-1)
         weightorigtmp(*,*,0:nz-1-lag)=weightorig(*,*,lag:nz-1)
      endif
      if lag lt 0 then begin
         radiotmp(*,*,-lag:nz-1)=radioran(*,*,0:nz-1+lag)
         weighttmp(*,*,-lag:nz-1)=weightran(*,*,0:nz-1+lag)
         weightorigtmp(*,*,-lag:nz-1)=weightorig(*,*,0:nz-1+lag)
      endif
;      for jj=0,nz-1 do begin
      ind=where(weighttmp*opterrmask gt 0,cind)
      if cind gt 0 then begin
         correlation_lag(ierr,ii)=total(radiotmp*opticalreal*weighttmp*opterrmask*weightx) $
                                  /total(weighttmp*opterrmask*weightx)
      endif
   endfor
   ;
   ;
if keyword_set(bootstrap) then begin

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
         booterr(imap,ii)=total(radiotmp*opticalerr*weighttmp*opterrmask*weightx) $
                          /total(weighttmp*opterrmask*weightx)
      endfor
   endfor
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
         wboottot(*,*,iz)=(1./variance(radiotmp(ind)));<thermal(*,*,iz)
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
         boottot_lag(ierr,ii)=total(radiotmp*opticalerr*weighttmp*opterrmask*weightx) $
                                  /total(weighttmp*opterrmask*weightx)
      endif
   endfor
   ;
   if ierr eq 0 then print,'submaptot correlation',boottot_lag(0,30)

endif

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
openw,1,'corrlag.radio.'+field+mode+'svd'+svdmode+'.dat' ;,/append
printf,1,corrlagvalue
close,1

openw,1,'correrr.radio.'+field+mode+'svd'+svdmode+'.dat' ;,/append
printf,1,corrlagerr
close,1

openw,1,'corrcovar.radio.'+field+mode+'svd'+svdmode+'.dat'; ,/append
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

      
      
