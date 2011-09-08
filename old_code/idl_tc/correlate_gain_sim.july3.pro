PRO correlate_gain_sim,corrv,corre,radio,weight,optical,field=field,fin=fin,radiocut=radiocut,fdir=fdir,noise=noise,sim=sim,daisy=daisy,outfile=outfile,dname=dname,svdfg=svdfg,mode=mode
seed=1
if not keyword_set(radiocut) then radiocut=0.1
if not keyword_set(fdir) then fdir=''
if not keyword_set(field) then field='f3'
if not keyword_set(dname) then dname='sept07'
if not keyword_set(mode) then mode=''

if field eq 'f3' then begin
   ;
   if keyword_set(sim) then thisdir='/cita/scratch/cottontail/tchang/GBT/Processed_sim/' else thisdir='/cita/scratch/cottontail/tchang/GBT/Processed/'
   print,'thisdir: ',thisdir
;   thisdir='/Users/tzu/Projects/GBT/pros/newsvd_day/'
   if dname eq 'aug29' then fdir=dname+'/f3_53-58/'
   if dname eq 'aug30' then fdir=dname+'/f3_45-50/'
   if dname eq 'sept07' then fdir=dname+'/f3_27-29/'
;   fdir=['aug29/f3_53-58/','aug30/f3_45-50/','sept07/f3_27-29/';,'sept17/f3_24-29/']
   ffdir=thisdir+fdir
   
   ;
endif

if field eq 'f4' then begin
   ;
   if keyword_set(sim) then thisdir='/cita/scratch/cottontail/tchang/GBT/Processed_sim/' else thisdir='/cita/scratch/cottontail/tchang/GBT/Processed/'
   print,'thisdir: ',thisdir
;   thisdir='/Users/tzu/Projects/GBT/pros/newsvd_day/'
   if dname eq 'aug29' then fdir=dname+'/f4_69-74/'
   if dname eq 'aug30' then fdir=dname+'/f4_67-72/'
   if dname eq 'sept06' then fdir=dname+'/f4_33-36/'
   if dname eq 'sept17' then fdir=dname+'/f4_48-53/'
   ffdir=thisdir+fdir
   ;
endif
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


npol=4
nx=fbin(0)
ny=fbin(1)
;nday=n_elements(ffdir)
;
radio=dblarr(npol,nx,ny,nz)
radiosim=dblarr(nx,ny,nz)
weight=radio
weightsim=radiosim
Da_z,oz,dist
distbin=(floor((dist-d0)/dd)) ;d0 is at bin0
;
; read in data
;for i=0,nday-1 do begin
;   tmp=dblarr(npol,nx,ny,nz)
;openr,1,ffdir+field+'.simradiotot_weighted'
   openr,1,ffdir+field+'.radiotot'+mode+'_weighted_aftsvd'
   readf,1,radio
   close,1
;   radioall(i,*,*,*)=tmp
   ;
;openr,1,ffdir+field+'.simweighttot'
   openr,1,ffdir+field+'.weighttot'+mode+'_aftsvd'
   readf,1,weight
   close,1
;   weightall(i,*,*,*)=tmp
   ;
   openr,1,ffdir+field+'.simradiotot_weighted'
   readf,1,radiosim
   close,1
   ;
   openr,1,ffdir+field+'.simweighttot'
   readf,1,weightsim
   close,1
;endfor
;radio=reform(radio,1,nx,ny,nz)
;weight=reform(weight,1,nx,ny,nz)
;
; combine data, do weighting or not?
; first get rid of the highest dec strip -- let's get rid of this earlier
; simply add up the days
;radio=dblarr(npol,nx,ny,nz)
;weight=radio
;for i=0,nday-1 do begin
;   radio=radio+radioall(i,*,*,*)/double(nday)
;   weight=weight+weightall(i,*,*,*)
;endfor

; apply calibratoin
; reform the arrays
radio=reform(radio,npol,nx*ny,nz)
weight=reform(weight,npol,nx*ny,nz)
;noiseradio=reform(noiseradio,Npol,nx*ny,nz)
help,weight
help,radio
;
device,filename='data.ps'
plt_image,transpose(reform(radio(0,*,*))),/scalable,/colbar
plt_image,transpose(reform(weight(0,*,*))),/scalable,/colbar
;
;
; read in gain_x and gain_z
gainx=dblarr(npol,nx*ny)
gainz=dblarr(npol,nz)
;
gdir='~/projects/GBT/pros/'+dname+'/'+field+'/'
;gdir='/Users/tzu/Projects/GBT/pros/'+dname+'/'
openr,1,gdir+'deep2_gain_x_b4svd.txt'
readf,1,gainx
close,1
;
openr,1,gdir+'deep2_gain_z_b4svd.txt'
readf,1,gainz
close,1
;
;Tsky=dblarr(Npol,nx*ny,nz)
; divide out g(z)g(x) from Psky(x,z)=g(z)g(x)Tsky(x,z) to get Tsky(x,z) 
;for ix=0,nx*ny-1 do begin
;   for iz=0,nz-1 do begin
;      for ipol=0,npol-1 do begin
;         if (gainx(ipol,ix)*gainz(ipol,iz) ne 0.) then begin
;            radio(ipol,ix,iz)=radio(ipol,ix,iz)/gainx(ipol,ix)/gainz(ipol,iz) 
;         endif else begin
;            weight(ipol,ix,iz)=0.
;         endelse
;      endfor
;   endfor
;endfor
;
;
radio=reform(radio,npol,nx,ny,nz)
weight=reform(weight,npol,nx,ny,nz)
help,radio
help,weight

;for field 3, chop off the top dec strip and some noisy strips
;radioall=radio
;weightall=weight
;
if keyword_set(whatever) then begin

if field eq 'f3' then begin
   radioall(*,*,3,*)=0
   weightall(*,*,3,*)=0 
;   radioall(*,*,0,*)=0   ;sept17
;   weightall(*,*,0,*)=0
   radioall(*,14:16,*,*)=0 ;??
   weightall(*,14:16,*,*)=0
;    radioall(*,1:2,0,*)=0  ;sept17
;   weightall(*,1:2,0,*)=0
;   radioall(*,13,*,*)=0   ;for aug29
;   weightall(*,13,*,*)=0
;   radioall(*,2,0,*)=0  ; for aug30
;   weightall(*,2,0,*)=0
endif

endif

; reform the gain too
gainz=gainz(*,2:nz-4)
nzz=555
gainz=gainz(*,nzz/3:nzz-1)
;gainxy=reform(gainx,npol,nx,ny)
;ny=3
;nxx=13
;gainxy=gainxy(*,0:nxx-1,0:ny-1)
;gainx=reform(gainxy,npol,nxx*ny)

; get rid of the empty dec strips
;ny=3
;radio=radio(*,*,0:ny-1,*)
;weight=weight(*,*,0:ny-1,*)
if field eq 'f3' then begin
;
radio(*,*,3,*)=0
weight(*,*,3,*)=0
;
; get rid of the ix=13 for every dec strip
;nx=13
;radio=radio(*,0:nx-1,*,*)
;weight=weight(*,0:nx-1,*,*)
radio(*,13:nx-1,*,*)=0
weight(*,13:nx-1,*,*)=0
;
endif

if field eq 'f4' then begin
;
radio(*,*,0,*)=0
weight(*,*,0,*)=0
;
endif

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

device,filename='radio.ps'

for ipol=0,Npol-1 do begin

;if keyword_set(nosim) then begin
   radio=reform(radioall(ipol,*,*,*))
   weight=reform(weightall(ipol,*,*,*))
   print,'max(radio),min(radio),mean(radio),variance(radio)'
   print,max(radio),min(radio),mean(radio),variance(radio)
   ind=where(weight gt 0,cind)

   ; do a svd on the spatial-redsfhit space
   nz=560
   zstart=2
   radio=radio(*,*,zstart:nz-4)
   weight=weight(*,*,zstart:nz-4)
   nz=555
   nnz=3
    if (field ne 'f1' and field ne '3c286') then begin
        help,radio
        help,weight
        radio=radio(*,*,nz/3:nz-1)
;        ;help,radio
        weight=weight(*,*,nz/3:nz-1)
        if ipol eq 0 then begin
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
        if ipol eq 0 then begin
           inddist=where(distbin ge zstart,count)
           print,'sources with z>zstart',count
           distbin=distbin(inddist)
           ora=ora(inddist)
           odec=odec(inddist)
           distbin=distbin-zstart
        endif
     endelse

     if ipol eq 0 then begin
        radioadd=dblarr(nx*ny,nz)
        weightadd=dblarr(nx*ny,nz)
     endif

     print,'mean HI emission:',mean(radio),min(radio),max(radio)
     ;
;     openw,1,ffdir(nday-1)+field+'.radiotot_weighted_comb_orig'
;     printf,1,radio
;     close,1
     ;
;     openw,1,ffdir(nday-1)+field+'.weighttot_comb_orig'
;     printf,1,weight
;     close,1


;     if keyword_set(sim) then begin
;        radiosim=dblarr(nx,ny,nz)
;        openr,1,ffdir+field+'.radiosim'
;        readf,1,radiosim
;        close,1
;        if n_elements(radio) ne n_elements(radiosim) then begin
;           print,'radio and sim mismatch!'
;        endif
;        radio=radio+radiosim
;     endif

;     radioadd=radioadd+radio*weight
;     weightadd=weightadd+weight
;  endfor



;ind=where(weightadd gt 0,cind)
;radioadd(ind)=radioadd(ind)/weightadd(ind)
;radio=radioadd
;weight=weightadd


if keyword_set(whatever) then begin

; what's these for???

;radio=radio(0:nx/2-1,*,*)
;weight=weight(0:nx/2-1,*,*)
if field eq 'f1' then begin
   nzero=n_elements(weight(1,ny/2:ny-1,*))
   print,'nzero',nzero
   weight(1,ny/2:ny-1,*)=dblarr(nzero)
   weight(*,*,211:219)=0.
   radio(*,*,211:219)=0.
   weight(*,*,260:265)=0.
   radio(*,*,260:265)=0.
endif else begin
   weight(*,*,211-185:222-185)=0.
   radio(*,*,211-185:222-185)=0.
   weight(*,*,260-185:266-185)=0.
   radio(*,*,260-185:266-185)=0.
endelse


radioorig=radio

;; these have been done earlier for the drift scans
radio1d=dblarr(nz)
for i=0,nz-1 do begin
   radioslice=reform(radio(*,*,i))
   weightslice=reform(weight(*,*,i))
   ind=where(weightslice gt 0,cind)
   if cind gt 0 then begin
      ;if cind eq 1 then radio1d(i)=radioslice(ind) else $
        ; radio1d(i)=median(radioslice(ind))
        radio1d(i)=total(radioslice(ind)*weightslice(ind))/total(weightslice(ind))
;       radio1d(i)=mean(radio2d(ind,i))
      radioslice2=dblarr(nx,ny)
      radioslice2(ind)=radioslice(ind)-radio1d(i)
      radio(*,*,i)=radioslice2
   endif
endfor

endif

plt_image,transpose(reform(radio,nx*ny,nz)),/scalable,/colbar
plt_image,transpose(reform(radio*weight,nx*ny,nz)),/scalable,/colbar
plt_image,transpose(reform(weight,nx*ny,nz)),/scalable,/colbar
;

if keyword_set(svdfg) then begin

radiosvd=dblarr(nx*ny,nz)
radiores=radiosvd
radiosvd1mode=dblarr(nx*ny,nz)
radiores1mode=radiosvd
radiosvd2mode=dblarr(nx*ny,nz)
radiores2mode=radiosvd
radio2dtot=reform(radio,nx*ny,nz)
weight2dtot=reform(weight,nx*ny,nz)

plt_image,transpose(weight2dtot),/scalable,/colbar
plt_image,transpose(radio2dtot),/scalable,/colbar

; do weighting
; weight by the mean weight per redshift
totweight=dblarr(nz)
for iz=0,nz-1 do begin
   totweight(iz)=total(weight2dtot(*,iz))
   radio2dtot(*,iz)=radio2dtot(*,iz)*totweight(iz)^2.
endfor

plt_image,transpose(radio2dtot),/scalable,/colbar

if keyword_set(whatever) then begin 
; subtract sturcture on 40 Mpc scales
npix=40./dd
ng=ceil(nz/npix)
npix=40/dd
for jj=0,nx*ny-1 do begin
   ind=where(weight2dtot(jj,*) gt 0,cind)
   if cind gt 0 then begin
      for kk=0,ng-1 do begin
         thisw=weight2dtot(jj,kk*npix:((kk+1)*npix-1)<(nz-1))
         thisr=radio2dtot(jj,kk*npix:((kk+1)*npix-1)<(nz-1))
         ;mask=where(thisw gt 1500,cind)
         mask=where(abs(thisr) gt 5e6,cind)
         if cind gt 0 then begin
            ;med=total(thisr(mask)*thisw(mask))/total(thisw(mask))
            med=mean(thisr(mask))
            thisr(mask)=thisr(mask)-med
            radio2dtot(jj,kk*npix:((kk+1)*npix-1)<(nz-1))=thisr
         endif
      endfor
   endif
endfor

plt_image,transpose(radio2dtot),/scalable,/colbar
endif


la_svd,radio2dtot, s,u,v,/double,status=status
print,'SVD status:',status
   
help,s
help,u
help,v

w=dblarr(n_elements(s))
help,w
w2=w
w(0)=s(0)
;w2(0:4)=s(0:4)  ;sept7
;w2(0:5)=s(0:5)  ;aug29
;w2(0:12)=s(0:12)  ;aug30
;w2(0:2)=s(0:2)  ;sep17
w2(0:0)=s(0:0)

radiosvd=u ## diag_matrix(s) ## transpose(v)
radiosvd1mode=u ## diag_matrix(w) ## transpose(v)
radiosvd2mode=u ## diag_matrix(w2) ## transpose(v)

radiores=radio2dtot-radiosvd
radiores1mode=radio2dtot-radiosvd1mode
radiores2mode=radio2dtot-radiosvd2mode

;plt_image,transpose(radiores),/scalable,/colbar
;plt_image,transpose(radiores1mode),/scalable,/colbar
plt_image,transpose(radiores2mode),/scalable,/colbar


; undo weighting
for iz=0,nz-1 do begin
   if totweight(iz) gt 0 then begin
      radiosvd(*,iz)=radiosvd(*,iz)/totweight(iz)^2.
      radiosvd1mode(*,iz)=radiosvd1mode(*,iz)/totweight(iz)^2.
      radiosvd2mode(*,iz)=radiosvd2mode(*,iz)/totweight(iz)^2.
      radio2dtot(*,iz)=radio2dtot(*,iz)/totweight(iz)^2.
   endif
endfor
;
radiores=radio2dtot-radiosvd
radiores1mode=radio2dtot-radiosvd1mode
radiores2mode=radio2dtot-radiosvd2mode

;plt_image,transpose(radiores),/scalable,/colbar
;plt_image,transpose(radiores1mode),/scalable,/colbar
plt_image,transpose(radiores2mode),/scalable,/colbar

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
print,'mean HI after sigmacut:',mean(radio),mean(radio(where(weight gt 0))),total(radio*weight)/total(weight)
;printf,5,'mean HI after sigmacut:',mean(radio),mean(radio(where(weight gt 0))),total(radio*weight)/total(weight)
print,'rms HI after sigmacut:',sqrt(total(radio^2.*weight^2.)/total(weight^2.))
;printf,5,'rms HI after sigmacut:',sqrt(total(radio^2.*weight^2.)/total(weight^2.))
;

;plt_image,radiores,/scalable,/colbar
plt_image,transpose(radio),/scalable,/colbar
;plt_image,radiores2mode,/scalable,/colbar

radio2modew=radio
for iz=0,nz-1 do begin
   if totweight(iz) gt 0 then begin
      radio2modew(*,iz)=radio2modew(*,iz)*totweight(iz)^2.
   endif
endfor

plt_image,transpose(radio2modew),/scalable,/colbar

endif else begin

radio=reform(radio,nx*ny,nz)
weight=reform(weight,nx*ny,nz)

endelse


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

;read,test

help,radio
help,weight
help,gainx
help,gainz
; set back to kelvin
for ix=0,nx*ny-1 do begin
   for iz=0,nz-1 do begin
;      for ipol=0,npol-1 do begin
         if (gainx(ipol,ix)*gainz(ipol,iz) ne 0.) then begin
            radio(ix,iz)=radio(ix,iz)/gainx(ipol,ix)/gainz(ipol,iz) 
         endif else begin
            weight(ix,iz)=0.
         endelse
;      endfor
   endfor
endfor


radioadd=radioadd+radio*weight
weightadd=weightadd+weight


endfor


; change the dimension of radio and do proper weight
; do weighted average?
ind=where(weightadd gt 0,cind)
radioadd(ind)=radioadd(ind)/weightadd(ind)

;read,test
help,radioadd
help,weightadd
help,gainx
help,gainz
;read,test
; set back to kelvin
;for ix=0,nx*ny-1 do begin
;   for iz=0,nz-1 do begin
;      for ipol=0,npol-1 do begin
;         if (gainx(ipol,ix)*gainz(ipol,iz) ne 0.) then begin
;            radioadd(ipol,ix,iz)=radioadd(ipol,ix,iz)/gainx(ipol,ix)/gainz(ipol,iz) 
;         endif else begin
;            weightadd(ipol,ix,iz)=0.
;         endelse
;      endfor
;   endfor
;endfor


radio=reform(radioadd,nx,ny,nz);   Tsys=40K - change the unit to Kelvin
weight=reform(weightadd,nx,ny,nz)
;radio=reform(radio,nx,ny,nz)
;weight=reform(weight,nx,ny,nz)



;radio=radio*40.d0
;
; write outputs
openw,1,ffdir+field+'.radio.processed.combined'
printf,1,radio
close,1
;
openw,1,ffdir+field+'.weight.processed.combined'
printf,1,weight
close,1

;
plt_image,transpose(reform(radio,nx*ny,nz)),/scalable,/colbar
plt_image,transpose(reform(radio*weight,nx*ny,nz)),/scalable,/colbar
plt_image,transpose(reform(weight,nx*ny,nz)),/scalable,/colbar
;
;read,test
device,/close
;
nerror=1000
correlation=dblarr(nerror+1,nz)
correlation_weight=dblarr(nerror+1,nz)
correlation2=dblarr(nerror+1)
correlation_weight2=correlation2
opticalreal=dblarr(nx,ny,nz)
optmaskreal=dblarr(nx,ny,nz)
lumweight=dblarr(nx,ny,nz)
nlag=61
lagarr=findgen(nlag)-30.
correlation_lag=dblarr(nerror+1,nz,nlag)

for ierr=0,nerror do begin

optical=dblarr(nx,ny,nz)
opticalerr=dblarr(nx,ny,nz)

; use magi as optical weight
magi=magi-max(magi)
lum=10^(-magi/2.5)

if ierr eq 0 then begin
; bin the optical survey
for i=0,n_elements(ora)-1 do begin
   ; figure out spatial bin
   x=(floor( (ora(i)-f0(0))*cos(odec(i)/!radeg) /dra)+nx/2) ;>0<(fbin(0)-1)
   y=(floor( (odec(i)-f0(1))/ddec) +ny/2) ;>0<(fbin(1)-1)
   if ierr eq 0 then zi=distbin(i) else begin
      testi=round(randomu(seed)*double(n_elements(ora)-1))
      zi=distbin(testi)
   endelse
   ;
   ; fill in the nearby cube slice
;   if (x ge 0 and x lt nx and y ge 0 and y lt ny and zi ge 0 and zi lt nz) then
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
               optical(xi,yi,zi)=optical(xi,yi,zi)+exp(-bdist/beam_sigma^2./2.);*lum(i)
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

if keyword_set(kernal) then begin
sigma12=6./dd  ;Mpc
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


;optical=smooth(optical,[1,1,2])

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
;optmask=lumweight
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
    if field eq 'f3' then thobar(i)=mean(optical(*,0:ny-1-1,ind))
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
;print,'mean optical tho:',thobar
;printf,5,'mean optical tho:',thobar
endif

if ierr eq 0 then optmaskreal=optmask
;print,'optmaskreal',max(optmaskreal),min(optmaskreal),mean(optmaskreal)

;read,test
; optical 2d projection
;print,'optical',mean(optical),min(optical),max(optical)
;
if ierr eq 0 then begin
    openw,1,ffdir+field+'.opticaltho'
    printf,1,optical
    close,1
endif
;
;
; subtract off spatial average                                                  
opt1d=dblarr(nz)
for i=0,nz-1 do begin
   ;
   ; ignore ny=3 dec strip here!!
   optslice=reform(optical(*,*,i))
   if field eq 'f3' then optslice=reform(optical(*,0:ny-1-1,i))
;   opt1d(i)=mean(optslice)
   opt1d(i)=total(optslice*optmask2d)/total(optmask2d)
   optical(*,*,i)=optical(*,*,i)-opt1d(i)
;   print,'optical density at:',i, mean(optical(*,*,i))
endfor

if ierr eq 0 then opticalreal=optical

;print,'opt1d:',opt1d

if ierr eq 0 then begin
    openw,1,ffdir+field+'.opticalthobar'
    printf,1,optical
    close,1
    openw,1,ffdir+field+'.opticalthobar_weight'
    printf,1,weight
    close,1
endif

;calculate the error by including the optical correlation
;opticalerr=dblarr(nx,ny,nz)
;opterrmask=dblarr(nx,ny,nz)
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
;   if ii eq noptz-1 then dnoptz=nz-noptz*ii
   radioerr(*,*,ii*dnoptz:((ii+1)*dnoptz-1))=radio(*,*,nstart:(nstart+dnoptz-1))
   weighterr(*,*,ii*dnoptz:((ii+1)*dnoptz-1))=weight(*,*,nstart:(nstart+dnoptz-1))
;   radioerrmask(*,*,ii*dnoptz:((ii+1)*dnoptz-1))=optmaskreal(*,*,nstart:(nstart+dnoptz-1))
endfor

;make a random shuffle of the radio cube
radioran=dblarr(nx,ny,nz)
weightran=dblarr(nx,ny,nz)
;for ix=0,nx-1 do begin
;   for iy=0,ny-1 do begin
      for iz=0,nz-1 do begin
         testi=round(randomu(seed)*double(nz-1))
;         radioran(ix,iy,iz)=radio(ix,iy,testi)
;         weightran(ix,iy,iz)=weight(ix,iy,testi)
;         print,'testi',testi
         radioran(*,*,iz)=radio(*,*,testi)
         weightran(*,*,iz)=weight(*,*,testi)
      endfor
;   endfor
;endfor

;if ierr eq 0 then opticalerr=opticalreal
opticalerr=opticalreal
opterrmask=optmaskreal
;if ierr eq 0 then opterrmask=optmaskreal
; change here for the simulation??
if ierr eq 0 then radioerr=radiosim
if ierr eq 0 then radioran=radiosim
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

; weight by redshift
for ii=0,nz-1 do begin
   ind=where(weighterr(*,*,ii)*opterrmask(*,*,ii) gt 0,cind1)
   ind=where(weightran(*,*,ii)*optmaskreal(*,*,ii) gt 0,cind2)
   if cind1 gt 0 then $
      correlation(ierr,ii)=total(corre(*,*,ii)*weighterr(*,*,ii)*opterrmask(*,*,ii)) $
                        /total(weighterr(*,*,ii)*opterrmask(*,*,ii))
  ; if (ierr ne 0 and ii eq 10) then print,'correlation(ierr,ii)',ierr,ii,correlation(ierr,ii)
   if cind2 gt 0 then $
      correlation_weight(ierr,ii)=total(correboots(*,*,ii)*weightran(*,*,ii)*optmaskreal(*,*,ii)) $
                               /total(weightran(*,*,ii)*optmaskreal(*,*,ii))
  ; if (ierr ne 0 and ii eq 10) then print,'correlation_weight(ierr,ii)',ierr,ii,correlation_weight(ierr,ii)
  ; if (ierr ne 0 and ii eq 10) then read,test
endfor

;print,'correlation_weight',correlation_weight

;read,test
if ierr eq 0 then begin
; calculate the old correlation                                                 
   corrold=total(corre*weight*optmaskreal)/total(weight*optmaskreal)
   print,'original correlation:',corrold
endif

;if ierr eq 0 then begin
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
   for jj=0,nz-1 do begin
      ind=where(weighttmp(*,*,jj)*opterrmask(*,*,jj) gt 0,cind)
      if cind gt 0 then begin
         correlation_lag(ierr,jj,ii)=total(radiotmp(*,*,jj)*opticalerr(*,*,jj)*weighttmp(*,*,jj)*opterrmask(*,*,jj)) $
                                     /total(weighttmp(*,*,jj)*opterrmask(*,*,jj))
      endif
   endfor
endfor


endfor


;print,'correlation_weight',correlation_weight

;read,test


; cooadd the correlation
correrr=dblarr(nz)
correrrslice=dblarr(nz)
for ii=0,nz-1 do begin
   weightslice=weighterr(*,*,ii)
   ind=where(weightslice gt 0,cind)
   if cind gt 0 then begin
      correrrslice(ii)=variance(correlation(1:nerror,ii))
   endif
   weightslice=weightran(*,*,ii)
   ind=where(weightslice gt 0,cind)
   print,'weightslice gt 0',cind
   if cind gt 0 then begin
      correrr(ii)=(variance(correlation_weight(1:nerror,ii)))
      print,'correrr(ii)',ii,correrr(ii)
   endif
endfor
;
;print,'correrr:',correrr
;plot,correrr
;
ind=where(correrr ne 0 and finite(correrr),cind)
print,'correrr:cind',cind
ind3=where(correrrslice ne 0 and finite(correrrslice),cind)
print,'correrrslice:cind',cind
;
;plot,1./correrr(ind)
ind2=where(1./correrr(ind) lt 1e9)
ind4=where(1./correrrslice(ind3) lt 1e9)
;print,'outlier weights:',ind(ind2),1./correrr(ind(ind2))
;
xcorr0=total(correlation_weight(0,ind(ind2))/correrr(ind(ind2)))/total(1./correrr(ind(ind2)))
xcorr=total(correlation(0,ind3(ind4))/correrrslice(ind3(ind4)))/total(1./correrrslice(ind3(ind4)))
corr=dblarr(nerror)
corr2=dblarr(nerror)
corr3=corr
for i=1,nerror do begin
   corr(i-1)=total(correlation(i,ind3(ind4))/correrrslice(ind3(ind4)))/total(1./correrrslice(ind3(ind4)))
   corr2(i-1)=total(correlation_weight(i,ind(ind2))/correrr(ind(ind2)))/total(1./correrr(ind(ind2)))
   ;corr2(i)=mean(correlation_weight(i,ind(ind2)))
endfor
;
;plot,1./correrrslice(ind3(ind4)),/ylog
;plot,corr
;plot,1./correrr(ind(ind2)),/ylog
;plot,corr2
;device,/close
;
;print,'variance(weight)',variance(1./correrrslice(ind3(ind4)))
;print,'variance(weight2)',variance(1./correrr(ind(ind2)))

corrlagvalue=dblarr(nlag)
corrlagerr=dblarr(nlag)
corrlagerr_mean=dblarr(nlag)
corrcovar=dblarr(nlag,nlag)
xcorrlag=dblarr(nerror+1,nlag)
for kk=0,nlag-1 do begin
   correrrlag=dblarr(nz)
;   xcorrlag=dblarr(nerror+1)
   for ii=0,nz-1 do begin
      ind=where(weight(*,*,ii) gt 0,cind)
      if cind gt 0 then begin
         correrrlag(ii)=variance(correlation_lag(1:nerror,ii,kk))
      endif
   endfor
   ind=where(1./correrrlag lt 1e9)
   for jj=0,nerror do begin
      xcorrlag(jj,kk)=total(correlation_lag(jj,ind,kk)/correrrlag(ind))/total(1./correrrlag(ind))
   endfor
   corrlagvalue(kk)=xcorrlag(0,kk)
;   corrlagerr(kk)=sqrt(variance(xcorrlag(1:nerror,kk)))
   corrlagerr(kk)=sqrt(total(xcorrlag(1:nerror,kk)^2.)/double(nerror))
   corrlagerr_mean(kk)=mean(xcorrlag(1:nerror,kk))
   ;
endfor
corrlagvalue=corrlagvalue-corrlagerr_mean
;
;calculate the noise covariance matrix
for jj=0,nlag-1 do begin
   for kk=0,nlag-1 do begin
      corrcovar(jj,kk)=total(xcorrlag(1:nerror,jj)*xcorrlag(1:nerror,kk))/(double(nerror))
   endfor
endfor

xerr=sqrt(variance(corr))
xerr_mean=mean(corr)
xerr2=sqrt(variance(correlation2(1:nerror)))
xerr2_mean=mean(correlation2(1:nerror))
;xerr2=sqrt(variance(corr2))
xerr3=sqrt(variance(corr2))
xerr3_mean=mean(corr2)
xerr4=sqrt(variance(correlation_weight2(1:nerror)))
xerr4_mean=mean(correlation_weight2(1:nerror))
;
;
print,'redshift-weighted correlation with fg-weight',xcorr,xcorr-xerr_mean,xerr_mean

corrv=xcorr-xerr_mean
corre=xerr

print,'redshift-weighted correlation',xcorr0,xcorr0-xerr3_mean,xerr3_mean
print,'error of correlation with fg-correlation and z-weight',xerr
print,'error of correlation with bootstrap and z-weight',xerr3
print,'original correlation:',corrold,corrold-xerr2_mean,corrold-xerr4_mean,xerr2_mean,xerr4_mean
print,'error of correlation with fg-correlation and no z-weight',xerr2
print,'error of correlation with bootstrap and no z-weight',xerr4
print,'rms HI after sigmacut:',sqrt(total(radio^2.*weight^2.)/total(weight^2.))
print,'3-sigma flag:',cind2,double(cind2)/double(nx*ny*double(nz))
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

      
      
