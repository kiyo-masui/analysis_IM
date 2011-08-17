pro mk_optmap,field=field

if not keyword_set(field) then field='f3'

set_plot,'ps'
device,/color, bits_per_pixel=8


; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
;dra=3./60.d0
;ddec=3./60.d0
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz
beam_sigma=beam/2.35482

cosmosdir='/cita/h/home-1/tchang/projects/GBT/zCOSMOS/'
f0=[150.11635,2.2247776]
fbin=[9,9]
readcol,cosmosdir+'zCOSMOS.ra',ora,format='d'
readcol,cosmosdir+'zCOSMOS.dec',odec,format='d'
readcol,cosmosdir+'zCOSMOS.z',oz,format='d'

print,'ra range',min(ora),max(ora)
print,'dec range',min(odec),max(odec)
print,'z range',min(oz),max(oz)
help,oz

nx=fbin(0)
ny=fbin(1)
nz=560
Da_z,oz,dist
distbin=(floor((dist-d0)/dd)) ;d0 is at bin0


nx=nx*2
ny1=ny*2
optical=dblarr(nx,ny1,nz)
opt2d=dblarr(nx,ny1,nz)
maxz=0.
minz=10.
for i=0L,n_elements(ora)-1 do begin
   ; figure out spatial bin
   x=(floor( (ora(i)-f0(0))*cos(odec(i)/!radeg) /dra)+nx/2) ;>0<(fbin(0)-1)
   y=(floor( (odec(i)-f0(1))/ddec) +ny1/2)                   ;>0<(fbin(1)-1)
   zi=distbin(i) 
   ;
   ; fill in the nearby cube slice
   if (x ge 0 and x lt nx and y ge 0 and y lt ny1 and zi ge 0 and zi lt nz) then begin
      optical(x,y,zi)=optical(x,y,zi)+1. ;*lum(i)
      opt2d(x,y)=opt2d(x,y)+1.
      if oz(i) gt maxz then maxz=oz(i)
      if oz(i) lt minz then minz=oz(i)
   endif
endfor
ind=where(optical eq 0, cind)
print,'empty bins',cind, double(cind)/double(nx*double(ny1)*nz)
print,'total number of galaxies',total(optical)
print,'redshift range',minz,maxz
;
; cut out the edges
mask=make_array(nx,ny1,nz,value=1.)
mask2d=make_array(nx,ny1,value=1.)
;
;opt2d=reform(optical(*,*,nz/2))
for x=0,nx-1 do begin
   for y=0,ny1-1 do begin
      if (opt2d(x,y) lt 1) then begin
         mask(x,y,*)=0
         mask2d(x,y)=0
      endif
   endfor
endfor
;
print,'optical',max(optical),min(optical),mean(optical)
ind=where(mask eq 0,cind)
print,'zero mask',cind
; 
; calculate delta
opticaltho=optical
noptz=4
thobar=dblarr(noptz)
doptz=nz/noptz
index=indgen(nz)
for i=0,noptz-1 do begin
   ind=where((index/doptz) eq i,cind)
   if cind ne doptz then print,'problem in optical rebinning!'
         ;
         ; ignore ny=3 dec strips here!!
   thobar(i)=total(optical(*,*,ind)*mask(*,*,ind))/total(mask(*,*,ind))
   opticaltho(*,*,ind)=opticaltho(*,*,ind)/thobar(i)-1.
endfor
;
opticaltho=opticaltho*mask
ind=where(opticaltho eq 0.,cind)
ind2=where(mask eq 0, cind2)
print,'zero optical bins:',cind
print,'zero mask:',cind2
print,'opticaltho',max(opticaltho),min(opticaltho),mean(opticaltho)
device,filename='zCOSMOS.opticaltho.ps'
help,opticaltho
print,nx,ny1,nz
optplot=reform(opticaltho,nx*ny1,nz)
plt_image,transpose(optplot),/scalable,/colbar
plt_image,opt2d,/scalable,/colbar
plt_image,mask2d,/scalable,/colbar
plt_image,transpose(reform(mask,nx*ny1,nz)),/scalable,/colbar
plt_image,transpose(reform(mask(*,ny1/2,*))),/scalable,/colbar
plt_image,transpose(reform(mask(nx/2,*,*))),/scalable,/colbar
plt_image,transpose(reform(opticaltho(nx/2,*,*))),/scalable,/colbar
;optmask=dblarr(nx,ny1,nz)
;optmask2d=dblarr(nx,ny1)
;for i=0,nx-1 do begin
;   for j=0,ny1-1 do begin
;      med=total(reform(opticaltho(i,j,*)))
;      if med ne 0 then med=1
;      optmask2d(i,j)=med
;      for k=0,nz-1 do begin
;         optmask(i,j,k)=med
;      endfor
;   endfor
;endfor
;
opt1d=dblarr(nz)
for i=0,nz-1 do begin
   ;
   optslice=reform(opticaltho(*,*,i))
   opt1d(i)=total(optslice*reform(mask(*,*,i)))/total(mask(*,*,i))
   opticaltho(*,*,i)=opticaltho(*,*,i)-opt1d(i)
endfor
;
opticaltho=opticaltho*mask
optplot=reform(opticaltho,nx*ny1,nz)
plt_image,transpose(optplot),/scalable,/colbar
device,/close

print,'thobar',thobar
print,min(opticaltho),max(opticaltho),mean(opticaltho)
; now convloved with radio beam
nx0=fbin(0)
density=dblarr(nx0,ny,nz)
dny=(ny1-ny)/2
dnx=(nx-nx0)/2
for x=0,nx0-1 do begin      
   for y=0,ny-1 do begin
      xcor=indgen(nx)-x-dnx
      ycor=indgen(ny1)-y-dny
      map=dblarr(nx,ny1,nz)
      for xi=0,nx-1 do begin
         for yi=0,ny1-1 do begin
            map(xi,yi,*)=exp(-((xcor(xi)*dra)^2.+ (ycor(yi)*ddec)^2.)/beam_sigma^2/2.)
         endfor
      endfor
      for z=0,nz-1 do begin
         density(x,y,z)=total(opticaltho(*,*,z)*map(*,*,z))/total(map(*,*,z));/(2*!dpi*beam_sigma^2.)
      endfor
   endfor
endfor

; write output
openw,1,'~/projects/GBT/zCOSMOS/zCOSMOS.opticaldensity.dat'
printf,1,density
close,1
;
help,density
;make a plot
device,filename='zCOSMOS.density_convolved.ps'
device,xsize=18.,ysize=6.
den2d=reform(density,nx0*ny,nz)
help,density
print,max(density),min(density),mean(density)
;radio1=reform(den2d(0:nx0-1,*))
;radio2=reform(den2d((ny-1)*nx0:nx0*ny-1,*))
;radio=radio1
;radio=[radio1,radio2]
;var=sqrt(variance(radio))
;print,'var:',var
;print,'mean',mean(radio)
;var=sqrt(total(radio^2.)/double(nx0*nz*2.))
;print,'var:',var

;plt_image,transpose(radio),/scalable,/colbar,frame=[1800,2600,0,120],xtitle='!6 Redshift Distance [Mpc]',ytitle='Spatial Distance [Mpc]'
plt_image,transpose(den2d),/scalable,/colbar
plt_image,transpose(reform(density(*,ny/2,*))),/scalable,/colbar
plt_image,transpose(reform(density(nx0/2,*,*))),/scalable,/colbar
device,/close


end
