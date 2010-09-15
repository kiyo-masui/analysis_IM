pro cmp_plots3,field=field,mode=mode,signal=signal,sigcut=sigcut


if not keyword_set(field) then field='f3'
if not keyword_set(signal) then signal='100muK'
if keyword_set(sigcut) then sigcut='_5sigcut' else sigcut=''

if field eq 'f3' then begin
   ssdir='sept07/f3_27-29/'
   gdir='~/projects/GBT/pros/sept07/f3/radiosim_100muK'
   nx=40
endif
;
if field eq 'f4' then begin
   ssdir='sept06/f4_33-36/'
   gdir='~/projects/GBT/pros/sept06/f4/radiosim_100muK'
   nx=30
endif

;mode=['5mode','10mode','','20mode','25mode']
mode=['1mode','2mode','3mode','4mode','5mode','6mode',$
     '7mode','8mode','9mode','10mode','11mode','12mode',$
     '15mode','20mode','25mode']
;mode=['0mode','1mode','2mode','3mode','4mode','6mode']
;mode=['1mode']
if keyword_set(sigcut) then begin
   sdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'_5sigcut/'
   fdir='/cita/d/raid-cita/tchang/GBT/Processed_5sigcut/'
   fmdir='/cita/d/raid-cita/tchang/GBT/Processed_5sigcut/'
endif else begin
   sdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'/'
   fdir='/cita/d/raid-cita/tchang/GBT/Processed/'
   fmdir='/cita/d/raid-cita/tchang/GBT/Processed/'
endelse
;
set_plot,'ps'
device,/color, bits_per_pixel=8

npol=4
ny=6
nz=370
nnz=560
nmode=n_elements(mode)
radio=dblarr(nmode,nx,ny,nz)
radiosim=radio
weight=dblarr(nmode,nx,ny,nz)
weightsim=weight
model=dblarr(nx,ny,nz)

;model
;openr,1,gdir
openr,1,fmdir+ssdir+field+'.radio.sim.combined_'+signal
readf,1,model
close,1
;
;readf,1,model
;close,1
;model=model(*,*,2+555/3:2+555-1)
help,model

for imode=0,nmode-1 do begin

rtmp=dblarr(nx,ny,nz)
rstmp=rtmp
wtmp=rtmp
wstmp=rtmp

;openr,1,fdir+ssdir+field+'.radio'+mode(imode)+'_weighted_aftsvd_combined2'
;readf,1,rtmp
;close,1

;openr,1,sdir+ssdir+field+'.radio'+mode(imode)+'_weighted_aftsvd_combined2'
;readf,1,rstmp
;close,1

;openr,1,fdir+ssdir+field+'.weight'+mode(imode)+'_aftsvd_combined2'
;readf,1,wtmp
;close,1

;openr,1,sdir+ssdir+field+'.weight'+mode(imode)+'_aftsvd_combined2'
;readf,1,wstmp
;close,1
;
mode1=mode(imode)
if mode1 eq '15mode' then mode1=''
openr,1,fdir+ssdir+field+'.radio.processed.combined'+mode1+'svd'
readf,1,rtmp
close,1
;
openr,1,fdir+ssdir+field+'.weight.processed.combined'+mode1+'svd'
readf,1,wtmp
close,1
;
openr,1,sdir+ssdir+field+'.radio.processed.combined'+mode1+'svd'
readf,1,rstmp
close,1
;
openr,1,sdir+ssdir+field+'.weight.processed.combined'+mode1+'svd'
readf,1,wstmp
close,1
;


radio(imode,*,*,*)=rtmp
radiosim(imode,*,*,*)=rstmp

weight(imode,*,*,*)=wtmp
weightsim(imode,*,*,*)=wstmp

endfor

; calculate weight
newweight=dblarr(nmode,nx,ny,nz)
newweightsim=newweight

; flags
ind=where(weight eq 0,cind)
if cind gt 0 then radio(ind)=0.
ind=where(weightsim eq 0,cind)
if cind gt 0 then radiosim(ind)=0.

for imode=0,nmode-1 do begin
   rtmp=reform(radio(imode,*,*,*))
   rtmp=reform(rtmp,nx*ny,nz)
   wtmp=reform(weight(imode,*,*,*))
   wtmp=reform(wtmp,nx*ny,nz)
   for iz=0,nz-1 do begin
      ind=where(wtmp(*,iz) gt 0,cind)
      rztmp=rtmp(*,iz)
      if cind gt 1 then newweight(imode,*,*,iz)=1./variance(rztmp(ind))
   endfor
endfor

for imode=0,nmode-1 do begin
   rtmp=reform(radiosim(imode,*,*,*))
   rtmp=reform(rtmp,nx*ny,nz)
   wtmp=reform(weightsim(imode,*,*,*))
   wtmp=reform(wtmp,nx*ny,nz)
   for iz=0,nz-1 do begin
      ind=where(wtmp(*,iz) gt 0,cind)
      rztmp=rtmp(*,iz)
      if cind gt 1 then newweightsim(imode,*,*,iz)=1./variance(rztmp(ind))
   endfor
endfor
device,filename=field+'.datamap'+sigcut+'.ps'

for imode=0,nmode-1 do begin
   rtmp=reform(radio(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
;   wtmp2=reform(weight(imode,*,*,*))
;   plt_image,transpose(reform(rtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(rtmp*wtmp2,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(rtmp*wtmp*wtmp2,nx*ny,nz)),/scalable,/colbar
   if imode eq 0 then begin
      openw,1,field+'.datamap'+sigcut+'.dat'
      printf,1,reform(rtmp*wtmp,nx*ny,nz)
      close,1
   endif
endfor

device,/close
;
;
device,filename=field+'.simmap'+sigcut+'.ps'


for imode=0,nmode-1 do begin
   rtmp=reform(radiosim(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
;   wtmp2=reform(weightsim(imode,*,*,*))
;   plt_image,transpose(reform(rtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(rtmp*wtmp2,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(rtmp*wtmp*wtmp2,nx*ny,nz)),/scalable,/colbar
   if imode eq 0 then begin
      openw,1,field+'.simmap'+sigcut+'.dat'
      printf,1,reform(rtmp*wtmp,nx*ny,nz)
      close,1
   endif
endfor

;for imode=0,nmode-1 do begin
;   rtmp=reform(radiosim(imode,*,*,*))
;   wtmp=reform(newweightsim(imode,*,*,*))
;   plt_image,transpose(reform(rtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(wtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar
;endfor

device,/close


diffmap=(radiosim-radio);<(0.0018)>(-0.0005)
ind=where(weight eq 0 or weightsim eq 0,cind)
if cind gt 0 then diffmap(ind)=0.
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   help,rtmp
   openw,1,field+'.diffmap.'+mode(imode)+sigcut+'.dat'
   printf,1,rtmp
   close,1
endfor
;
diffweight=dblarr(nmode,nx,ny,nz)
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   rtmp=reform(rtmp,nx*ny,nz)
   wtmp=reform(weight(imode,*,*,*))
   wtmp=reform(wtmp,nx*ny,nz)
   for iz=0,nz-1 do begin
      ind=where(wtmp(*,iz) gt 0,cind)
      rztmp=rtmp(*,iz)
      if cind gt 1 then diffweight(imode,*,*,iz)=1./variance(rztmp(ind))
   endfor
endfor

device,filename=field+'.diffmap'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   openw,1,field+'.diffmask.'+mode(imode)+sigcut+'.dat'
   printf,1,(rtmp-model)*wtmp
   close,1
   openw,1,field+'.diffnewweight.'+mode(imode)+sigcut+'.dat'
   printf,1,wtmp
   close,1
   openw,1,field+'.diffweight.'+mode(imode)+sigcut+'.dat'
   printf,1,wtmp2
   close,1
;   wtmp2=reform(weight(imode,*,*,*)*weightsim(imode,*,*,*))
;   plt_image,transpose(reform(rtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(model,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform((model-rtmp),nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(rtmp*wtmp*wtmp2,nx*ny,nz)),/scalable,/colbar
;    plt_image,transpose(reform((rtmp*wtmp)<max(model*wtmp)>min(model*wtmp),nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(model*wtmp,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform((rtmp-model)*wtmp,nx*ny,nz)),/scalable,/colbar
;   if imode eq 0 then begin
;      openw,1,field+'.diffmap.dat'
;      printf,1,reform((rtmp-model)*wtmp,nx*ny,nz)
;      close,1
;   endif
endfor
;rtmp=model
;plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar

device,/close

device,filename=field+'.diffmapmasked'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   mask=(rtmp-model)*wtmp
   ind=where(mask ge 100. or mask lt -25.,cind)
;   ind=where(abs(mask) ge 20.,cind)
   print,'masked points:',cind
   mask=(model-rtmp)
   if cind gt 0 then mask(ind)=0.
   openw,1,field+'.diffzeromap.'+mode(imode)+sigcut+'.dat'
   printf,1,mask
   close,1
   mask=(model-rtmp)*wtmp
   if cind gt 0 then mask(ind)=0.
   plt_image,transpose(reform(mask,nx*ny,nz)),/scalable,colbar=[-20,20]
endfor
device,/close

device,filename=field+'.diffmapmasked2'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   mask=(rtmp-model)*wtmp
   ind=where(mask ge 100. or mask lt -25.,cind)
;   ind=where(abs(mask) ge 20.,cind)
   print,'masked points:',cind
   mask=(model-rtmp*1.5)*wtmp
   if cind gt 0 then mask(ind)=0.
   plt_image,transpose(reform(mask,nx*ny,nz)),/scalable,colbar=[-20,20]
endfor
device,/close

device,filename=field+'.diffmapmasked3'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   mask=(rtmp-model)*wtmp
   ind=where(mask ge 100. or mask lt -25.,cind)
;   ind=where(abs(mask) ge 20.,cind)
   print,'masked points:',cind
   mask=(model-rtmp)*wtmp
   if cind gt 0 then mask(ind)=0.
   plt_image,transpose(reform(mask*3.8,nx*ny,nz)),/scalable,colbar=[-20,20]
endfor
device,/close

device,filename=field+'.diffmodelmasked'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   mask=(rtmp-model)*wtmp
   ind=where(mask ge 100. or mask lt -25.,cind)
;   ind=where(abs(mask) ge 20.,cind)
   print,'masked points:',cind
   mask2=model
   if cind gt 0 then mask2(ind)=0.
   openw,1,field+'.modelmap.'+mode(imode)+sigcut+'.dat'
   printf,1,mask2
   close,1
   plt_image,transpose(reform(mask2*wtmp,nx*ny,nz)),/scalable,colbar=[-20,20]
endfor
plt_image,transpose(reform(model*wtmp,nx*ny,nz)),/scalable,colbar=[-20,20]

device,/close


device,filename=field+'.diffmap_noweight'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   mask=(rtmp-model)<(0.002)>(-0.001)
   plt_image,transpose(reform(mask,nx*ny,nz)),/scalable,/colbar
endfor
plt_image,transpose(reform(model,nx*ny,nz)),/scalable,/colbar
device,/close

device,filename=field+'.diffmapsky_noweight'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   mask=(rtmp)<(0.002)>(-0.001)
   plt_image,transpose(reform(mask,nx*ny,nz)),/scalable,/colbar
endfor
plt_image,transpose(reform(model,nx*ny,nz)),/scalable,/colbar
device,/close

device,filename=field+'.diffmapsky'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   mask=(rtmp-model)*wtmp
   ind=where(mask ge 100. or mask lt -25.,cind)
;   mask=(rtmp)<(0.002)>(-0.001)
   mask=rtmp
   if cind gt 0 then mask(ind)=0.
   plt_image,transpose(reform(mask*wtmp,nx*ny,nz)),/scalable,colbar=[-20,20]
endfor
;plt_image,transpose(reform(model,nx*ny,nz)),/scalable,/colbar
device,/close

device,filename=field+'.diffmapsky2'+sigcut+'.ps'
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   wtmp2=reform(weight(imode,*,*,*))
   mask=(rtmp-model)*wtmp
   ind=where(mask ge 100. or mask lt -25.,cind)
;   mask=(rtmp)<(0.002)>(-0.001)
   mask=rtmp
   if cind gt 0 then mask(ind)=0.
   plt_image,transpose(reform(mask*wtmp*1.5,nx*ny,nz)),/scalable,colbar=[-20,20]
endfor
;plt_image,transpose(reform(model,nx*ny,nz)),/scalable,/colbar
device,/close

if keyword_set(corr) then begin

; calculate the correlation
optical=dblarr(nx,ny,nz)

openr,1,fdir+ssdir+field+'.opticalthobar'
readf,1,optical
close,1

;corre=dblarr(nmode)
corresim=corre
corresim2=corre
corrediff=corre

openw,1,field+'corr.comp'+sigcut+'.dat'
for imode=0,nmode-1 do begin

corre(imode)=total(optical*radio(imode,*,*,*)*newweight(imode,*,*,*))/total(newweight(imode,*,*,*))
corresim(imode)=total(optical*radiosim(imode,*,*,*)*newweightsim(imode,*,*,*))/total(newweightsim(imode,*,*,*))
corresim2(imode)=total(optical*radiosim(imode,*,*,*)*newweight(imode,*,*,*))/total(newweight(imode,*,*,*))
corrediff(imode)=total(optical*diffmap(imode,*,*,*)*newweight(imode,*,*,*))/total(newweight(imode,*,*,*))

print,mode(imode)
print,corre(imode),corresim(imode),corresim2(imode),corrediff(imode)
printf,1,mode(imode)
printf,1,corre(imode),corresim(imode),corresim2(imode),corrediff(imode)

endfor
close,1

endif


end
