pro cmp_plots2,field=field,mode=mode


if not keyword_set(field) then field='f3'

if field eq 'f3' then begin
   ssdir='sept07/f3_27-29/'
   gdir='~/projects/GBT/pros/sept07/f3/radiosim_200muK'
   nx=18
endif
;
if field eq 'f4' then begin
   ssdir='sept06/f4_33-36/'
   gdir='~/projects/GBT/pros/sept06/f4/radiosim_200muK'
   nx=12
endif

mode=['5mode','10mode','','20mode','25mode']
sdir='/cita/d/raid-cita/tchang/GBT/Processed_sim/'
fdir='/cita/d/raid-cita/tchang/GBT/Processed/'

set_plot,'ps'
device,/color, bits_per_pixel=8

npol=4
ny=4
nz=370
nnz=560
nmode=n_elements(mode)
radio=dblarr(nmode,nx,ny,nz)
radiosim=radio
weight=dblarr(nmode,nx,ny,nz)
weightsim=weight
model=dblarr(nx,ny,nnz)

;model
openr,1,gdir
readf,1,model
close,1
model=model(*,*,2+555/3:2+555-1)
help,model

for imode=0,nmode-1 do begin

rtmp=dblarr(nx,ny,nz)
rstmp=rtmp
wtmp=rtmp
wstmp=rtmp

openr,1,fdir+ssdir+field+'.radio'+mode(imode)+'_weighted_aftsvd_combined2'
readf,1,rtmp
close,1

openr,1,sdir+ssdir+field+'.radio'+mode(imode)+'_weighted_aftsvd_combined2'
readf,1,rstmp
close,1

openr,1,fdir+ssdir+field+'.weight'+mode(imode)+'_aftsvd_combined2'
readf,1,wtmp
close,1

openr,1,sdir+ssdir+field+'.weight'+mode(imode)+'_aftsvd_combined2'
readf,1,wstmp
close,1

radio(imode,*,*,*)=rtmp
radiosim(imode,*,*,*)=rstmp
weight(imode,*,*,*)=wtmp
weightsim(imode,*,*,*)=wstmp

endfor

; calculate weight
newweight=dblarr(nmode,nx,ny,nz)
newweightsim=newweight

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
device,filename=field+'.datamap.ps'

for imode=0,nmode-1 do begin
   rtmp=reform(radio(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   plt_image,transpose(reform(rtmp,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform(wtmp,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar
endfor

device,/close
;
;
device,filename=field+'.simmap.ps'

for imode=0,nmode-1 do begin
   rtmp=reform(radiosim(imode,*,*,*))
   wtmp=reform(newweightsim(imode,*,*,*))
   plt_image,transpose(reform(rtmp,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform(wtmp,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar
endfor

device,/close

device,filename=field+'.diffmap.ps'

diffmap=(radiosim-radio);<0.00006>(-0.00001)
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
for imode=0,nmode-1 do begin
   rtmp=reform(diffmap(imode,*,*,*))
   wtmp=reform(newweight(imode,*,*,*))
   plt_image,transpose(reform(rtmp,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform(wtmp,nx*ny,nz)),/scalable,/colbar
   plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar
endfor
rtmp=model
plt_image,transpose(reform(rtmp*wtmp,nx*ny,nz)),/scalable,/colbar

device,/close

; calculate the correlation
optical=dblarr(nx,ny,nz)

openr,1,fdir+ssdir+field+'.opticalthobar'
readf,1,optical
close,1

corre=dblarr(nmode)
corresim=corre
corresim2=corre
corrediff=corre

openw,1,field+'corr.comp.dat'
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


end
