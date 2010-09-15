pro plt_maps,field=field


if not keyword_set(field) then field='f3'
if not keyword_set(signal) then signal='100muK'
mode='1mode'

set_plot,'ps'
device,/color, bits_per_pixel=8

;nx=40
nx=30
ny=6
nz=370
npol=4

;device,filename='field3.maps.ps'

radio=dblarr(nx,ny,nz)
weight=radio
neweight=radio

openr,1,'~/projects/GBT/pros/new_weight/'+field+'.diffnewweight.'+mode+'.dat'
readf,1,neweight
close,1

openr,1,'~/projects/GBT/pros/new_weight/'+field+'.diffweight.'+mode+'.dat'
readf,1,weight
close,1

weightx=dblarr(nx*ny,nz)
weightadd=reform(weight,nx*ny,nz)
for ix=0,nx*ny-1 do begin
   weightx(ix,*)=total(weightadd(ix,*))
endfor
;
;
weight=neweight
weight=reform(weight,nx*ny,nz)
;
;
; read in data
radiopol1=dblarr(npol,nx,ny,560)
radiopol2=dblarr(npol,nx,ny,560)
radiopol3=dblarr(npol,nx,ny,560)
; 0mode
;openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_mode_simonly_100muK/aug29/f3_05-10/f3.radio0mode_weighted_aftsvd'
openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_mode_simonly_100muK/sept06/f4_33-36/f4.radiotot0mode_weighted_aftsvd'
readf,1,radiopol3
close,1
;
;openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_mode/aug29/f3_05-10/f3.radio0mode_weighted_aftsvd'
openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_5sigcut/sept06/f3_33-36/f4.radio0mode_weighted_aftsvd'
readf,1,radiopol2
close,1
;
openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_mode_sim_100muK/aug29/f3_05-10/f3.radio0mode_weighted_aftsvd'
readf,1,radiopol1
close,1
;
radiopol1=radiopol1(*,*,*,2+185:555+2-1)
radio1=reform(radiopol1(0,*,*,*))
radiopol2=radiopol2(*,*,*,2+185:555+2-1)
radio2=reform(radiopol2(0,*,*,*))
radiopol3=radiopol3(*,*,*,2+185:555+2-1)
radio3=reform(radiopol3(0,*,*,*))
;
; 1mode
radiopol1=dblarr(npol,nx,ny,560)
radiopol2=dblarr(npol,nx,ny,560)
radiopol3=dblarr(npol,nx,ny,560)
radiopol4=dblarr(npol,nx,ny,560)
; 
openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_mode/aug29/f3_05-10/f3.svdresidual_weighted_'+mode
readf,1,radiopol3
close,1
;
openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_mode/aug29/f3_05-10/f3.radio'+mode+'_weighted_aftsvd'
readf,1,radiopol2
close,1
;
openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_mode_sim_100muK/aug29/f3_05-10/f3.radio'+mode+'_weighted_aftsvd'
readf,1,radiopol1
close,1
;
openr,1,'/cita/d/raid-cita/tchang/GBT/Processed_mode/aug29/f3_05-10/f3.svddiff_weighted_'+mode
readf,1,radiopol4
close,1
;
radiopol1=radiopol1(*,*,*,2+185:555+2-1)
radio4=reform(radiopol1(0,*,*,*))
radiopol2=radiopol2(*,*,*,2+185:555+2-1)
radio5=reform(radiopol2(0,*,*,*))
radiopol3=radiopol3(*,*,*,2+185:555+2-1)
radio6=reform(radiopol3(0,*,*,*))
radiopol4=radiopol4(*,*,*,2+185:555+2-1)
radio7=reform(radiopol4(0,*,*,*))
;
gainx=dblarr(npol,nx*ny)
gainz=dblarr(npol,nz)
   ;
gdir='~/projects/GBT/pros/aug29/'+field+'/'
openr,1,gdir+'deep2_gain_x_b4svd.txt'
readf,1,gainx
close,1
   ;
openr,1,gdir+'deep2_gain_z_b4svd.txt'
readf,1,gainz
close,1
   ;
gainx=reform(gainx(0,*))
gainz=reform(gainz(0,*))
   ;
radio1=reform(radio1,nx*ny,nz)
radio2=reform(radio2,nx*ny,nz)
radio3=reform(radio3,nx*ny,nz)
radio4=reform(radio4,nx*ny,nz)
radio5=reform(radio5,nx*ny,nz)
radio6=reform(radio6,nx*ny,nz)
radio7=reform(radio7,nx*ny,nz)
for ix=0,nx*ny-1 do begin
   for iz=0,nz-1 do begin
      if (gainx(ix)*gainz(iz) ne 0.) then begin
         radio1(ix,iz)=radio1(ix,iz)/gainx(ix)/gainz(iz) 
         radio2(ix,iz)=radio2(ix,iz)/gainx(ix)/gainz(iz) 
         radio3(ix,iz)=radio3(ix,iz)/gainx(ix)/gainz(iz) 
         radio4(ix,iz)=radio4(ix,iz)/gainx(ix)/gainz(iz) 
         radio5(ix,iz)=radio5(ix,iz)/gainx(ix)/gainz(iz) 
         radio6(ix,iz)=radio6(ix,iz)/gainx(ix)/gainz(iz) 
         radio7(ix,iz)=radio7(ix,iz)/gainx(ix)/gainz(iz) 
      endif
   endfor
endfor
;
;
!p.multi=[0,1,3]
device,filename=field+'.maps.0mode.ps'
;plt_image,transpose((radio1-radio2)*weight*weightx),/scalable,colbar=[-2e7,6e7]
;plt_image,transpose(radio3),/colbar,/scalable,$
;          frame=[1800,2600,0,360],xtitle='!6 Redshift Distance [Mpc]',$
;          ytitle='Spatial Distance [Mpc]'
plt_image,transpose(radio3*weight*weightx),/colbar,/scalable,$
          frame=[1800,2600,0,360],xtitle='!6 Redshift Distance [Mpc]',$
          ytitle='Spatial Distance [Mpc]'
;plt_image,transpose((radio1-radio2-radio3)*weight*weightx),/scalable,colbar=[-2e7,6e7]
;
device,/close
!p.multi=0

device,filename=field+'.maps.10mode.ps'
;plt_image,transpose((radio4-radio5)*weight*weightx),/scalable,colbar=[-2e7,6e7]
plt_image,transpose(radio6*weight*weightx),/scalable,colbar=[-2e7,6e7]
;plt_image,transpose((radio4-radio5-radio6)*weight*weightx),/scalable,colbar=[-2e7,6e7]
;
device,/close

device,filename=field+'.maps.1modediff.ps'
plt_image,transpose((radio1-radio4-radio2+radio5)*weight*weightx),/scalable,colbar=[-2e7,6e7]
plt_image,transpose((radio3-radio6)*weight*weightx),/scalable,colbar=[-2e7,6e7]
plt_image,transpose((radio7)*weight*weightx),/scalable,colbar=[-2e7,6e7]
plt_image,transpose((radio3-radio7)*weight*weightx),/scalable,colbar=[-2e7,6e7]
plt_image,transpose((-(radio1-radio4-radio2+radio5)+(radio3-radio6))*weight*weightx),/scalable,colbar=[-2e7,6e7]
plt_image,transpose(((radio4-radio5)-(radio3-radio7))*weight*weightx),/scalable,colbar=[-2e7,6e7]
plt_image,transpose((-(radio1-radio4-radio2+radio5)+(radio7))*weight*weightx),/scalable,colbar=[-2e7,6e7]
;
;
openw,1,field+'.residual1mode.dat'
printf,1,radio3-radio7
close,1

;device,/close


device,filename='field3.maps.ps'
;
radio=dblarr(nx,ny,nz)
radio0=radio
openr,1,field+'.svdtest.diff.dat'
readf,1,radio
close,1

openr,1,field+'.svdtest.residual.dat'
readf,1,radio0
close,1

radio=reform(radio,nx*ny,nz)
radio0=reform(radio0,nx*ny,nz)

plt_image,transpose(radio*weight*weightx),/scalable,colbar=[-2e7,5e7]
plt_image,transpose(radio0*weight*weightx),/scalable,colbar=[-2e7,5e7]
plt_image,transpose((radio-radio0)*weight*weightx),/scalable,colbar=[-2e7,5e7]
;
device,/close

end

