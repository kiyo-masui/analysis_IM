pro cmp_plots

ssdir=['f4_09-14/','f4_15-20/','f4_21-26/','f4_27-32/','f4_33-36/']
sdir='/cita/scratch/cottontail/tchang/GBT/Processed_sim/sept06/'
fdir='/cita/scratch/cottontail/tchang/GBT/Processed/sept06/'
set_plot,'ps'
device,filename='comp.ps',/color, bits_per_pixel=8
for ifile=0,n_elements(ssdir)-1 do begin

npol=4
nx=12
ny=4
nz=560
radio=dblarr(npol,nx,ny,nz)
radiosim=radio

openr,1,fdir+ssdir(ifile)+'f4.radio2mode_weighted_aftsvd'
readf,1,radio
close,1

openr,1,sdir+ssdir(ifile)+'f4.radio2mode_weighted_aftsvd'
readf,1,radiosim
close,1

radio=reform(radio(0,*,*,*))
radiosim=reform(radiosim(0,*,*,*))

plt_image,transpose(reform(radio,nx*ny,nz)),/scalable,/colbar
plt_image,transpose(reform(radiosim,nx*ny,nz)),/scalable,/colbar
plt_image,transpose(reform((radiosim-radio)<0.00006>(-0.00001),nx*ny,nz)),/scalable,/colbar

endfor
device,/close

;
radio=dblarr(npol,nx,ny,nz)
radiosim=radio

openr,1,fdir+ssdir(ifile-1)+'f4.radiotot2mode_weighted_aftsvd'
readf,1,radio
close,1

openr,1,sdir+ssdir(ifile-1)+'f4.radiotot2mode_weighted_aftsvd'
readf,1,radiosim
close,1

radio=reform(radio(0,*,*,*))
radiosim=reform(radiosim(0,*,*,*))

;set_plot,'ps'
device,filename='comptot.ps',/color, bits_per_pixel=8
plt_image,transpose(reform(radio,nx*ny,nz)),/scalable,/colbar
plt_image,transpose(reform(radiosim,nx*ny,nz)),/scalable,/colbar
plt_image,transpose(reform((radiosim-radio)<0.00006>(-0.00001),nx*ny,nz)),/scalable,/colbar
device,/close
end
