pro docalib_gensim,weight,t_off,radio_on,weight_on,field=field,fn=fn,sigmacut=sigmacut,dname=dname,svdfg=svdfg,fitsfile=fitsfile,signal=signal,daisy=daisy,drift=drift

; field:  data field
; fn:  fits file name

if not keyword_set(field) then field='zcosmos'  ; for the radio/optical cube setup
if not keyword_set(dname) then dname='02'
if not keyword_set(signal) then signal='100muK'
if not keyword_set(daisy) then drift=drift

if dname eq '02' then begin
;   ffdir=['30','scan9-17', 'scan18-26','scan27-28'] 
   if keyword_set(daisy) then ffdir=['30']
   if keyword_set(drift) then ffdir=['scan9-17', 'scan18-26','scan27-28'] ;
endif

set_plot,'ps'
!p.multi=[0,8,1]

; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz
nfreq=2048
Npol=4

; field center
f0=[150.11635,2.2247776]
fbin=[9,9]
;
;
; parameters and matrices
nx=fbin(0)
ny=fbin(1)
npol=4
simradio=dblarr(npol,nx,ny,nz)
radio=dblarr(npol,nx,ny,nz)
weight=radio

; parameters
nfile=n_elements(ffdir)
nn=16320L  ;255.*8*8  8 freq windows and 4 pol x 2 on-off
nscan=9
;nscan=1
ntime=255
npol=4
nwin=8
nfreq=2048
;
; read in simulated HI file (in z=560 format and without unit P=T /gain)
gdir='~/projects/GBT/pros/cosmos/Calibration/'+dname+'/'
openr,1,gdir+'radiosim_'+signal
readf,1,simradio
close,1
;
for ifile=0,nfile-1 do begin
   ;
   nscan=9
   if ffdir(ifile) eq '30' then nscan=1
   if ffdir(ifile) eq 'scan27-28' then nscan=2
   fdir='/cita/d/raid-project/gmrt/tchang/GBT/zCOSMOS/Processed/'+dname+'/'+ffdir(ifile)+'/'
   ;
   ;
   ;
   print,'fdir: ',fdir
      ;
      ; collect the frequencies,ra,dec as well
   radiotf=dblarr(ntime*nscan,npol,nfreq*3/5*nwin)
   freqtot=dblarr(nfreq*3/5*nwin)
   ra=dblarr(ntime*nscan)
   dec=dblarr(ntime*nscan)
   sim_Psky=dblarr(ntime*nscan,npol,nfreq*3/5*nwin)
   sim_count=sim_Psky
      ;
      ; read in the data files
   openr,1,fdir+'freqtot'
   readf,1,freqtot
   close,1
      ;
   openr,1,fdir+'ra'
   readf,1,ra
   close,1
   ;
   openr,1,fdir+'dec'
   readf,1,dec
   close,1
      ;   
      ; convert freq to redshift
   Da_z,(freq0/freqtot-1.),dist
   distbin=(floor((dist-d0)/dd))                       ;d0 is at bin0
      ;
      ; ra and dec to x y
   x=(floor( (ra-f0(0))*cos(dec/!radeg) /dra)+nx/2) ;>0<(fbin(0)-1)
   y=(floor( (dec-f0(1))/ddec) +ny/2)               ;>0<(fbin(1)-1)
      ;
;   print,'x',x
;   print,'y',y
;   return

                                ; from ra and dec and freq work out
                                ; the corresponding density field and
                                ; thus HI temperature field it
                                ; corresponds to
   for ipol=0,npol-1 do begin
      for it=0,ntime*nscan-1 do begin
         for ifreq=0,nfreq*3/5*nwin-1 do begin
;            print,ipol,x(it),y(it),distbin(ifreq)
            if x(it) ge 0 and x(it) lt fbin(0) and y(it) ge 0 and y(it) lt fbin(1) then $
               radiotf(it,ipol,ifreq)=simradio(ipol,x(it),y(it),distbin(ifreq))
         endfor
      endfor
   endfor
      ;
      ;
      ; output file
   openw,1,fdir+'simradio_tf_b4svd_'+signal
   printf,1,radiotf
   close,1
                                ;
   device,filename=fdir+'simradio_tf_'+signal+'.ps',/color, bits_per_pixel=8
   plt_image,transpose(reform(radiotf(*,0,*))),/scalable,/colbar
   device,/close
   ;
   flagtot=make_array(ntime*nscan,npol,nfreq*3/5*nwin,value=1.)
   mk_hicube_sim_fg,radiotf,flagtot,ra,dec,freqtot,pol,parallatic,$
                    radio=radio,weight=weight,field=field
      ;
   device,filename=fdir+'radiosim_'+signal+'.ps',/color, bits_per_pixel=8
   radiotmp=reform(radio,npol,nx*ny,nz)
   plt_image,transpose(reform(radiotmp(0,*,*))),/scalable,/colbar
   device,/close
endfor
;
openw,1,fdir+field+'.simradiotot_'+signal
printf,1,radio
close,1
;
openw,1,fdir+field+'.simweighttot_'+signal
printf,1,weight
close,1
;
radio2=radio
ind=where(weight gt 0, cind)
radio(ind)=radio(ind)/weight(ind)
;
openw,1,fdir+field+'.simradiotot_weighted_'+signal
printf,1,radio
close,1
;
set_plot,'ps'
!p.multi=[0,2,1]
device,filename='radiosim_'+signal+'.ps',/color, bits_per_pixel=8
radio=reform(radio,npol,nx*ny,nz)
plt_image,transpose(reform(radio(0,*,*))),/scalable,/colbar
radio2=reform(radio2,npol,nx*ny,nz)
plt_image,transpose(reform(radio2(0,*,*))),/scalable,/colbar
device,/close
!p.multi=0
;
;
;
end
