pro combine_scans,radio,weight,t_off,radio_on,weight_on,field=field,fn=fn,sigmacut=sigmacut,dname=dname,svdfg=svdfg,fitsfile=fitsfile,sim=sim,signal=signal,daisy=daisy,drift=drift

; field:  data field
; fn:  fits file name

if not keyword_set(field) then field='1hr'  ; for the radio/optical cube setup
if not keyword_set(signal) then signal='100muK'
if not keyword_set(dname) then dname='june14'
;if keyword_set(daisy) then obsmode='daisy' else obsmode='drift' ;file IO only
obsmode = 'azel'
if not keyword_set(daisy) then drift=drift

if dname eq 'june14' and field eq '1hr' then begin ; changed all this -km
   ffdir1=['49-56','57-64','65-72','73-80','81-88','89-96','97-104']
   ffdir2=['105-112','113-120','121-128','129-136','137-144','145-152','153-160']
endif
if dname eq 'june14' and field eq '22hr' then begin ; changed all this -km
   ffdir1=['9-16','17-24','25-32']
   ffdir2=['33-40','41-48']
endif

;raidgmrt = '/mnt/raid-project/gmrt' ; added this -km
; raidgmrt = /cita/d/raid-project/gmrt
raidgmrt = '/cita/d/raid-cita/tchang/'

set_plot,'ps'
device,/color, bits_per_pixel=8

; GBT parameters
nx=40
ny=20
nz=560
Npol=4
ffdir=ffdir1
nfile=n_elements(ffdir)
;
radio=dblarr(Npol,nx,ny,nz)
radiotot=radio
weight=radio
weighttot=weight
radiotot20mode=radio
weighttot20mode=weight
radiotot25mode=radio
weighttot25mode=weight
radiotot10mode=radio
weighttot10mode=weight
radiotot5mode=radio
weighttot5mode=weight
radiotot0mode=radio
weighttot0mode=weight
;
for ifile=0,nfile-1 do begin
   ; only changed the root of this line -km
   fdir=raidgmrt+'/wiggleZ/processed/'+dname+'/'+ffdir(ifile)+'/'
   ;
   radio=dblarr(Npol,nx,ny,nz)
   weight=radio
   radio20mode=dblarr(Npol,nx,ny,nz)
   weight20mode=radio
   radio25mode=dblarr(Npol,nx,ny,nz)
   weight25mode=radio
   radio10mode=dblarr(Npol,nx,ny,nz)
   weight10mode=radio
   radio5mode=dblarr(Npol,nx,ny,nz)
   weight5mode=radio
   radio0mode=dblarr(Npol,nx,ny,nz)
   weight0mode=radio
   filenum=strtrim(string(ifile),2)
   ;
   ;
   ;
   openr,1,fdir+field+'.radio_aftsvd'
   readf,1,radio
   close,1
   ;
   openr,1,fdir+field+'.weight_aftsvd'
   readf,1,weight
   close,1
   ;
   openr,1,fdir+field+'.radio20mode_aftsvd'
   readf,1,radio20mode
   close,1
   ;
   openr,1,fdir+field+'.weight20mode_aftsvd'
   readf,1,weight20mode
   close,1
   ;
   openr,1,fdir+field+'.radio25mode_aftsvd'
   readf,1,radio25mode
   close,1
   ;
   openr,1,fdir+field+'.weight25mode_aftsvd'
   readf,1,weight25mode
   close,1
   ;
   openr,1,fdir+field+'.radio10mode_aftsvd'
   readf,1,radio10mode
   close,1
   ;
   openr,1,fdir+field+'.weight10mode_aftsvd'
   readf,1,weight10mode
   close,1
   ;
   openr,1,fdir+field+'.radio5mode_aftsvd'
   readf,1,radio5mode
   close,1
   ;
   openr,1,fdir+field+'.weight5mode_aftsvd'
   readf,1,weight5mode
   close,1
   ;
   openr,1,fdir+field+'.radio0mode_aftsvd'
   readf,1,radio0mode
   close,1
   ;
   openr,1,fdir+field+'.weight0mode_aftsvd'
   readf,1,weight0mode
   close,1
   ;
   radiotot=radiotot+radio
   weighttot=weighttot+weight
   radiotot20mode=radiotot20mode+radio20mode
   weighttot20mode=weighttot20mode+weight20mode
   radiotot25mode=radiotot25mode+radio25mode
   weighttot25mode=weighttot25mode+weight25mode
   radiotot10mode=radiotot10mode+radio10mode
   weighttot10mode=weighttot10mode+weight10mode
   radiotot5mode=radiotot5mode+radio5mode
   weighttot5mode=weighttot5mode+weight5mode
   radiotot0mode=radiotot0mode+radio0mode
   weighttot0mode=weighttot0mode+weight0mode
   ;
endfor
;
;
; total file outputs per day
openw,1,fdir+field+'.radiotot_aftsvd_half'
printf,1,radiotot
close,1
;
openw,1,fdir+field+'.weighttot_aftsvd_half'
printf,1,weighttot
close,1
;
ind=where(weighttot gt 0,cind)
radiotot(ind)=radiotot(ind)/weighttot(ind)
;
openw,1,fdir+field+'.radiotot_weighted_aftsvd_half'
printf,1,radiotot
close,1
;
;

openw,1,fdir+field+'.radiotot20mode_aftsvd_half'
printf,1,radiotot20mode
close,1
;
openw,1,fdir+field+'.weighttot20mode_aftsvd_half'
printf,1,weighttot20mode
close,1
;
ind=where(weighttot20mode gt 0,cind)
radiotot20mode(ind)=radiotot20mode(ind)/weighttot20mode(ind)
;
openw,1,fdir+field+'.radiotot20mode_weighted_aftsvd_half'
printf,1,radiotot20mode
close,1
;
openw,1,fdir+field+'.radiotot25mode_aftsvd_half'
printf,1,radiotot25mode
close,1
;
openw,1,fdir+field+'.weighttot25mode_aftsvd_half'
printf,1,weighttot25mode
close,1
;
ind=where(weighttot25mode gt 0,cind)
radiotot25mode(ind)=radiotot25mode(ind)/weighttot25mode(ind)
;
openw,1,fdir+field+'.radiotot25mode_weighted_aftsvd_half'
printf,1,radiotot25mode
close,1
;
;
openw,1,fdir+field+'.radiotot10mode_aftsvd_half'
printf,1,radiotot10mode
close,1
;
openw,1,fdir+field+'.weighttot10mode_aftsvd_half'
printf,1,weighttot10mode
close,1
;
ind=where(weighttot10mode gt 0,cind)
radiotot10mode(ind)=radiotot10mode(ind)/weighttot10mode(ind)
;
openw,1,fdir+field+'.radiotot10mode_weighted_aftsvd_half'
printf,1,radiotot10mode
close,1
;
;
openw,1,fdir+field+'.radiotot5mode_aftsvd_half'
printf,1,radiotot5mode
close,1
;
openw,1,fdir+field+'.weighttot5mode_aftsvd_half'
printf,1,weighttot5mode
close,1
;
ind=where(weighttot5mode gt 0,cind)
radiotot5mode(ind)=radiotot5mode(ind)/weighttot5mode(ind)
;
openw,1,fdir+field+'.radiotot5mode_weighted_aftsvd_half'
printf,1,radiotot5mode
close,1
;
;
openw,1,fdir+field+'.radiotot0mode_aftsvd_half'
printf,1,radiotot0mode
close,1
;
openw,1,fdir+field+'.weighttot0mode_aftsvd_half'
printf,1,weighttot0mode
close,1
;
ind=where(weighttot0mode gt 0,cind)
radiotot0mode(ind)=radiotot0mode(ind)/weighttot0mode(ind)
;
openw,1,fdir+field+'.radiotot0mode_weighted_aftsvd_half'
printf,1,radiotot0mode
close,1
;

;
;
end

