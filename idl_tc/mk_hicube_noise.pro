pro mk_hicube_noise,Tsky,flagtot,noisecal,ra,dec,freq,pol,parallatic,radio=radio,weight=weight,fin=fin,field=field,kstart=kstart,kend=kend,f0=f0,flagnum=flagnum,Npol=Npol,noiseradio=noiseradio

; bin the GBT data into cubes for cross-correlation
;
; must be run after calibration.pro
;
; Tsky:  output from calibration.pro, or read from file
; Radio:  input and output -- the binned radio(ra,dec,redshift) cube

; pixsel size
; should be half-beam spatially, and 2Mpc in redshift space

; for Deep2 
; f2 ra=[253.30618       251.55840]
;    dec=[35.188347  34.688137]
;
; f3  ra=[353.61133       351.37021]
;     dec=[0.43494931   -0.10343208]
;
; f4  ra=[38.051064       36.510399]
;     dec=[0.88520736    0.34644267]
;
;   
;parameters
;radio cube (ra,dec,redshift)

if not keyword_set(flagnum) then flagnum=0L

; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz

print,'ra range',min(ra),max(ra)
print,'dec range',min(dec),max(dec)

; field center
if (field eq 'f1') then begin
   ;f0=[214.250,52.50]
   f0=[96.47,59.75]
   fbin=[2,14]
endif 
if (field eq 'f2') then begin
   f0=[252.43229,34.938242]     ;field center coordinate
   fsize=[1.43867,0.500210]    ;apparent size on the sky
   fbin=[12,4]
;   radio=dblarr(12,4,nz)
endif
if (field eq 'f3') then begin
   f0=[352.49077,0.13736937]
   fsize=[2.24112, 0.538381]
   fbin=[18,4]
;   fbin=[18,5]
;   radio=dblarr(18,5,nz)
endif
if (field eq 'f4') then begin
   f0=[37.255199,0.58632926]
   fbin=[12,4]
endif
if (field eq '3c286') then begin
;   f0= [202.78453, 30.509155] 
   fbin=[4,4]
endif
;3c67       02:24:12.2973    +27:50:12.016
if (field eq '3c324') then begin
   f0 = [237.45427,21.427460]
endif

; input parameters
; dimension of Tsky(nt,nrec,nfreq)
Nif=1
if not keyword_set(Npol) then Npol=4
nwin=8
nfreq=2048*3/5*nwin
nrec=Nif*Npol


fdir='/home/scratch/tchang/'
if keyword_set(fin) then begin
   readcol,fdir+fin+'.params',nt,nrec,nfreq,format='i,i,i'
   Tsky=dblarr(nt,nrec,nfreq)
   ra=dblarr(nt,nrec)
   dec=ra
   pol=intarr(nt,nrec)
   freq=dblarr(nt,nrec,3)
   openr,1,fdir+fin+'.t'
   readf,1,Tsky
   close,1
   help,Tsky
   openr,1,fdir+fin+'.ra'
   readf,1,ra
   close,1
   openr,1,fdir+fin+'.dec'
   readf,1,dec
   close,1
   openr,1,fdir+fin+'.pol'
   readf,1,pol
   close,1
   openr,1,fdir+fin+'.freq'
   readf,1,freq
   close,1
endif else begin
   Tsky=Tsky
   ra=ra
   dec=dec
   freq=freq
   nt=n_elements(Tsky)/nfreq/nrec 
endelse

help,Tsky
ind=where(Tsky ne 0,cind)
print,'inside mk_hicube:  non-zero T elements:',cind
print,''
ind=where(radio ne 0,cind)
print,'inside mk_hicube:  non-zero radio elements:',cind
print,''
;
nt=nt(0)
nfreq=nfreq(0)
Nif=Nif(0)
help,Nif
help,nfreq
help,nt
nx=fbin(0)
ny=fbin(1)
dnu=50e6/double(nfreq)
dt=0.5
sigma=5.

; for every timestamp, cal on-off runs first, then Npol, then IF
;
; first decode the frequency and compute the distance
Da_z,(freq0/freq-1.),distslice
distbin=(floor((distslice-d0)/dd)) ;d0 is at bin0
;
print,'distbin',max(distbin),min(distbin)
;
; figure out which freq channels are overlapping
; each IF has 2048 channels, 50 MHz, only the center 10 MHz does not
; overlap
;overlap=make_array(nfreq,/double,value=0.25d)  ;consider LL+RR
;overlap(nfreq/2-nfreq/10:nfreq/2+nfreq/10-1)=0.5d0 
 
for i=0,nt-1 do begin
   x=(floor( (ra(i)-f0(0))*cos(dec(i)/!radeg)/dra)+nx/2)
   y=(floor( (dec(i)-f0(1))/ddec ) +ny/2)
   if (x ge 0 and x lt nx and y ge 0 and y lt ny) then begin
       for j=0,Npol-1 do begin
           zslice=distbin
           flagslice=reform(flagtot(i,j,*))
                                ;
                                ; compute the corresponding frequency and distance
                                ;
           
           for k=0,nz-1 do begin
         ;
         ; find out the range that goes into the same redshift bin
         ; and do a median on that
               tcind=0
               tvalue=0.
;               ind=where(zslice eq k and flagslice eq 1 and finite(Tsky(i,j,*)),cind)
               ind=where(zslice eq k and flagslice eq 1,cind)
               ;if j eq 0 then print,'k,cind:',k,cind
               if (cind gt 0) then begin
                  tslice=Tsky(i,j,ind)
                  tvalue=total(tslice)
                  nslice=noisecal(i,j,ind)
                  nvalue=total(nslice)
                  ;
                  if (abs(tvalue) lt sigma*sqrt(double(cind)/dnu*dt)) then begin
                  ; do a sigmacut
                     radio(j,x,y,k)=radio(j,x,y,k)+tvalue
                     weight(j,x,y,k)=weight(j,x,y,k)+double(cind)
                     noiseradio(j,x,y,k)=noiseradio(j,x,y,k)+nvalue
                  endif
               endif
            endfor
        endfor
   endif
endfor

; subtract off the foreground
; calculate the median in some freq range and subtract it off
; not counting the flagged points
;dnfreq=40
;dfreq=nz/dnfreq
;lrfreqmedian=dblarr(nx,ny,dnfreq)
;rlfreqmedian=dblarr(nx,ny,dnfreq)
;for x=0,nx-1 do begin
;   for y=0,ny-1 do begin
;      for k=0,dnfreq-1 do begin
;         medd=medianradio(1,x,y,k*dfreq:(k+1)*dfreq-1)

;         if (cind gt 1) then begin
;            med=median(medd(ind),/double)
;;            for ll=0,dfreq-1 do begin
;;               Tsky(i,j,k*dfreq+ll)=Tsky(i,j,k*dfreq+ll)-med
;;            endfor
;;         endif
;;      endfor
;;   endfor
;;endfor

; do flagging
;nfreq=2048
;dnu=50e6/double(nfreq)
;dt=0.5d
;for x=0,nx-1 do begin
;   for y=0,ny-1 do begin
;      ind=where(reform(weight(0,x,y,*)) gt 0,cind)
;      if cind gt 0 then begin
;         r=(radio(1,x,y,ind)^2+radio(2,x,y,ind)^2.)/(radio(0,x,y,ind)*radio(3,x,y,ind))
;         print,'** r :',r
;         print,'** sigma :',25.*2./(reform(weight(0,x,y,ind))*dnu*dt)
;         ind2=where(r gt 25.*2./(reform(weight(0,x,y,ind))*dnu*dt),cind2)
;         if cind2 gt 0 then begin
;            weight(*,x,y,ind(ind2))=0.
;            radio(*,x,y,ind(ind2))=0.
;            flagnum=flagnum+cind2
;         endif
;      endif
;   endfor
;endfor
;
;print,'** flagged points',flagnum
;print,'** flagged percentage ',double(flagnum)/double(n_elements(radio(0,*,*,*)))



;openw,1,fdir+fin+'.radio'
;printf,1,radio
;close,1


end

