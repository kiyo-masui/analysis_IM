pro mk_hicube_sim_fg,Tsky,flagtot,ra,dec,freq,pol,parallatic,radio=radio,weight=weight,fin=fin,field=field,kstart=kstart,kend=kend,f0=f0,flagnum=flagnum,Npol=Npol,noiseradio=noiseradio

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
;dra=3./60.d0
;ddec=3./60.d0
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz

print,'ra range',min(ra),max(ra)
print,'dec range',min(dec),max(dec)

; field center
f0=[150.11635,2.2247776]
fbin=[9,9]



; input parameters
; dimension of Tsky(nt,nrec,nfreq)
Nif=1
if not keyword_set(Npol) then Npol=4
nwin=8
nfreq=2048*3/5*nwin
nrec=Nif*Npol


help,Tsky
ind=where(Tsky ne 0,cind)
print,'inside mk_hicube:  non-zero T elements:',cind
print,''
ind=where(radio ne 0,cind)
print,'inside mk_hicube:  non-zero radio elements:',cind
print,''
;
nt=n_elements(ra)
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
                  ;nslice=noisecal(i,j,ind)
                  ;nvalue=total(nslice)
                  ;
                  ;if (abs(tvalue) lt sigma*sqrt(double(cind)/dnu*dt)) then begin
                  ; do a sigmacut
                     radio(j,x,y,k)=radio(j,x,y,k)+tvalue
                     weight(j,x,y,k)=weight(j,x,y,k)+double(cind)
                  ;   noiseradio(j,x,y,k)=noiseradio(j,x,y,k)+nvalue
                  ;endif
               endif
            endfor
        endfor
   endif
endfor



end

