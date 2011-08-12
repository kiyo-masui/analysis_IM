pro mk_hicube_half,Tsky,flagtot,noisecal,ra,dec,freq,pol,parallatic,radio=radio,weight=weight,fin=fin,field=field,kstart=kstart,kend=kend,flagnum=flagnum,Npol=Npol,noiseradio=noiseradio,f0=f0,fbin=fbin,dra=dra,ddec=ddec


; Tsky has dimension of (time, npol, redshift)
; binning into radio=(x,y,z)
;
;
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
;beam=15./60.d0 ; beam size
;dra=beam/2.d
;ddec=beam/2.d
;dra=3./60.d0
;ddec=3./60.d0
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz

; field center
;f0=[150.11667,2.2269444]
;fbin=[40,40]
;f0=[150.11635,2.2247776]
;fbin=[9,9]

print,'ra range',min(ra),max(ra)
print,'dec range',min(dec),max(dec)


; input parameters
; dimension of Tsky(nt,nrec,nfreq)
if not keyword_set(Npol) then Npol=4
nt=n_elements(ra) 


help,Tsky
ind=where(Tsky ne 0,cind)
print,'inside mk_hicube:  non-zero T elements:',cind
print,''
ind=where(radio ne 0,cind)
print,'inside mk_hicube:  non-zero radio elements:',cind
print,''
;
nx=fbin(0)
ny=fbin(1)

; for every timestamp, cal on-off runs first, then Npol, then IF

 
for i=0,nt-1 do begin
   x=(floor( (ra(i)-f0(0))*cos(dec(i)/!radeg)/dra)+nx/2)
   y=(floor( (dec(i)-f0(1))/ddec ) +ny/2)
   if (x ge 0 and x lt nx and y ge 0 and y lt ny) then begin
       for j=0,Npol-1 do begin
           flagslice=reform(flagtot(i,j,*))
           ind=where(flagslice gt 0,cind)
           if (cind gt 0) then begin
              radio(j,x,y,ind)=radio(j,x,y,ind)+Tsky(i,j,ind)*flagslice(ind)
              weight(j,x,y,ind)=weight(j,x,y,ind)+flagslice(ind)
             if keyword_set(noiseradio) then noiseradio(j,x,y,ind)=noiseradio(j,x,y,ind)+noisecal(i,j,ind)*flagslice(ind)
           endif
        endfor
    endif
endfor




end

