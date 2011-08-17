pro dosvd_gain,radio,weight,t_off,radio_on,weight_on,field=field,dname=dname,svdfg=svdfg

; field:  data field
; fn:  fits file name

if not keyword_set(field) then field='f3'  ; for the radio/optical cube setup
if not keyword_set(dname) then dname='aug29'

if dname eq 'sept07' then begin
;   ffdir=['05-08']
   ffdir=['09-14', '15-20', '21-26','27-29'] ;sept07
   dnum='d4'
endif
if dname eq 'aug29' then begin
   ffdir=['05-10', '15-20', '21-26', '27-32', '37-42', '47-52', '53-58'] ;aug29
   dnum='d1'
endif
if dname eq 'aug30' then begin
   ffdir=['9-14', '15-20', '21-26', '27-32','33-38', '39-44', '45-50'] ;aug30
   dnum='d2'
endif
if dname eq 'sept17' then begin
   ffdir=['12-17','18-23','24-29'];,'30-35','36-41','42-47','48-53'] ;sept 17
   dnum='d5'
endif

set_plot,'ps'

; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq_rest=1420.405751 ;MHz
nfreq=2048
Npol=4

; field center
if (field eq 'nissim') then begin
   f0=[40.801455,-5.8477362]
   fbin=[12,4]
;   radio=dblarr(12,4,nz)
endif
if (field eq 'f1') then begin
;   f0=[214.250,52.50]
   ; use galactic coordiate
   f0=[96.47,59.75]
   fbin=[2,14]
endif 
if (field eq 'f2') then begin
   f0=[252.43229,34.938242]     ;field center coordinate
   fsize=[1.43867,0.500210]    ;apparent size on the sky
   fbin=[12,4]
endif
if (field eq 'f3') then begin
   f0=[352.49077,0.13736937]
   fsize=[2.24112, 0.538381]
   fbin=[18,5]
endif
if (field eq 'f4') then begin
   f0=[37.255199,0.58632926]
   fbin=[12,4]
endif
if (field eq '3c286') then begin
   f0= [202.78453, 30.509155] 
endif
;3c67       02:24:12.2973    +27:50:12.016
if (field eq '3c324') then begin
   f0 = [237.45427,21.427460]
endif
;
;
; parameters and matrices
Npol=4
nx=fbin(0)
ny=fbin(1)
radio=dblarr(Npol,nx,ny,nz)
weight=radio
noiseradio=dblarr(Npol,nx,ny,nz)

; parameters
nfile=n_elements(ffdir)
npol=4
;
; read in data file
;fdir='/cita/scratch/cottontail/tchang/GBT/Processed/'+dname+'/f3_'+ffdir(nfile-1)+'/'
fdir='/Users/tzu/Projects/GBT/Processed/'+dname+'/f3_'+ffdir(nfile-1)+'/'
;
openr,1,fdir+field+'.radiotot_weighted_b4svd'
readf,1,radio
close,1
;
openr,1,fdir+field+'.weighttot_b4svd'
readf,1,weight
close,1
;
openr,1,fdir+field+'.noisetot_weighted_b4svd'
readf,1,noiseradio
close,1
;
; reform the arrays
radio=reform(radio,Npol,nx*ny,nz)
weight=reform(weight,Npol,nx*ny,nz)
noiseradio=reform(noiseradio,Npol,nx*ny,nz)
help,weight
help,radio
;
; read in gain_x and gain_z
gainx=dblarr(Npol,nx*ny)
gainz=dblarr(Npol,nz)
;
;gdir='~/projects/GBT/pros/'+dname+'/'
gdir='/Users/tzu/Projects/GBT/pros/'+dname+'/'
openr,1,gdir+'deep2_gain_x_b4svd.txt'
readf,1,gainx
close,1
;
openr,1,gdir+'deep2_gain_z_b4svd.txt'
readf,1,gainz
close,1
;
Tsky=dblarr(Npol,nx*ny,nz)
; divide out g(z)g(x) from Psky(x,z)=g(z)g(x)Tsky(x,z) to get Tsky(x,z) 
for ix=0,nx*ny-1 do begin
   for iz=0,nz-1 do begin
      for ipol=0,npol-1 do begin
         if (gainx(ipol,ix)*gainz(ipol,iz) ne 0.) then begin
            Tsky(ipol,ix,iz)=radio(ipol,ix,iz)/gainx(ipol,ix)/gainz(ipol,iz) 
         endif else begin
            weight(ipol,ix,iz)=0.
         endelse
      endfor
   endfor
endfor
;
; get rid of the empty dec strips
ny=3
Tsky=Tsky(*,0:nx*ny-1,*)
weight=weight(*,0:nx*ny-1,*)
help,Tsky
print,max(Tsky),min(Tsky)
;read,test
help,weight
;
; write outputs
openw,1,fdir+'Tskytot_weighted_b4svd'
printf,1,Tsky
close,1
;
device,filename=fdir+'Tskytot_weighted_b4svd.ps'
plt_image,transpose(reform(Tsky(0,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(1,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(2,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(3,*,*))),/scalable,/colbar
device,/close
; 

; for svd, do -1 (here, is it equivalent to -1/(g(z)*g(x))...)
; first calculate mean Tsky(z)
ind=where(weight gt 0,cind)
print,'non-zero weight',cind
;read,test
Tmean=dblarr(npol,nz)
for ipol=0,npol-1 do begin
   for iz=0,nz-1 do begin
      ind=where(weight(ipol,*,iz) gt 0,cind)
      if (cind gt 0) then begin
         Tmean(ipol,iz)=median(Tsky(ipol,ind,iz),/double)
         ;print,tmean
         if Tmean(ipol,iz) ne 0. then Tsky(ipol,ind,iz)=Tsky(ipol,ind,iz)/Tmean(ipol,iz)-1.
      endif
   endfor
endfor

; write outputs
openw,1,fdir+'Tskytot_weighted_b4svd_sub1'
printf,1,Tsky
close,1
;
openw,1,fdir+'Tskytot_mean'
printf,1,Tmean
close,1
;
device,filename=fdir+'Tskytot_weighted_b4svd_sub1.ps'
plt_image,transpose(reform(Tsky(0,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(1,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(2,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(3,*,*))),/scalable,/colbar
device,/close
; 
;
; do sigma flagging
if keyword_set(doscut) then begin
dnu=50e6/double(nfreq)
dt=0.5d
sigma=1./sqrt(dnu*dt);*30  ;Tsys=30K
sigma_th=12.
ind=where(weight gt 0,cind)
sigma2d=dblarr(npol,nx*ny,nz)
sigma2d(ind)=sigma_th*sigma/sqrt(weight(ind))
;
for ipol=0,npol-1 do begin
   for ix=0,nx*ny-1 do begin
      for iz=0,nz-1 do begin
         if (abs(Tsky(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then begin
            Tsky(ipol,ix,iz)=0
            weight(ipol,ix,iz)=0
         endif
      endfor
   endfor
endfor
;
device,filename=fdir+'Tskytot_weighted_b4svd_sub1_scut.ps'
plt_image,transpose(reform(Tsky(0,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(1,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(2,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(3,*,*))),/scalable,/colbar
device,/close
;
endif
; 
;
; fill in the holes in frequency (two freq columns)
; first calculate P(x)
;psub=dblarr(npol,nx*ny)
;for ipol=0,npol-1 do begin
;   for ix=0,nx*ny-1 do begin
;      psub(ipol,ix)=median(Tsky(ipol,ix,*),/double)
;   endfor
;endfor
;
zlist=0
for iz=10,nz-10 do begin
   ind=where(weight(*,*,iz) gt 0,cind)
   if (cind lt 100) then begin
;      Tsky(*,*,iz)=psub
      zlist=[zlist,iz]
   endif
endfor
zlist=zlist(1:n_elements(zlist)-1)
zlist=[zlist,263,262,261,260,259]
if dname eq 'aug30' then zlist=[zlist,findgen(20)+nz-20]
;
help,zlist
print,zlist
for iz=0,n_elements(zlist)-1 do begin
   list=zlist(iz)+findgen(40)-20.<(nz-1)
   for ipol=0,npol-1 do begin
      for ix=0,nx*ny-1 do begin
         Tsky(ipol,ix,zlist(iz))=median(Tsky(ipol,ix,list),/double)
      endfor
   endfor
endfor
;
device,filename=fdir+'Tskytot_weighted_b4svd_sub1_fill.ps'
plt_image,transpose(reform(Tsky(0,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(1,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(2,*,*))),/scalable,/colbar
plt_image,transpose(reform(Tsky(3,*,*))),/scalable,/colbar
device,/close
;
;
; svd
radiosvd=dblarr(Npol,nx*ny,nz)
radiores=radiosvd
radiosvd1mode=radiosvd
radiores1mode=radiosvd
radiosvd2mode=radiosvd
radiores2mode=radiosvd
radio2dtot=radiosvd
radiosvd5mode=radiosvd
radiores5mode=radiosvd
radiosvd10mode=radiosvd
radiores10mode=radiosvd
;
for ii=0,npol-1 do begin
   ;
   radio2d=reform(Tsky(ii,*,*))
   ;
   la_svd,radio2d, s,u,v,/double,status=status
   print,'SVD status:',status
   ;
;   help,s
;   help,u
;   help,v
   ;
   w=dblarr(n_elements(s))
   w2=w
   w5=w
   w10=w
   w(0)=s(0)
   w2(0:1)=s(0:1)
   w5(0:4)=s(0:4)
   w10(0:9)=s(0:9)
   ;
   radiosvd(ii,*,*) = u ## diag_matrix(s) ## transpose(v)
   radiosvd1mode(ii,*,*) =u ## diag_matrix(w) ## transpose(v)
   radiosvd2mode(ii,*,*) = u ## diag_matrix(w2) ## transpose(v)
   radiosvd5mode(ii,*,*) =u ## diag_matrix(w5) ## transpose(v)
   radiosvd10mode(ii,*,*) = u ## diag_matrix(w10) ## transpose(v)
   ;
endfor
;
radiores = Tsky-radiosvd
radiores1mode = Tsky-radiosvd1mode
radiores2mode = Tsky-radiosvd2mode
radiores5mode = Tsky-radiosvd5mode
radiores10mode = Tsky-radiosvd10mode
;
device,filename=fdir+'Tskytot_weighted_svd1.ps'
;plt_image, transpose(reform(Tsky(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiosvd(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiores(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiosvd1mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores1mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiosvd2mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores2mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiosvd5mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores5mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiosvd10mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores10mode(2,*,*))),/scalable,/colbar
device,/close
;
; do sigma flagging
dnu=50e6/double(nfreq)
dt=0.5d
sigma=1./sqrt(dnu*dt);*30  ;Tsys=30K
sigma_th=8.
ind=where(weight gt 0,cind)
sigma2d=dblarr(npol,nx*ny,nz)
sigma2d(ind)=sigma_th*sigma/sqrt(weight(ind))
;
fcount=0L
for ipol=0,npol-1 do begin
   for ix=0,nx*ny-1 do begin
      for iz=0,nz-1 do begin
         if (abs(radiores5mode(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then begin
            radiores5mode(ipol,ix,iz)=0
            weight(ipol,ix,iz)=0
            Tsky(ipol,ix,iz)=0.
            fcount=fcount+1
         endif
         if (abs(radiores1mode(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then $
            radiores1mode(ipol,ix,iz)=0
         if (abs(radiores2mode(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then $
            radiores2mode(ipol,ix,iz)=0
         if (abs(radiores10mode(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then $
            radiores10mode(ipol,ix,iz)=0
      endfor
   endfor
endfor
;
print,'flagged ponits after SVD1:',fcount,double(fcount)/n_elements(sigma2d)
ind=where(weight eq 0,cind)
print,'total masked:',cind,double(cind)/n_elements(sigma2d)
;
device,filename=fdir+'Tskytot_weighted_svd1_scut.ps'
plt_image, transpose(reform(radiores1mode(2,*,*)*weight(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores2mode(2,*,*)*weight(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores5mode(2,*,*)*weight(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores10mode(2,*,*)*weight(2,*,*))),/scalable,/colbar
device,/close

; fill in the svd 5-mode values to the holes
ind=where(weight eq 0, cind)
Tsky(ind)=radiosvd5mode(ind)

device,filename=fdir+'Tskytot_weighted_svd1_scut_fill.ps'
plt_image, transpose(reform(Tsky(0,*,*)*weight(0,*,*))),/scalable,/colbar
plt_image, transpose(reform(Tsky(2,*,*)*weight(1,*,*))),/scalable,/colbar
plt_image, transpose(reform(Tsky(2,*,*)*weight(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(Tsky(3,*,*)*weight(3,*,*))),/scalable,/colbar
device,/close

; do frequency weighting
;if keyword_set(freqw) then begin
freqweight=dblarr(npol,nz)
for ii=0,npol-1 do begin
   for kk=0,nz-1 do begin
      ind=where(weight(ii,*,kk) gt 0,cind)
      if cind gt 1 then begin
         freqweight(ii,kk)=variance(radiores5mode(ii,ind,kk))
         Tsky(ii,ind,kk)=Tsky(ii,ind,kk)/freqweight(ii,kk)
      endif
   endfor
endfor
;endif
;
device,filename=fdir+'Tskytot_weighted_svd1_scut_fill_fweight.ps'
plt_image, transpose(reform(Tsky(0,*,*))),/scalable,/colbar
plt_image, transpose(reform(Tsky(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(Tsky(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(Tsky(3,*,*))),/scalable,/colbar
device,/close
;
;
;for ipol=0,npol-1 do begin
;   for ix=0,nx*ny-1 do begin
;      psub(ipol,ix)=median(Tsky(ipol,ix,*),/double)
;   endfor
;endfor
;for iz=0,n_elements(zlist)-1 do begin
;   Tsky(*,*,zlist(iz))=psub
;endfor
;zlist=zlist(1:n_elements(zlist)-1)
;for iz=0,n_elements(zlist)-1 do begin
;   list=zlist(iz)+findgen(20)-10.
;   for ipol=0,npol-1 do begin
;      for ix=0,nx*ny-1 do begin
;         Tsky(ipol,ix,zlist(iz))=median(Tsky(ipol,ix,list),/double)
;      endfor
;   endfor
;endfor
;
;plt_image, transpose(reform(Tsky(2,*,*))),/scalable,/colbar
;
for ii=0,npol-1 do begin
            ;
   radio2d=reform(Tsky(ii,*,*))
   la_svd,radio2d, s,u,v,/double,status=status
   print,'SVD status:',status
;   help,s
;   help,u
;   help,v
   w=dblarr(n_elements(s))
   w2=w
   w5=w
   w10=w
   w(0)=s(0)
   w2(0:1)=s(0:1)
   w5(0:4)=s(0:4)
   w10(0:9)=s(0:9)
                                ;
   radiosvd(ii,*,*) = u ## diag_matrix(s) ## transpose(v)
   radiosvd1mode(ii,*,*) =u ## diag_matrix(w) ## transpose(v)
   radiosvd2mode(ii,*,*) = u ## diag_matrix(w2) ## transpose(v)
   radiosvd5mode(ii,*,*) =u ## diag_matrix(w5) ## transpose(v)
   radiosvd10mode(ii,*,*) = u ## diag_matrix(w10) ## transpose(v)
   ;
endfor
;
radiores = Tsky-radiosvd
radiores1mode = Tsky-radiosvd1mode
radiores2mode = Tsky-radiosvd2mode
radiores5mode = Tsky-radiosvd5mode
radiores10mode = Tsky-radiosvd10mode
;
; de-freq weight
;if keyword_set(freqw) then begin
for ii=0,npol-1 do begin
   for kk=0,nz-1 do begin
      ind=where(weight(ii,*,kk) gt 0,cind)
      if cind gt 1 then begin
         radiores5mode(ii,ind,kk)=radiores5mode(ii,ind,kk)*freqweight(ii,kk)
         radiores10mode(ii,ind,kk)=radiores10mode(ii,ind,kk)*freqweight(ii,kk)
         radiores1mode(ii,ind,kk)=radiores1mode(ii,ind,kk)*freqweight(ii,kk)
         radiores2mode(ii,ind,kk)=radiores2mode(ii,ind,kk)*freqweight(ii,kk)
      endif
   endfor
endfor
;endif
;
device,filename=fdir+'Tskytot_weighted_svd2.ps'
plt_image, transpose(reform(radiores1mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores2mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores5mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores10mode(2,*,*))),/scalable,/colbar
device,/close
;
;
; do sigma flagging
dnu=50e6/double(nfreq)
dt=0.5d
sigma=1./sqrt(dnu*dt);*30  ;Tsys=30K
sigma_th=5.
ind=where(weight gt 0,cind)
sigma2d=dblarr(npol,nx*ny,nz)
sigma2d(ind)=sigma_th*sigma/sqrt(weight(ind))
;
fcount=0L
for ipol=0,npol-1 do begin
   for ix=0,nx*ny-1 do begin
      for iz=0,nz-1 do begin
         if (abs(radiores5mode(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then begin
            radiores5mode(ipol,ix,iz)=0
            weight(ipol,ix,iz)=0
            fcount=fcount+1
         endif
         if (abs(radiores1mode(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then $
            radiores1mode(ipol,ix,iz)=0
         if (abs(radiores2mode(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then $
            radiores2mode(ipol,ix,iz)=0
         if (abs(radiores10mode(ipol,ix,iz)) gt sigma2d(ipol,ix,iz)) then $
            radiores10mode(ipol,ix,iz)=0
      endfor
   endfor
endfor
;
print,'flagged points after SVD2:',fcount,double(fcount)/n_elements(weight)
ind=where(weight eq 0,cind)
print,'total masked:',cind,double(cind)/n_elements(sigma2d)
;
;
device,filename=fdir+'Tskytot_weighted_svd2_scut.ps'
plt_image, transpose(reform(radiores1mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores2mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores5mode(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores10mode(2,*,*))),/scalable,/colbar
device,/close
;
device,filename=fdir+'Tskytot_weighted_svd2_scut_weighted.ps'
plt_image, transpose(reform(radiores1mode(2,*,*)*weight(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores2mode(2,*,*)*weight(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores5mode(2,*,*)*weight(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(radiores10mode(2,*,*)*weight(2,*,*))),/scalable,/colbar
device,/close
;
device,filename=fdir+'Tskytot_weighted_svd2_scut_weight.ps'
plt_image, transpose(reform(weight(0,*,*))),/scalable,/colbar
plt_image, transpose(reform(weight(1,*,*))),/scalable,/colbar
plt_image, transpose(reform(weight(2,*,*))),/scalable,/colbar
plt_image, transpose(reform(weight(3,*,*))),/scalable,/colbar
device,/close
;
;
;
; put the holes back
Tsky(*,*,zlist)=0.

         ;
;         plt_image,v,/scalable,/colbar
;         plot,v(0,*)
;         plot,v(*,0)
;         plt_image,reform(v(0:10,*)),/scalable,/colbar
;         plot,w(0:10)
         ;
;plt_image, transpose(reform(radiosvd1mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiores1mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiosvd2mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiores2mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiosvd5mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiores5mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiosvd10mode(2,*,*))),/scalable,/colbar
;plt_image, transpose(reform(radiores10mode(2,*,*))),/scalable,/colbar
         ;
;zdevice,/close
         ;
         ;
         ;
openw,1,fdir+field+'.Tsky_svd_5mode'
printf,1,radiores5mode
close,1
;
openw,1,fdir+field+'.Tsky_svd_1mode'
printf,1,radiores1mode
close,1
;
openw,1,fdir+field+'.Tsky_svd_2mode'
printf,1,radiores2mode
close,1
;
openw,1,fdir+field+'.Tsky_svd_10mode'
printf,1,radiores10mode
close,1
;
openw,1,fdir+field+'.weight_svd_5mode'
printf,1,weight
close,1
;
;
;
;
end
