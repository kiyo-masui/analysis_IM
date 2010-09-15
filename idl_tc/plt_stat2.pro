pro plt_stat2,covar2,mode1=mode1,mode2=mode2

;if not keyword_set(mode2) then mode2=mode1
mode1='10mode'
mode2='10mode'
mode=mode1+mode2

; plot the correlation function

nlag=61
nmid=(nlag-1)/2.
nf=2
nbin=6
;factor=0.000233280/0.000143717
;factor1=2./2.2
;factor2=2./3.5
;factor1=1.0964533
;factor2= 0.56113702
;factor1=0.94454721 ;10mode Oct21
;factor2=0.47637239
factor1=1.2681108 ;10modesvd3
factor2=0.61357939
;
;dir='Results/data/'
dir='Oct21/'
;
fn1=dir+'corrlag.radio.f3'+mode1+'svd3.dat'
fn2=dir+'corrlag.radio.f4'+mode2+'svd3.dat'
fn3=dir+'correrr.radio.f3'+mode1+'svd3.dat'
fn4=dir+'correrr.radio.f4'+mode2+'svd3.dat'
fn5=dir+'corrcovar.radio.f3'+mode1+'svd3.dat'
fn6=dir+'corrcovar.radio.f4'+mode2+'svd3.dat'
;
corrall=dblarr(nlag,nf)
corrall2=dblarr(nlag)
openr,1,fn1
readf,1,corrall2
close,1
corrall(*,0)=corrall2*factor1
openr,1,fn2
readf,1,corrall2
close,1
corrall(*,1)=corrall2*factor2
;
correrr=dblarr(nlag,nf)
correrr2=dblarr(nlag)
                                                                                                      
openr,1,fn3
readf,1,correrr2
close,1
correrr(*,0)=correrr2*factor1
openr,1,fn4
readf,1,correrr2
close,1
correrr(*,1)=correrr2*factor2
;
covar=dblarr(nlag,nlag,nf)
covar2=dblarr(nlag,nlag)
openr,1,fn5
readf,1,covar2
close,1
covar(*,*,0)=covar2*factor1^2.
openr,1,fn6
readf,1,covar2
close,1
covar(*,*,1)=covar2*factor2^2.
;
acorr=corrall
aerr=correrr^2.
;
;
corr=dblarr(nlag)
err=dblarr(nlag)
allcovar=dblarr(nlag,nlag)
for i=0,nlag-1 do begin
   corr(i)=total(acorr(i,*)/aerr(i,*))/total(1./aerr(i,*))
   err(i)=1./total(1./aerr(i,*))
   for j=0,nlag-1 do begin
      allcovar(i,j)=total(covar(i,j,*)/aerr(i,*)/aerr(j,*))/total(1./aerr(i,*)/aerr(j,*))
   endfor
endfor
;
print,'corr',corr
print,'err',err
;

lag=findgen(nlag)-nmid
;
set_plot,'ps'
device,filename='corr'+mode+'.ps'

; calculate the folded correlation fn -lag + lag
nmid=(nlag-1)/2
newcorr=dblarr(nmid+1)
newerr=newcorr
newcorr(0)=corr(nmid)
newerr(0)=err(nmid)
for i=1,nmid do begin
   newcorr(i)=(corr(nmid+i)+corr(nmid-i))/2.
   newerr(i)=allcovar(nmid+i,nmid+i)/4.+allcovar(nmid-i,nmid-i)/4.+$
             allcovar(nmid+i,nmid-i)/2.
endfor

print,'newcorr',newcorr
print,'newerr',newerr

; correlation binning
lag=findgen(nmid+1)*2+1e-20
xlag=alog(lag)
;take log bins
maxv=alog(max(lag))
nbin=6
;a0=sqrt(2.)
a0=1.
a=a0
factor=0.
for i=0,nbin-1 do begin
   factor=factor+a
   a=a*a0
endfor
bsize=maxv/factor
print,bsize,maxv,factor
print,a0^(findgen(nbin)+1)
binvalue=dblarr(nbin)
temp=0
for i=0,nbin-1 do begin
   temp=temp+bsize*a0^(i+1)
   binvalue(i)=temp
endfor

print,'binvalue',binvalue
;read,test
corr=dblarr(nbin)
covar=dblarr(nbin,nbin)
covar2=covar
err=dblarr(nbin)
for i=0,nbin-1 do begin
   if i eq 0 then ind=where(xlag le binvalue(i),cind) 
   if i gt 0 then ind=where(xlag le binvalue(i) and xlag gt binvalue(i-1),cind) 
   print,'i,cind',i,cind
   if cind gt 0 then begin
      corr(i)=mean(newcorr(ind))
      err(i)=total(newerr(ind))/double(cind)^2.
   endif
endfor

err=sqrt(err)

covar2=dblarr(nbin,nbin)
for i=0,nbin-1 do begin
   if i eq 0 then ind=where(xlag le binvalue(i),cind) 
   if i gt 0 then ind=where(xlag le binvalue(i) and xlag gt binvalue(i-1),cind) 
   if cind gt 0 then begin
      for j=0,nbin-1 do begin
         if j eq 0 then ind2=where(xlag le binvalue(j),cind2) 
         if j gt 0 then ind2=where(xlag le binvalue(j) and xlag gt binvalue(j-1),cind2) 
         if cind2 gt 0 then begin
            for ii=0,cind-1 do begin
               for jj=0,cind2-1 do begin
                  covar2(i,j)=covar2(i,j)+allcovar(ind(ii),ind2(jj))/double(cind)/double(cind2)
               endfor
            endfor
         endif
      endfor
   endif
endfor

for i=0,nbin-1 do begin
   if i eq 0 then ind=where(xlag le binvalue(i),cind) 
   if i gt 0 then ind=where(xlag le binvalue(i) and xlag gt binvalue(i-1),cind) 
   if cind gt 0 then begin
      a1=1./allcovar(ind,ind)/total(1./allcovar(ind,ind))
      for j=0,nbin-1 do begin
         if j eq 0 then ind2=where(xlag le binvalue(j),cind2) 
         if j gt 0 then ind2=where(xlag le binvalue(j) and xlag gt binvalue(j-1),cind2) 
         if cind2 gt 0 then begin
            a2=1./allcovar(ind2,ind2)/total(1./allcovar(ind2,ind2))
            for ii=0,cind-1 do $
               allcovar(ind(ii),ind2)=a1(ii)*allcovar(ind(ii),ind2)
            for jj=0,cind2-1 do begin
               allcovar(ind,ind2(jj))=a2(jj)*allcovar(ind,ind2(jj))
               ;print,'i,j',i,j,cind,cind2
               ;print,covar(i,j)
               ;print,allcovar(ind,ind2(jj))
               covar(i,j)=covar(i,j)+total(allcovar(ind,ind2(jj)))
            endfor
         endif
      endfor
   endif
endfor


err2=dblarr(nbin)
for i=0,nbin-1 do begin
   if i eq 0 then ind=where(xlag le binvalue(i),cind) 
   if i gt 0 then ind=where(xlag le binvalue(i) and xlag gt binvalue(i-1),cind) 
   if cind gt 0 then begin
      err2(i)=err2(i)+total(allcovar(ind,ind)/double(cind)^2.)
      if cind gt 1 then begin
         for ii=0,cind-1 do begin
            for jj=ii+1,cind-1 do begin
               err2(i)=err2(i)+allcovar(ii,jj)*2./double(cind)^2.
            endfor
         endfor
      endif
   endif
endfor

var2=covar2(indgen(nbin),indgen(nbin))
meancorr=total(corr/var2)/total(1./var2)

noise=0.
for i=0,nbin-1 do begin
   a1=1./var2(i)/total(1./var2)
   for j=0,nbin-1 do begin
      a2=1./var2(j)/total(1./var2)
      noise=noise+a1*a2*covar2(i,j)
   endfor
endfor

print,'meancorr',meancorr
print,'meanvar',noise,sqrt(noise)
print,'S/N',meancorr/sqrt(noise)

binvalue=exp(binvalue)
xlag=binvalue
xlag(0)=binvalue(0)/2.
for i=1,nbin-1 do begin
   xlag(i)=binvalue(i)-(binvalue(i)-binvalue(i-1))/2.
endfor
print,'binvalue',(binvalue)
print,'xlag',(xlag)
print,'xlag^2',xlag^2
err2=sqrt(err2)
ind=indgen(nbin)

; read in the model
fn='corrmodel.dat'
readcol,fn,x,model
;x=x*4.1
x=x*3.53
;x=x*6.
;model=model/((model(2)+model(3))/2.)*corr(0)

; simulation
;simv=0.000233280
simv=corr(0)
model=model/model(0)*simv
;model=model/model(2)*simv
help,x
help,model


covar3=covar2
for i=0,nbin-1 do begin
   for j=0,nbin-1 do begin
      covar3(i,j)=covar3(i,j)/sqrt(abs(covar2(i,i)))/sqrt(abs(covar2(j,j)))
      if j gt i then covar3(i,j)=0.
   endfor
endfor

bmodel=dblarr(nbin)
for i=0,nbin-1 do begin
   ind=where(abs(x-xlag(i)) eq min(abs(x-xlag(i))),cind)
   if cind eq 1 then $
      bmodel(i)=(model(ind))
endfor


covar2=covar2*(simv/corr(0))^2.
corr=corr*(simv/corr(0))
corr2=corr
covar4=covar2
for i=0,nbin-1 do begin
   corr(i)=corr(i)/bmodel(i)
   for j=0,nbin-1 do begin
      covar4(i,j)=covar4(i,j)/bmodel(i)/bmodel(j)
   endfor
endfor
;covar4(2,2)=1e20
;covar4(3,3)=1e20
;covar4(4,4)=1e20
;covar4(5,5)=1e20
;
covar5=covar4
LA_TRIRED, covar5, d, e,/double  
LA_TRIQL, d, e, covar5, /double

print,'eigenvalues:',d
print,'eigenvectors:',covar5

A=[[bmodel/bmodel]]
;A(2:*)=0.
print,'A',A
A=transpose(A)
help,A

print,'corr2',corr2
print,'covar2',covar2
;print,'icovar',icovar
print,'N^-1/2',sqrt(covar2)##invert(covar2)
print,'N^-1/2',invert(covar2)##sqrt(covar2)
print,'N^-1/2 b',sqrt(covar2)##invert(covar2)##transpose(corr2)

bg4=invert(transpose(A)##invert(covar4)##A)##transpose(A)##invert(covar4)##transpose(corr)
print,'bg4',bg4
A=[bmodel]
A=transpose(A)
bg2=invert(transpose(A)##invert(covar2)##A)##transpose(A)##invert(covar2)##transpose(corr2)

smatrix=((invert(transpose(A)##invert(covar2)##A)));*mse(0)
dx=sqrt(smatrix(indgen(1),indgen(1)))
help,smatrix
;fit2=xlag*bg2(1)+bg2(0)
;std2=xlag*dx(1)+dx(0)
fit2=model*bg2(0)
std2=model*dx(0)
print,'bg2',bg2
;print,'mse',mse
print,'s',smatrix
print,'dx',dx
print,'fit S/N',bg2/dx

;oplot,(xlag),xlag*bg(1)+bg(0)
;oplot,xlag,fit0,linestyle=1
;oploterr,xlag,fit,std,4
help,x,fit2
;oplot,xlag,fit2
;oploterr,xlag,fit2,std2,5

; calculate the residual of model and data
;dcorr=corr2-bg2#bmodel
;help,dcorr
covar6=covar4
LA_TRIRED, covar6, d, e,/double  
LA_TRIQL, d, e, covar6, /double
help,d
help,covar6
print,'eigenvalues:',d
print,'eigenvectors:',covar6
;print,covar6
print,'covar6',covar6(Nbin-1,*)
;eigenvector=covar6(4,*)
eigenvector=covar6(*,Nbin-1)
help,eigenvector
;print,'projection',eigenvector##transpose(eigenvector)##transpose(dcorr)
print,'projection',transpose(eigenvector)##(eigenvector)##transpose(corr)
print,'eigenvalues',d
;print,'projection',(transpose(eigenvector)##(eigenvector))
print,'corr',corr
print,'projection',(transpose(eigenvector)##(eigenvector)##transpose(corr))*corr2(0)
;eigenvector=covar6(Nbin-1,*)
;help,eigenvector
;print,'projection',((eigenvector)##transpose(eigenvector)##transpose(corr))*corr2(0)
;projection=eigenvector##transpose(eigenvector)##transpose(dcorr)
projection=(transpose(eigenvector)##(eigenvector)##transpose(corr))
projection=0
eigenvector2=covar6(*,Nbin-2)
projection2=(transpose(eigenvector2)##(eigenvector2)##transpose(corr))
projection2=0
corr_correct=(corr-projection-projection2)*bmodel
;corr_correct=(corr-projection)*bmodel
print,'corr_correct',corr_correct
;corr_correct=(corr-(transpose(eigenvector)##(eigenvector)##transpose(corr)))*corr2(0)

factor=1
;factor=2./3.5
plot,(xlag),corr_correct,/xlog,xtitle='!6 Lag [h!u-1!nMpc]',ytitle='!6 Correlation [K]',yran=[-5e-5,2e-4],$
     xran=[0.8,50],/xstyle,/ystyle,psym=3
err=reform(covar2(indgen(nbin),indgen(nbin)))
err=sqrt(err)
;print,'corr_correct',corr_correct
errplt=err
ind=where(err gt corr_correct,cind)
;if cind gt 0 then errplt(ind)=corr_correct(ind)-1e-20
oploterr,(xlag),corr_correct*factor,errplt*factor
;oploterr,(xlag(4:5)),corr_correct(4:5)+errplt(4:5),errplt(4:5)
oplot,x(1:*),fit2(1:*)
;oplot,xlag,(corr_correct+opterr)*factor,linestyle=2
;oplot,xlag,(corr_correct-opterr)*factor>1e-6,linestyle=2
;oplot,xlag*1.05,corr_correct,psym=4
;err=reform(covar2(indgen(nbin),indgen(nbin)))
;err=sqrt(err)
;print,'corr_correct',corr_correct
;errplt=err
;ind=where(err gt corr_correct,cind)
;if cind gt 0 then errplt(ind)=corr_correct(ind)-1e-20
;oploterr,xlag*1.05,corr_correct,errplt

;plt_image,covar3,colbar=[0.,0.2,0.4,0.6,0.8,1.0],/scalabl

device,/close

;print,binvalue,total(binvalue)
print,'bg',bg2
print,'dx',dx
print,'fit S/N',bg2/dx
print,'zero lag correlation',corr2(0)*bg2
print,'omega HI',corr2(0)*bg2/0.0003
print,'xlag',xlag
print,'corr',corr2
print,'err',err
print,'err2',sqrt(reform(covar(indgen(nbin),indgen(nbin))))
print,'covar2',covar2

openw,1,'corrdata'+mode+'.dat'

printf,1,'bg',bg2
printf,1,'dx',dx
printf,1,'fit S/N',bg2/dx
printf,1,'zero lag correlation',corr2(0)*bg2
printf,1,'omega HI',corr2(0)*bg2/0.00028
printf,1,'xlag',xlag
printf,1,'corr',corr2
printf,1,'err',err
printf,1,'covar2',covar2
close,1

;openw,1,'corrdata.opt.err.dat'
;printf,1,err
;close,1

end
