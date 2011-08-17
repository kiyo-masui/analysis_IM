pro combine_simcorrmode,field=field,mask=mask

if not keyword_set(field) then field='f4'

;
if field eq 'f3' then begin

; field3
openw,5,'field3.optsim.samew.mode.mask.corr.dat'
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f3',mode='0mode'
printf,5,'0mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f3',mode='1mode',mask=mask
printf,5,'1mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f3',mode='2mode',mask=mask
printf,5,'2mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f3',mode='3mode',mask=mask
printf,5,'3mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f3',mode='4mode',mask=mask
printf,5,'4mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
;correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f3',mode='5mode',mask=mask
;printf,5,'5mode'
;printf,5,'corr:',corrv1,corrv2
;printf,5,'error:',corre1,corre2
;printf,5,'zero:',corrm1,corrm2
;printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f3',mode='6mode',mask=mask
printf,5,'6mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
close,5
;
;
endif


;
;
if field eq 'f4' then begin

openw,5,'field4.optsim.samew.mode.mask.corr.dat'
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f4',mode='0mode'
printf,5,'0mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f4',mode='1mode',mask=mask          
printf,5,'1mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f4',mode='2mode',mask=mask
printf,5,'2mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f4',mode='3mode',mask=mask
printf,5,'3mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f4',mode='4mode',mask=mask
printf,5,'4mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f4',mode='5mode',mask=mask
printf,5,'5mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f4',mode='6mode',mask=mask
printf,5,'6mode'
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
;
close,5
;
;
endif



end
