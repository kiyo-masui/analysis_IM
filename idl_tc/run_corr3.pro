pro run_corr3,field=field


if field eq 'f3' then begin

correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f3',mode='5mode'
correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f3',mode='10mode'
correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f3',mode=''
correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f3',mode='20mode'
correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f3',mode='25mode'

endif

if field eq 'f4' then begin

correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f4',mode='5mode'
correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f4',mode='10mode'
correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f4',mode=''
correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f4',mode='20mode'
correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f4',mode='25mode'

endif


end
