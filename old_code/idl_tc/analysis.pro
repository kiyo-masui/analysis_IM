pro analysis

;mk_optmap,dname='02',signal='100muK'
;gen_radiosims2,/daisy,dname='02',signal='100muK'
;gen_radiosims2,/drift,dname='02',signal='100muK'
;docalib_gensim,/daisy,dname='02',signal='100muK'
;docalib_gensim,/drfit,dname='02',signal='100muK'

;docalib_svdfill_gain,/daisy,dname='02',/svdfg,/fitsfile
;docalib_svdfill_gain,/drift,dname='02',/svdfg,/fitsfile
docalib_svdfill_gain,/daisy,/sim,dname='02',signal='100muK',/svdfg
docalib_svdfill_gain,/drift,/sim,dname='02',signal='100muK',/svdfg

combine_corrmode2,/drift
combine_corrmode2,/daisy
combine_corrmode2,/drift,/sim
combine_corrmode2,/daisy,/sim
combine_corrmode2,/drift,/onlysim
combine_corrmode2,/daisy,/onlysim

end
