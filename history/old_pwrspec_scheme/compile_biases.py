import copy
from correlate import compile_single_crosspwr as csc
from correlate import pwrspec_estimation as pe

cleaned_fgcorr_ini = "input/ers/spectra/GBT_15hr_map_fdgcal_cleanedplussim_cleaned_noconv_combined_xspecsimsingle-nocomp.ini"
native_fgcorr_ini = "input/ers/spectra/GBT_15hr_map_fdgcal_cleaned_noconv_combined_xspecsimsingle-nocomp.ini"
cleaned_fgcorr_all = csc.wrap_batch_single_crosspwr(cleaned_fgcorr_ini)
native_fgcorr_all = csc.wrap_batch_single_crosspwr(native_fgcorr_ini)

for treatment in cleaned_fgcorr_all:
    cleaned_fgcorr = cleaned_fgcorr_all[treatment]
    native_fgcorr = native_fgcorr_all[treatment]

    pwrdiff = copy.deepcopy(cleaned_fgcorr)
    pwrdiff[0]['binavg'] = cleaned_fgcorr[0]['binavg'] + native_fgcorr[0]['binavg']
    pwrdiff[1]['binavg'] = cleaned_fgcorr[1]['binavg'] + native_fgcorr[1]['binavg']
    pwrdiff[2]['binavg'] = cleaned_fgcorr[2]['binavg'] + native_fgcorr[2]['binavg']

    mtag = "correlated_fg_%s" % treatment
    pe.summarize_pwrspec(pwrdiff[0], pwrdiff[1], pwrdiff[2],
                         mtag, outdir="./plots/correlated_fg/")
