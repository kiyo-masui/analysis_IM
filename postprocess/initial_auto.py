""""This script post processes auto-correlations into final results."""

import os


file_root = os.getenv("GBT_YL")

autopower_root = file_root + "ps_result/power_auto_GBT_15hr_map_oldcal_"
transfers_root = file_root + "ps_result/bias/auto_GBT_15hr_map_oldcal_"

dir_middle = "_legendre_modes_0gwj_20"


