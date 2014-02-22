import h5py

def convert_trans(file_in, file_out):
    trans_in = h5py.File(file_in, "r")
    trans_out = h5py.File(file_out, "w")

    for treatment in trans_in:
        trans_cross = trans_in[treatment].value
        trans_out[treatment] = trans_cross ** 2.

    trans_in.close()
    trans_out.close()

#GBT_1hr_map_oldcalkiyo_x_WiggleZ_1hr_v2_order1_modetransfer.hd5
#GBT_1hr_map_oldcalkiyo_x_WiggleZ_1hr_v2_order1_beamtransfer.hd5
trans_root = "/mnt/raid-project/gmrt/eswitzer/GBT/bulksim/"

trans_15hr_beam_in = "GBT_15hr_map_oldcal_x_WiggleZ_15hr_blackman_order1_beamtransfer.hd5"
trans_15hr_mode_in = "GBT_15hr_map_oldcal_x_WiggleZ_15hr_blackman_order1_modetransfer.hd5"
trans_1hr_beam_in = "GBT_1hr_map_oldcal_x_WiggleZ_1hr_blackman_order1_beamtransfer.hd5"
trans_1hr_mode_in = "GBT_1hr_map_oldcal_x_WiggleZ_1hr_blackman_order1_modetransfer.hd5"

trans_15hr_beam_out = "GBT_15hr_map_oldcal_blackman_order1_one-sided_beamtransfer.hd5"
trans_15hr_mode_out = "GBT_15hr_map_oldcal_blackman_order1_one-sided_modetransfer.hd5"
trans_1hr_beam_out = "GBT_1hr_map_oldcal_blackman_order1_one-sided_beamtransfer.hd5"
trans_1hr_mode_out = "GBT_1hr_map_oldcal_blackman_order1_one-sided_modetransfer.hd5"

convert_trans(trans_root + trans_15hr_mode_in, trans_root + trans_15hr_mode_out)
convert_trans(trans_root + trans_15hr_beam_in, trans_root + trans_15hr_beam_out)
convert_trans(trans_root + trans_1hr_mode_in, trans_root + trans_1hr_mode_out)
convert_trans(trans_root + trans_1hr_beam_in, trans_root + trans_1hr_beam_out)
