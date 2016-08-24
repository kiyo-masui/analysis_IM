from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import astropy
import astropy.units as u
import astropy.coordinates as ac
from astropy.time import Time
import datetime
#import scipy as sp
import scipy
import scipy.optimize

data_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cal_sdfits/B0043/'
#data_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cal_sdfits/1934/'

el_in = 90 - 29.1/60.
el_out = 90 - 50.8/60.
beam_rel_angles = {'1': [0,90],
                   '2': [120, el_in],
                   '3': [180, el_in],
                   '4': [240, el_in],
                   '5': [300, el_in],
                   '6': [0, el_in],
                   '7': [60, el_in],
                   '8': [90, el_out],
                   '9': [150, el_out],
                   '10': [210, el_out],
                   '11': [270, el_out],
                   '12': [330, el_out],
                   '13': [30, el_out]}

def load_data(dir, fname):
    hdu = fits.open(dir + fname)
    return hdu[1].data

def fits_list(dir):
    hdu_list = []
    for el in os.walk(dir):
        for fname in el[2]:
            hdu_list.append(load_data(dir, fname))
    return hdu_list

def combine_data(hdu_list, type):
    #CRVAL3 is RA.  CRVAL4 is dec
    array_list = []
    for data in hdu_list:
        array_list.append(data[type])
    return np.concatenate(array_list)

def make_hist(ra, dec, beam_bool, pix = 0.08, weight = None, polar = 'False', angle_fix = 0, ra_fix = 'True'):
    ra = ra[beam_bool]
    dec = dec[beam_bool]
    if polar == 'False':
        range = [np.max(ra) - np.min(ra), np.max(dec) - np.min(dec)]
        dec_bins = range[1]//pix
        if ra_fix == 'True':
            ra_bins = range[0]*np.cos(np.deg2rad(np.mean(dec) - angle_fix))//pix
        else:
            print 'No ra fix'
            ra_bins = range[0]//pix
        print str(range[0]) + ' in ' + str(ra_bins) + ' bins.' 
        if weight != None:
            weight = weight[beam_bool]
        info = np.histogram2d(ra, dec, bins = [ra_bins, dec_bins], weights=weight)
    else:
        x = np.multiply(dec, np.cos(ra))
        y = np.multiply(dec, np.sin(ra))
        range = [np.max(x) - np.min(x), np.max(y) - np.min(y)]
        dec_bins = range[1]//pix
        if ra_fix:
            ra_bins = range[0]*np.cos(np.deg2rad(np.mean(dec) - angle_fix))//pix
        else:
            ra_bins = range[0]//pix
        if weight != None:
            weight = weight[beam_bool]
        info = np.histogram2d(x, y, bins = [ra_bins, dec_bins], weights=weight)
    return info 

def plot_hist(info, filename, transpose= True, interpolation = 'none', vmax = 0):
    plt.clf()
    hist = info[0]
    if transpose:
        hist = np.transpose(hist)
    if vmax == 0:
        im = plt.imshow(hist, origin='lower', extent = [info[1][0], info[1][-1], info[2][0],info[2][-1]], interpolation = interpolation)
    else:
        im = plt.imshow(hist, origin='lower', extent = [info[1][0], info[1][-1], info[2][0],info[2][-1]], interpolation = interpolation,vmax=vmax)
    plt.colorbar()
    plt.savefig(filename)
    plt.clf()

def make_utc_times(times, dates):
    print dates.shape
    print times.shape
    final_times = [date_increment(dates[el],int((times[el]//3600)//24) ) + 'T' + str(int((times[el]//3600)%24)) + ':' + str(int((times[el]%3600)//60)) + ':' + str(times[el]%60) for el in range(times.shape[0])]
    return final_times

def date_increment(date, checker):
    if checker >= 1:
        date = datetime.datetime(int(date[0:4]),int(date[5:7]),int(date[8:10]),0,0,0)
        date += datetime.timedelta(days=checker)
        date = date.strftime("%Y-%m-%d")
    return date

def gauss_2d(grid, x_cen, y_cen, amp):
    sigma = 0.106
    x = grid[0]
    y = grid[1]
    #amp is the amplitude at the peak
    r = (x-x_cen)**2 + (y-y_cen)**2 
    val = r/(2*(sigma**2)*1.0)
    return amp*np.exp(-val)

def find_center(ra, dec, temp):
    #fitting to gaussian to find the center
    max = np.argmax(temp)
    ra_guess = ra[max]
    dec_guess = dec[max]
    print ra_guess
    print dec_guess
    print temp[max]
    x_vals = np.row_stack((ra,dec))
    print x_vals.shape
    print temp.shape
    #initial_guess = (ra_guess, dec_guess, 0.25, temp[max])
    initial_guess = (ra_guess, dec_guess, temp[max])
    print initial_guess
    bounds = ([ra_guess - 0.25, dec_guess - 0.25, 0.125, temp[max]/1.4],[ra_guess + 0.25, dec_guess + 0.25, 0.375, temp[max]*1.4])
    popt, pcov = scipy.optimize.curve_fit(gauss_2d, x_vals, temp, p0=initial_guess)
    print x_vals[0][max]
    print x_vals[1][max]
    return popt

def gauss(x,y,x_cen,y_cen, sigma, amp):
    r = (x-x_cen)**2 + (y-y_cen)**2
    val = r/(2*(sigma**2)*1.0)
    return amp*np.exp(-val)

def fit_lsq(az, el, temp):
    #fit to constant plus gaussian
    def model(pars):
        #pars[0] is background Tsys guess, pars[1] is cos(el)daz cent
        #pars[2] is del center, pars[3] is gaussian sigma,
        #pars[4] is gaussian amplitude at max
        model = pars[0] + gauss(az,el,pars[1],pars[2],pars[3],pars[4])
        return model
    max_temp = np.argmax(temp)
    init_pars = [np.min(temp), az[max_temp], el[max_temp], 0.106, np.max(temp) - np.min(temp)]
    fit_pars, err = scipy.optimize.leastsq(lambda q: model(q) - temp, init_pars, ftol=1e-5, xtol=1e-5)
    return fit_pars

def azel_to_xyz(az, el):
    x = np.cos(np.deg2rad(el))*np.cos(np.deg2rad(az))
    y = np.cos(np.deg2rad(el))*np.sin(np.deg2rad(az))
    z = np.sin(np.deg2rad(el))
    #ans = np.array([[el[0],el[1],el[2]] for el in zip(x,y,z)])
    ans = np.array([x,y,z])
    return ans

def xyz_to_azel(b):
    x = b[0]
    y=b[1]
    z=b[2]
    az = np.rad2deg(np.arctan(1.0*y/x))
    el = np.rad2deg(np.arcsin(z))
    if abs(x) > 10**-12 and x<0:
        az += 180
    if x>0 and y<0 and abs(x) > 10**-12 and abs(y) > 10**-12:
        az += 360
    return az,el

def xy_to_theta(x,y):
    theta = np.rad2deg(np.arctan(1.0*y/x))
    if abs(x) > 10**-12 and x<0:
        theta += 180
    if x>0 and y<0 and abs(x) > 10**-12 and abs(y) > 10**-12:
        theta += 360
    return theta

def rotate_xy(x,y,theta):
    angle = np.deg2rad(theta)
    rot = [np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang), np.cos(ang)]]) for ang in angle]
    ans = np.array([np.dot(np.array(el[2]), np.array([el[0],el[1]])) for el in zip(x,y,rot)])
    return ans

def rotate_abt_x(b,theta):
    angle = np.deg2rad(theta)
    rot = [np.array([[1,0,0],[0,np.cos(ang),-np.sin(ang)],[0,np.sin(ang), np.cos(ang)]]) for ang in angle]
    ans = np.array([np.dot(el[0],el[1]) for el in zip(rot,b)])
    return ans

def rotate_abt_z(b,theta):
    angle = np.deg2rad(theta)
    rot = [np.array([[np.cos(ang),-np.sin(ang),0],[np.sin(ang), np.cos(ang),0], [0,0,1]]) for ang in angle]
    ans = np.array([np.dot(el[0],el[1]) for el in zip(rot,b)])
    return ans

def rotate_abt_y(b,theta):
    angle = np.deg2rad(theta)
    rot = [np.array([[np.cos(ang),0,np.sin(ang)],[0,1,0],[-np.sin(ang),0, np.cos(ang)]]) for ang in angle]
    ans = np.array([np.dot(el[0],el[1]) for el in zip(rot,b)])
    return ans

def az_el_beam(az, el, beams):
    #finds apparent az el of a beam
    #Try adding 15 degrees for receiver angle?
    xyz = np.array([azel_to_xyz(beam_rel_angles[str(beam)][0]+15,beam_rel_angles[str(beam)][1]) for beam in beams])
    xyz = rotate_abt_y(xyz, 90 - el)
    xyz = rotate_abt_z(xyz, az)
    print xyz
    ans = np.array([xyz_to_azel(el) for el in xyz])
    return ans

def angle_av_profile(cent, bins, rmax, ra, dec, vals):
    ra -= cent[0]
    dec -= cent[1]
    r = ((1.0*ra)**2 + (1.0*dec)**2)**0.5
    range = (0,rmax)
    hits = np.histogram(r, bins=bins, range=range)
    prof = np.histogram(r, bins=bins, range=range, weights=vals)
    prof = [prof[0]/hits[0], prof[1]]
    return prof, hits

if __name__ == '__main__':
    #dir_out = '/scratch2/p/pen/andersoc/second_parkes_pipe/cal_sdfits/B0043_hists/'
    #dir_out = '/scratch2/p/pen/andersoc/second_parkes_pipe/cal_sdfits/1934_hists/'
    dir_out = '/scratch2/p/pen/andersoc/second_parkes_pipe/cal_sdfits/B0043_angle_avgs/'
    #dir_out = '/scratch2/p/pen/andersoc/second_parkes_pipe/cal_sdfits/B0043_hists/beam_adjust/'
    hdu_list = fits_list(data_dir)
    ra = combine_data(hdu_list, 'CRVAL3')
    dec = combine_data(hdu_list, 'CRVAL4')
    beams = combine_data(hdu_list, 'BEAM')
    data = combine_data(hdu_list, 'DATA')
    fav_data = np.mean(data, axis = 4)
    xx_dat = fav_data[:,0,0,0]
    i_dat = np.mean(fav_data, axis = (1,2,3))

    #For making daz del maps.
    times = combine_data(hdu_list, 'TIME')
    dates = combine_data(hdu_list, 'DATE-OBS')
    utc_times = make_utc_times(times, dates)
    print utc_times[-200:] 
    print type(utc_times)
    utc_times = Time(utc_times, format = 'isot', scale= 'utc')
    #Parkes location

    #According to pkscat90, B0043-424 is at RA 43 min, 54.3 sec, dec -42 deg, 24min  
    source_ra = 15*(43.0 + (54.3/60.0))/60.0
    source_dec = -42.0 - 24.0/60.0

    observatory = ac.EarthLocation(lat=-32.9983*u.deg, lon=148.2636*u.deg, height=64*u.m)
    print utc_times
    frame = ac.AltAz(obstime=utc_times,location=observatory)
    beam1bool = beams == 1
    ra1 = ra[beam1bool]
    dec1 = dec[beam1bool]
    i_dat1 = i_dat[beam1bool]
    #Getting coordinates of brightest point in beam 1
    radecs = ac.SkyCoord(ra = ra1[np.argmax(i_dat1)]*u.deg, dec = dec1[np.argmax(i_dat1)]*u.deg)
    #radecs = ac.SkyCoord(ra = source_ra*u.deg, dec = source_dec*u.deg)
    print radecs
    altaz_source = radecs.transform_to(frame)
    print altaz_source
    alt_source = np.array(altaz_source.alt.deg)
    az_source = np.array(altaz_source.az.deg)
    ra = combine_data(hdu_list, 'AZIMUTH')
    dec = combine_data(hdu_list, 'ELEVATIO')

    #Attempt to calculate az and el that each beam is effectively looking at
    radec = az_el_beam(ra, dec, beams)
    az_diff = ra - radec[:,0]
    el_diff = dec - radec[:,1]
    ra = radec[:,0]
    dec = radec[:,1]
    #Now ra and dec are az el.  Next, make them daz, del.
    #Actually, let's also try multiplying daz by cos(el), so both are on equal footing
    ra -= az_source
    ra *= np.cos(np.deg2rad(dec))
    az_diff *= np.cos(np.deg2rad(dec))
    dec -= alt_source
    print 'az range ' + str(np.max(ra)) + ' ' + str(np.min(ra))
    print 'alt range ' + str(np.max(alt_source)) + ' ' + str(np.min(alt_source))

    #Rotate cos(el)daz and del such that direction to middle beam (1) is Right
    angle = np.array([xy_to_theta(el[0],el[1]) for el in zip(az_diff,el_diff)])
    angle *= -1.0
    print angle
    print zip(az_diff,el_diff)
    adj = rotate_xy(ra,dec,angle)
    ra = adj[:,0]
    dec = adj[:,1]

    for beam in range(1,14):
        bool = beams == beam

        #Code below is just to line up az, el of all beams, if desired.
        az_off = ra[bool][np.argmax(i_dat[bool])]
        alt_off = dec[bool][np.argmax(i_dat[bool])]
        print 'Beam ' + str(beam) + ' offset by' + str(az_off) + ', ' + str(alt_off)
        #ra[bool] -= az_off
        #dec[bool] -= alt_off

        #Optionally, also subtract off min intensity from each beam
        print np.min(i_dat[bool])
        print np.max(i_dat[bool])
        i_dat[bool] -= np.min(i_dat[bool])
  
        print 'Fit parameters are ' 
        #cent = find_center(ra[bool],dec[bool],i_dat[bool])
        params = fit_lsq(ra[bool], dec[bool], i_dat[bool])
        cent = [params[1],params[2]]
        print params
        r = 1
        bins = 30
        dr = r/(bins*1.0)
        filename = dir_out + 'beam' + str(beam) + 'avg_lsqcent'
        prof, hist = angle_av_profile(cent[0:2], bins, r, ra[bool], dec[bool], i_dat[bool]) 
        plt.clf()
        plt.plot(np.arange(dr/2.0,r,dr), prof[0])
        plt.savefig(filename)
        

        xx_tot = make_hist(ra, dec, bool, pix= 0.04, weight = i_dat, angle_fix = 75, ra_fix = 'False')
        xx_hits = make_hist(ra, dec, bool, pix = 0.04, angle_fix = 75, ra_fix = 'False')
        xx_beam = [xx_tot[0]/xx_hits[0], xx_tot[1], xx_tot[2]]
        #xx_beam = xx_tot
        print str(beam)
        print xx_beam[0].shape
        print xx_beam[0]
        max = np.unravel_index(np.nanargmax(xx_beam[0]), xx_beam[0].shape)
        print max
        max_pos = (np.min(xx_tot[1]) + (max[0]/xx_beam[0].shape[0])*(np.max(xx_tot[1])-np.min(xx_tot[1])),np.min(xx_tot[2]) + (max[1]/xx_beam[0].shape[0])*(np.max(xx_tot[2])-np.min(xx_tot[2])))
        print max_pos
        print np.nanmax(xx_beam[0])
        print [xx_beam[1][0],xx_beam[1][-1],xx_beam[2][0],xx_beam[2][-1]]
        #plot_hist(xx_beam, dir_out + 'beam' + str(beam) + '_cos(el)daz_del_nearest_' + 'I' + '_rot', interpolation = 'nearest')
'''        if beam == 1:
            np.save(dir_out +'data/'+ 'beam' + str(beam) + '_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_hist',xx_beam[0])
            np.save(dir_out +'data/'+ 'beam' + str(beam) + '_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_hitmap',xx_hits[0])
            np.save(dir_out +'data/' +'beam' + str(beam) + '_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_xedges',xx_beam[1])
            np.save(dir_out +'data/' + 'beam' + str(beam) + '_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_yedges',xx_beam[2])
   
    #Combine beams 2 through 7 and beams 8 through 13, excluding 10 and 13
    inner_bool = np.logical_and(beams>1, beams<8)
    outer_bool = np.logical_and(np.logical_and(beams>7, beams<13), beams != 10)
    xx_tot = make_hist(ra, dec, inner_bool, pix=0.04, weight = i_dat, angle_fix = 75, ra_fix = 'False')
    xx_hits = make_hist(ra, dec, inner_bool, pix=0.04, angle_fix = 75, ra_fix = 'False')
    xx_beam = [xx_tot[0]/xx_hits[0], xx_tot[1], xx_tot[2]]
    plot_hist(xx_beam, dir_out + 'beam2through7_cos(el)daz_del_nearest_I_stacked_rot_vmax3', interpolation = 'nearest',vmax=3)
    plot_hist(xx_beam, dir_out + 'data/' + 'beam2through7_cos(el)daz_del_nearest_' + 'I' + '_rot', interpolation = 'nearest')
    np.save(dir_out +'data/' + 'beam2through7_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_hist',xx_beam[0])
    np.save(dir_out +'data/' + 'beam2through7_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_hitmap',xx_hits[0])
    np.save(dir_out +'/data/' + 'beam2through7_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_xedges',xx_beam[1])
    np.save(dir_out +'/data/' + 'beam2through7_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_yedges',xx_beam[2])

    xx_tot = make_hist(ra, dec, outer_bool, pix=0.04, weight = i_dat, angle_fix = 75, ra_fix = 'False')
    xx_hits = make_hist(ra, dec, outer_bool, pix=0.04, angle_fix = 75, ra_fix = 'False')
    xx_beam = [xx_tot[0]/xx_hits[0], xx_tot[1], xx_tot[2]]
    plot_hist(xx_beam, dir_out + 'beam8_9_11_12_cos(el)daz_del_nearest_I_stacked_rot_vmax3', interpolation = 'nearest',vmax=3)
    np.save(dir_out +'data/' + 'beam8_9_11_12_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_hist',xx_beam[0])
    np.save(dir_out +'data/' + 'beam8_9_11_12_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_hitmap',xx_hits[0])
    np.save(dir_out +'/data/' + 'beam8_9_11_12_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_xedges',xx_beam[1])
    np.save(dir_out +'/data/' + 'beam8_9_11_12_cos(el)daz_del_nearest_' + 'I' + '_rot' + '_yedges',xx_beam[2])
'''
''' 
    #Beam 10 and 13 look kind of bad, calculate map from all other beams.
    bool = np.logical_and(beams != 10, beams != 13)
    xx_tot = make_hist(ra, dec, bool, pix=0.04, weight = i_dat, angle_fix = 75, ra_fix = 'False')
    xx_hits = make_hist(ra, dec, bool, pix=0.04, angle_fix = 75, ra_fix = 'False')
    xx_beam = [xx_tot[0]/xx_hits[0], xx_tot[1], xx_tot[2]]
    max = np.unravel_index(np.nanargmax(xx_beam[0]), xx_beam[0].shape)
    plot_hist(xx_beam, dir_out + 'all_except_beam10and13_cos(el)daz_del_nearest_I_stacked', interpolation = 'nearest')
    print combine_data(hdu_list, 'PARANGLE')
'''
