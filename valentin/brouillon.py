import numpy as np
import matplotlib.pyplot as plt

import core.algebra as al


file_name = "/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_15hr_map_autopaper/sec_C_cleaned_clean_map_I_with_B_20modes.npy"

map = al.load(file_name)
map = al.make_vect(map)

#---------------------------------------------------------------------
#DFT : versions naives, complexite quadratique
def iexp(n):
	return complex(np.cos(n), np.sin(n))

def dft(vec):
	n = len(vec)
	res = np.zeros(n)
	for k in range(n):
		for i in range(n):
			res[k] += vec[i] * iexp(-2*np.pi*k*i/n)
	return res	

def idft(vec):
	n = len(vec)
	res = np.zeros(n)
	for i in range(n):
		for k in range(n):
			res[i] += vec[k] * iexp(2*np.pi*i*k/n)
	return res/n
	
def deriv(vec):
	indx = np.array(range(len(vec)))
	return idft(1j*indx*dft(vec))

#---------------------------------------------------------------------
#Essais de fenetre gaussienne	
def gaussian_window(map, sigma=1.):
	
	return np.exp(-norm(x)**2/(2*sigma*sigma))
	
def part_deriv_gaussian(x, i, sigma=1.):
	"""partial derivative with respect to i of th gaussian window
	i is the index of the direction : 
	i=0 (resp. 1) for  x (resp. y) direction"""
	return -x[i]/sigma*sigma * gaussian_window(x,sigma)
	
#Verification de l'impact de sigma
def show_sigma(map, sigma, i = 101):
	ws = fast_gaussian_window(map, sigma)
	map_s = fftconvolve(map, ws, 'same')
	plt.figure(i)
	plt.imshow(map_s[30,:,:].T)
	del ws
	del map_s

#---------------------------------------------------------------------
#Calcul du kernel
def kernel(k,i,j):
	k_2 = k[0]**2 + k[1]**2 + k[2]**2
	delta = 1 if i==j else 0
	return k_2/(k[0]**2 + k[1]**2)*(k[i]*k[j]/k_2 - delta/2.)
	
def large_map_kernel(map):
	t = time.time()
	l,m,n = map.shape
	ker_xx = np.zeros((2*l,2*m,2*n))
	ker_yy = np.zeros((l,m,n))
	ker_xy = np.zeros((l,m,n))
	for k in range(l):
		for i in range(m):
			for j in range(n):
				cst2 = i*i + j*j + 0.0
				if (cst2 != 0):
					cst = (cst2 + k*k)/cst2
					ker_xx[k,i,j] = cst*(i*i/cst2 - 0.5)
					ker_xx[k,2*m-i-1,j] = ker_xx[k,i,j]
					ker_xx[k,i,2*n-j-1] = ker_xx[k,i,j]
					ker_xx[k,2*m-i-1,2*n-j-1] = ker_xx[k,i,j]
					ker_xx[2*l-k-1,i,j] = cst*(i*i/cst2 - 0.5)
					ker_xx[2*l-k-1,2*m-i-1,j] = ker_xx[k,i,j]
					ker_xx[2*l-k-1,i,2*n-j-1] = ker_xx[k,i,j]
					ker_xx[2*l-k-1,2*m-i-1,2*n-j-1] = ker_xx[k,i,j]
					ker_yy[k,i,j] = cst*(j*j/cst2 - 0.5)
					ker_xy[k,i,j] = cst*(i*j/cst2)
#	ker_xx = np.fft.ifftn(ker_xx)
#	ker_yy = np.fft.ifftn(ker_yy)
#	ker_xy = np.fft.ifftn(ker_xy)
	t = time.time() - t
	print t
	return ker_xx,ker_yy,ker_xy

def density_field(t_set, k_set):
	t = time.time()
	kappa_xx = fftconvolve(t_set[0], k_set[0], 'same')
	kappa_yy = fftconvolve(t_set[1], k_set[1], 'same')
	kappa_xy = fftconvolve(t_set[2], k_set[2], 'same')
	kappa = kappa_xx + kappa_yy + 2*kappa_xy
	t = time.time() - t
	print t
	return kappa
	
tmp = np.zeros((2*l, 2*m, 2*n))
tmp[:l,:m,:n] = t_xx
l_t_xx = np.zeros((2*l-1, 2*m-1, 2*n-1))
l_t_xx[:] = tmp[:2*l-1, :2*m-1, :2*n-1]


def map_kernel(map):
	t = time.time()
	l,m,n = map.shape
	ker_xx = np.zeros((l,m,n))
	ker_yy = np.zeros((l,m,n))
	ker_xy = np.zeros((l,m,n))
	
	for k in range(l):
		for i in range(m):
			for j in range(n):
				x = (i - m/2) #/((m-1)*ra_delta)
				y = (j - n/2) #/((n-1)*dec_delta)
				z = (k - l/2) /((l-1)*freq_delta)
				cst2 = x*x + y*y + 0.0
				if (cst2 != 0):
					cst = (cst2 + z*z)/cst2
					ker_xx[k,i,j] = cst*(x*x/cst2 - 0.5)
					ker_yy[k,i,j] = cst*(y*y/cst2 - 0.5)
					ker_xy[k,i,j] = cst*(x*y/cst2)
					
	#Shift so that the fonction ifftn work
	ker_xx = np.fft.ifftshift(ker_xx)
	ker_yy = np.fft.ifftshift(ker_yy)
	ker_xy = np.fft.ifftshift(ker_xy)

#	#Apply the inverse dft
#	ker_xx = np.fft.ifftn(ker_xx)
#	ker_yy = np.fft.ifftn(ker_yy)
#	ker_xy = np.fft.ifftn(ker_xy)
#
#	#Keep only the real part (imaginary part should be negligible)
#	ker_xx = ker_xx.real
#	ker_yy = ker_yy.real
#	ker_xy = ker_xy.real
#
#	#Shift again so that x=0 is at the centre of the matrix
#	ker_xx = np.fft.fftshift(ker_xx)
#	ker_yy = np.fft.fftshift(ker_yy)
#	ker_xy = np.fft.fftshift(ker_xy)
	
	t = time.time() - t
	print t
	return ker_xx,ker_yy,ker_xy

	
k_xx = k_set[0]
t_xx = t_set[0]

if_t_xx = np.fft.ifftshift(t_xx)
if_t_xx = np.fft.fftn(if_t_xx)
if_kap_xx = if_t_xx * k_xx
if_kap_xx = np.fft.ifftn(if_kap_xx)
if_kap_xx = if_kap_xx.real
if_kap_xx = np.fft.fftshift(if_kap_xx)

f_k_xx = np.fft.ifftn(k_xx)
f_k_xx = f_k_xx.real
f_k_xx = np.fft.fftshift(f_k_xx)
#f_kap_xx = fftconvolve(t_xx, f_k_xx,'same')

ff_kap_xx = fftconvolve(t_xx, f_k_xx,'full')
for k in range(l):
    for i in range(m):
        for j in range(n):
            test[k,i,j] = f_kap_xx[(l/2 + k), (m/2 + i), (n/2 + j)]


            
def save_map(map, name):
	file_name = "valentin/" + name
	file_w = open(file_name,'w')
	l,m,n = map.shape
	file_w.write("{} {} {} ".format(l,m,n))
	for k in range(l):
		for i in range(m):
			for j in range(n):
				file_w.write("{} ".format(map[k,i,j]))
	file_w.close()


offset = 10
nb_maps = 60
path = 'valentin/random_maps/'
name = 'kappa_gaussian_'	
for i in range(nb_maps):
	print 'Step',i+1,':'
	kappa = full_process()[1]
	np.save(path + name + "%s" %(offset+i+1),kappa)
	print 'Step',i+1,'complete'
	print '-----------------------------------------------'
offset += nb_maps

nb_maps = offset
l,m,n = np.load(path + 'test.npy').shape
maps_set = np.zeros((nb_maps,l,m,n))
for i in range(nb_maps):
	maps_set[i] = np.load(path + name + "%s" %(i+1) + '.npy')



#for i in range(nb_maps):
#	tmp = np.load(path + name + "%s" %(i+1) + '.npy')
#	tmp = (4.*np.sqrt(2*np.pi))**(-6) * tmp
#	np.save(path + name + "%s" %(i+1),tmp)




#----------------------------------------
#Useful stuff
#----------------------------------------

#animation to show every slices z
def show_seq(mat,step=0.2, ref=0):
	n = mat.shape[0]
	fig, ax = plt.subplots()

	im = ax.imshow(mat[ref,:,:].T)
	fig.show()
	for i in range(n):
		im.set_data(mat[i,:,:].T)
		fig.canvas.draw()
		time.sleep(.2)

#compact a map by a factor k
def low_resolution(map, k=2):
	l,m,n = map.shape
	map_red = map.reshape(l/k,k,m/k,k,n/k,k)
	map_red = np.sum(map_red, axis=5)
	map_red = np.sum(map_red, axis=3)
	map_red = np.sum(map_red, axis=1)
	return map_red/k**3

#short vesion of density field calculation (not tested yet)	
def short_density_field(t_set,k_set):
	t = time.time()
	f_t = np.fft.fftn(t_set,axes=[1,2,3])
	f_t = f_t * k_set
	kappa = f_t[0] + f_t[1] + 2*f_t[2]
	kappa = np.fft.ifftn(kappa).real
	t = time.time() - t
	print t
	return kappa


#Attempt to create a fast version of kernel calculation, avoiding triple for loop
def fast_kernel(k_axes):
	kz = k_axes[0]
	kx = k_axes[1]
	ky = k_axes[2]
	k_xx = (kx*kx)[:,None] * np.ones((1,len(ky)))
	k_yy = np.ones((len(kx),1)) * (ky*ky)[None,:]
	k_xy = kx[:,None] * ky[None,:]
	rad_2d = k_xx + k_yy #calculate kx^2 + ky^2
	k_xx = k_xx/rad_2d - 0.5
	k_yy = k_yy/rad_2d - 0.5
	k_xy = k_xy/rad_2d
	k_z_3d = (kz*kz)[:,None,None] * (1/rad_2d)[None,:,:] #a tester!!
	k_xx = (1 + k_z_3d)*k_xx[None,:,:]
	k_yy = (1 + k_z_3d)*k_yy[None,:,:]
	k_xy = (1 + k_z_3d)*k_xy[None,:,:]
	k_xx[:,len(kx)/2,len(ky)/2] = 0. #Temporary solution, ok only if x_centre,y_centre = 0
	k_yy[:,len(kx)/2,len(ky)/2] = 0.
	k_xy[:,len(kx)/2,len(ky)/2] = 0.
		#Shift so that the fonction ifftn work
	k_xx = np.fft.ifftshift(k_xx)
	k_yy = np.fft.ifftshift(k_yy)
	k_xy = np.fft.ifftshift(k_xy)
	return k_xx,k_yy,k_xy
	
def kernel_compare():
	map = sim_random_map((128,128,128),(100.,100.,100.))
	map_set = MapBuilder(map, full=False)
	t = time.time()
	kz = map_set.get_k_axis('z')
	kx = map_set.get_k_axis('x')
	ky = map_set.get_k_axis('y')
	k2 = fast_kernel([kz,kx,ky])
	t = time.time() - t
	print t

#after comparison, the fast version is really faster (more than 100x)
	

#Two ideas for implementation of gaussian window calculation
def gaussian_window(map, sigma):
	z = map.get_axis('z')
	x = map.get_axis('x')
	y = map.get_axis('y')
	mesh = np.indices(map.shape)
	sph = np.zeros(mesh.shape)
	sph[0] = z[mesh[0]] - map.info['z_centre']
	sph[1] = x[mesh[1]] - map.info['x_centre']
	sph[2] = y[mesh[2]] - map.info['y_centre']
	rad = sph[0]**2 + sph[1]**2 + sph[2]**2
	return (sigma*np.sqrt(2*np.pi))**(-3) * np.exp(-rad/(2.*sigma*sigma))	
	
def alt_gaussian(map,sigma):
	l,m,n = map.shape
	z = map.get_axis('z')
	x = map.get_axis('x')
	y = map.get_axis('y')
	rad_3d = (z*z)[:,None,None] * np.ones((1,m,n))
	rad_3d += (x*x)[None,:,None] * np.ones((l,1,n))
	rad_3d += (y*y)[None,None,:] * np.ones((l,m,1))
	return (sigma*np.sqrt(2*np.pi))**(-3) * np.exp(-rad_3d/(2.*sigma*sigma))

def gaus_compare():
	map = sim_random_map((512,512,512),(100.,100.,100.))
	t = time.time()
	gaussian_window(map, 8.)
	t = time.time() - t
	print 'v1 : %.3f sec' %t
	t = time.time()
	alt_gaussian(map, 8.)
	t = time.time() - t
	print 'v2 : %.3f sec' %t
	
#after comparison, alt_version is faster (about 2x)



#----------------------------------------

def gaussian_noise(n,px=64):
	rd_set = []
	for i in range(n):
		print "----------------------"
		print "Simulation",i
		t = time.time()
		rdm = sim_random_map((px,px,px),(322.,322.,322.))
		rds = MapBuilder(rdm, full=False, verbose=False, log_smooth=True)
		f_rdk = np.fft.fftn(rds.kappa)
		rd_set.append(f_rdk)
		t = time.time() - t
		print t
		print " "
	return rd_set

def pix_statistics(map_set):
	n = len(map_set)
	x,y,z = map_set[0].shape
	std = (1.+0.j)*np.zeros((map_set[0].shape))
	mn = (1.+0.j)*np.zeros((map_set[0].shape))
	for i in range(n):
		map = map_set[i]
		mn += map
		std += map*map.conjugate()
	mn = mn/(x*y*z)
	std = std/(x*y*z)
	return mn,np.sqrt(std.real)
	
def show_hist(std,rad,z,nb_bins=100):
	a = np.histogram(rad, nb_bins, weights=std)[0]
	b,bins = np.histogram(rad,nb_bins)
	axis = (2*np.pi/322.)*np.arange(nb_bins)*128/nb_bins
	plt.plot(axis,(a/b))
	return a/b, bins
	
def r_current(kap,map):
	m1m1 = np.sum(kap*kap)
	m2m2 = np.sum(map*map)
	m1m2 = np.sum(kap*map)
	return m1m2/np.sqrt(m1m1*m2m2)

#A tester + apporter modif necessaires
def fast_quad_est(maps):
	l,m,n = maps.info['shape']
	sigma_k = maps.info['sigma_recon']**(-1)
	kz = maps.get_k_axis('z')
	kx = maps.get_k_axis('x')
	ky = maps.get_k_axis('y')
	w = maps.gaussian_window(sigma_k, axes=[kz,kx,ky])
	
	d_x = 1.j*kx[None,:,None] * w
	d_y = 1.j*ky[None,None,:] * w
	d_x = np.fft.ifftshift(d_x)
	d_y = np.fft.ifftshift(d_y)

	f_map = np.fft.ifftshift(maps.map)
	f_map = np.fft.fftn(f_map)
	d_x = d_x * f_map
	d_y = d_y * f_map
	d_x = np.fft.ifftn(d_x).real
	d_y = np.fft.ifftn(d_y).real
	d_x = np.fft.fftshift(d_x)
	d_y = np.fft.fftshift(d_y)

	map_log = f_map*np.fft.ifftshift(w)
	map_log = np.fft.ifftn(map_log)
	map_log = np.fft.fftshift(map_log.real)
	d_x = d_x/(1. + map_log)
	d_y = d_y/(1. + map_log)
	
	return d_x*d_x,d_y*d_y,d_x*d_y


def complex_mult():
	a = np.random.normal(size=(512,512,512),scale=10.)
	b = 1.j*np.random.normal(size=(512,512,512),scale=10.)
	c = np.random.normal(size=(512,512,512),scale=10.)
	d = a + b
	t = time.time()
	e = a*c
	t1 = time.time() - t
	t = time.time()
	e = d*c	
	t2 = time.time() - t
	t = time.time()
	e = d*b	
	t3 = time.time() - t
	print "Real x Real:",t1
	print "Complex x Real:",t2
	print "Complex x Complex:",t3
	

def cut(map, epsilon=1e-5):
	valid = map>=(-1+epsilon)
	return valid*map - (1-epsilon)*(1-valid)
	
def rescale(map, epsilon=1e-5):
	return map*(-1+epsilon)/map.min()

def log_smoothing(delta, w):
	t = time.time()
	delta1 = rescale(delta)
	del delta
	delta2 = own_fftconvolve(delta1,w)
	del delta1,w
	res = np.log(1+delta2)
	del delta2
	t = time.time() - t
	print t
	return res

	
	
def iter_mask(n_slice, px_ref=(635,88)):
	mn = n_slice[px_ref]
	mask0 = np.zeros_like(n_slice)
	mask1 = n_slice>(mn/2)
	diff = mask1 - mask0
	c=0
	while (diff.max()!=0)&(c<20):
		mn = np.sum(n_slice*mask1)/np.sum(mask1)
		mask0 = np.copy(mask1)
		mask1 = n_slice>(mn/2)
		diff = mask1 - mask0
		c+=1
	if c<20:	
#		print "converge in", c, "steps"	
		return mask1
	else:
		print "ne converge pas"
		return mask1, diff
	
def mask(noise):
	full_mask = np.zeros_like(noise)
	for z in range(noise.shape[0]):
		px_ref=(692,88)
#		print "Slice", z, ":",
		full_mask[z] = iter_mask(noise[z],px_ref=px_ref)
	return full_mask
	
def compact_x(map):
	l,m,n = map.shape
	comp = np.ones((l,m/2,n))
	for j in range(m/2):
		comp[:,j,:] = (map[:,2*j,:]+map[:,2*j+1,:])/2.
	return comp
	
def semi_log(map):
	return (map<=0)*map + (map>0)*np.log(1+map)

#def find_bin(r, bins):
#	found = False
#	current = 0
#	while (not(found)) & (current< bins.shape[0]-1):
#		current += 1
#		found = r < bins[current]
#	return current, found
#
#def correlation(data,bins,rad):
#	value = np.zeros(len(bins))
#	count = np.zeros(len(bins))
#	for i in range(data.shape[0]):
#		for j in range(data.shape[1]):
#			for k in range(data.shape[2]):
#				r = rad[i,j,k]
#				ind,valid = find_bin(r, bins)
#				if valid :
#					value[ind] += data[i,j,k]
#					count[ind] += 1
#	return value[1:],count[1:]

	
def dec_symetry(map):
	l,m,n = map.shape
	map_sym = np.zeros_like(map)
	for i in range(m):
		map_sym[:,i,:] = map[:,m-1-i,:]
	return map_sym
	
def process_cats(path, n, bins):
	sel_func = al.load(path+'sele_map_2df.npy')
	info = cp.deepcopy(sel_func.info)
	temp = al.make_vect(np.ones((1,1,1)),('z','x','y'))
	temp.set_axis_info('z',info['freq_centre'],info['freq_delta'])
	temp.set_axis_info('x',info['ra_centre'],info['ra_delta'])
	temp.set_axis_info('y',info['dec_centre'],info['dec_delta'])
	info = cp.deepcopy(temp.info)
	p_specs = {}
	for i in range(n):
		print '-------------------Step %s----------------'%i
		cat = al.load(path+'real_map_2df_ycli_50_%s.npy'%i)
		delta = 2*cat/sel_func - 1
		del cat
		delta[np.isnan(delta)] = 0
		delta = al.make_vect(delta)
		delta_xyz = physical_grid_largeangle(delta,order=1,refinement=1)[0]
		del delta
		delta_xyz = al.make_vect(delta_xyz)
		delta_xyz = rescale(delta_xyz,epsilon=1e-2)
		delta_xyz.info = cp.deepcopy(info)
		map_s = MapBuilder(delta_xyz[3:103])
		del delta_xyz
		p_specs[i] = power_specs(map_s,bins)
		del map_s
	return p_specs
	
try:
	print max([])
except (MemoryError,ValueError):
	print 'Error raised'
		
	