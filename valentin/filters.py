import numpy as np

#First version, cubic maps
def cubic_separation_filter(map, hi_k, lo_k):
	l,m,n = map.shape
	f_map = np.fft.fftn(map)
	hi_map = np.ones(map.shape)
	lo_map = np.zeros(map.shape)
	
	hi_map[(l/2-hi_k):(l/2 + hi_k + 1),(l/2-hi_k):(l/2 + hi_k + 1),(l/2-hi_k):(l/2 + hi_k + 1)] = np.zeros((2*hi_k+1,2*hi_k+1,2*hi_k+1))
	lo_map[(l/2-lo_k):(l/2 + lo_k + 1),(l/2-lo_k):(l/2 + lo_k + 1),(l/2-lo_k):(l/2 + lo_k + 1)] = np.ones((2*lo_k+1,2*lo_k+1,2*lo_k+1))
	
	hi_map = np.fft.ifftshift(hi_map)
	lo_map = np.fft.ifftshift(lo_map)
	
	hi_map = f_map*hi_map
	lo_map = f_map*lo_map
	
	hi_map = np.fft.ifftn(hi_map)
	lo_map = np.fft.ifftn(lo_map)
	
	hi_map = hi_map.real
	lo_map = lo_map.real
	return hi_map,lo_map


def spherical_separation_filter(map, k_axes, hi_k, lo_k=None, high=True, low=True, isReal=True):
	if low & (lo_k==None): lo_k=hi_k	
	l,m,n = map.shape
	if isReal: 
		f_map = np.fft.ifftshift(map)
		f_map = np.fft.fftn(f_map)
	else : 
		f_map = map
	
#	#Change for alternate version (cf gaussian window)
#	mesh = 1.0*np.mgrid[0:l,0:m,0:n]
#	mesh[0] -= l/2
#	mesh[1] -= m/2
#	mesh[2] -= n/2
#	mesh[0] *= d_k[0]
#	mesh[1] *= d_k[1]
#	mesh[2] *= d_k[2]
#	rad_3d = mesh[0]*mesh[0] + mesh[1]*mesh[1] + mesh[2]*mesh[2]
	assert len(k_axes)==3, "Probleme avec les axes dans separation_filter"
	kz,kx,ky = k_axes
#	kz = d_k[0]*(np.arange(l) - l/2)
#	kx = d_k[1]*(np.arange(m) - n/2)
#	ky = d_k[2]*(np.arange(n) - m/2)
	rad_3d = (kz*kz)[:,None,None] * np.ones((1,m,n))
	rad_3d += (kx*kx)[None,:,None] * np.ones((l,1,n))
	rad_3d += (ky*ky)[None,None,:] * np.ones((l,m,1))
	
	res=[]
	
	if high:
		hi_map = rad_3d>(hi_k*hi_k)
		hi_map = np.fft.ifftshift(hi_map)
		hi_map = f_map*hi_map
		hi_map = np.fft.ifftn(hi_map)
		hi_map = np.fft.fftshift(hi_map.real)
		res.append(hi_map)
	
	if low:
		lo_map = rad_3d<=(lo_k*lo_k)
		lo_map = np.fft.ifftshift(lo_map)
		lo_map = f_map*lo_map
		lo_map = np.fft.ifftn(lo_map)
		lo_map = np.fft.fftshift(lo_map.real)
		res.append(lo_map)
	
	if len(res)==1: return res[0]
	else: return res

def kappa_filter(kappa, hi_k, axes):
	assert len(axes)==3, "Probleme avec les axes dans kappa_filter"
	return spherical_separation_filter(kappa, axes, hi_k, hi_k, high=False)


def autocorr(map_s, n=20):
	r = np.zeros(n)
	r_ref = np.zeros(n)
	w = map_s.gaussian_window(8.)
	for j in range(n):
		k = j*0.48/30.
		print "step",j,": k=",k
		map_s.modify_filters(k,k, reconstruct=True, kz_filt=True)
		map_s.kappa_hi = own_fftconvolve(map_s.kappa_hi,w)
		if map_s.mask!=None:
			map_s.kappa_hi = map_s.mask*map_s.kappa_hi
#			map_s.lo_map = map_s.mask*map_s.lo_map
		m1m1 = np.sum(map_s.kappa_hi * map_s.kappa_hi)
		m2m2 = np.sum(map_s.lo_map * map_s.lo_map)
		m1m2 = np.sum(map_s.kappa_hi * map_s.lo_map)
		r[j] = m1m2/np.sqrt(m1m1*m2m2)
		print r[j]
		m1m1b = np.sum(map_s.hi_map * map_s.hi_map)
		m1m2b = np.sum(map_s.hi_map * map_s.lo_map)
		r_ref[j] = m1m2b/np.sqrt(m1m1b*m2m2)
	return r, r_ref

def kappa_autocorr(map_s, dk, mask=None, kappa=None):
	if kappa==None: kappa = map_s.kappa
	kz = map_s.get_k_axis('z')
	kx = map_s.get_k_axis('x')
	ky = map_s.get_k_axis('y')
	r = np.zeros(20)
	for j in range(20):
		k = j*dk
		print "step",j,": k=",k
		map_s.modify_filters(k,k)
#		kap_tmp = map_s.kappa_hi
		kap_tmp = kappa_filter(kappa,k,[kz,kx,ky])
		if mask!=None:
			kap_tmp = mask*kap_tmp
#			map_s.lo_map = mask*map_s.lo_map
		m1m1 = np.sum(kap_tmp * kap_tmp)
		m2m2 = np.sum(map_s.lo_map * map_s.lo_map)
		m1m2 = np.sum(kap_tmp * map_s.lo_map)
		r[j] = m1m2/np.sqrt(m1m1*m2m2)
		print r[j]
	return r


def kz_filter(kappa, k_axes):
	l,m,n = kappa.shape
	kz = np.fft.ifftshift(k_axes[0])
	kx = np.fft.ifftshift(k_axes[1])
	ky = np.fft.ifftshift(k_axes[2])
	rad_2d = (kx*kx)[None,:,None] * np.ones((l,1,n))
	rad_2d += (ky*ky)[None,None,:] * np.ones((l,m,1))
	filt = rad_2d - (kz*kz)[:,None,None]*np.ones((1,m,n))
	filt = filt > 0
	f_kappa = np.fft.fftn(kappa)
	f_kappa *= filt
	return np.fft.ifftn(f_kappa).real
	

