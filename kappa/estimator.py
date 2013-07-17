import numpy as np
import core.algebra as al
from kappa.my_utils import *
import kappa.filters as filters
from utils import batch_handler
import time
import copy as cp

class KappaEstimator:
	"""Class where every useful tools for tidal field reconstruction is stored

		attributes :
			-map : initial map
			-w : gaussian window
			-map_smooth : map convoluted with w
			-t_set : quadratic estimators of map
			-k_set : kernel (in Fourier space)
			-kappa : reconstructed density map			
			-hi_map : high modes of the original map
			-lo_map : low modes of the original map
			-kappa_hi : reconstructed density map from hi_map
			
			-sigma_recon : sigma of the gaussian window used in reconstruction
			-sigma_smooth : sigma of the gaussian window used in display
			-k_hi : lowest mode in hi_map
			-k_lo : highest mode in lo_map
	"""	
	
	def __init__(self, original_map, params, mask=None): # sigma_gaussian=1.25, sigma_smooth=8., full=True, verbose=True, fast=False, auto_rescale=True):
		self.map = np.copy(original_map)
		self.info = cp.deepcopy(original_map.info)
		self.params = cp.deepcopy(params)
		
		self.info['shape'] = original_map.shape
		self.info['sigma_recon'] = params['sigma_recon']
		self.info['sigma_smooth'] = params['sigma_smooth']
		self.z = original_map.get_axis('z') - self.info['z_centre']
		self.x = original_map.get_axis('x') - self.info['x_centre']
		self.y = original_map.get_axis('y') - self.info['y_centre']
		self.mask = mask
		
		if params['fast']:
			sigma = self.params['sigma_recon']
			if (sigma<=self.info['z_delta'])|(sigma<=self.info['x_delta'])|(sigma<=self.info['y_delta']):
				print "Warning : sigma for Gaussian window smaller than pixel size\nWith fast version, problem in the reconstruction might occur, consider turning flag 'fast' to False"
		
	def run(self):
		sigma1 = self.params['sigma_recon']
		self.w = self.gaussian_window(sigma1)
		self.info['w_sum'] = np.sum(self.w)
		if self.info['w_sum']>1:
			self.w = self.w/self.info['w_sum']
		if self.params['smoothing']:
			sigma2 = self.params['sigma_smooth']
			w_smooth = self.gaussian_window(sigma2)
		if self.params['verbose']: print "Gaussian window generated\n"
		
		if self.params['smoothing']:
			self.map_smooth = own_fftconvolve(self.map,w_smooth)
			if self.params['verbose']: print "Smoothing with gaussian window complete\n"
		
		if self.params['fast']: self.t_set = self.fast_quad_est(auto_rescale=self.params['auto_rescale'])
		else: self.t_set = self.quadratic_estimators(self.map, auto_rescale=self.params['auto_rescale'])
		if self.params['verbose']: print "Quadratic estimators complete\n"

		self.k_set = self.kernel()
		if self.params['verbose']: print "Kernel generated\n"

		self.kappa = self.density_field(self.t_set)
		if self.params['verbose']: print "Density field reconstructed\n"

		if not self.params['save_quad_est']: del self.t_set
		if not self.params['save_kernel']: del self.k_set

		self.noisy_modes()
		if self.params['clean']:
			self.apply_weight()
			if self.params['verbose']: print "Noisy modes filtered\n"
			
			if self.params['smoothing']:
				self.kappa_smooth = own_fftconvolve(self.clean, w_smooth)
				if self.params['verbose']: print "Kappa smoothing complete\n"

		print "Done"
	
	
	def get_k_axis(self, axis_name):
		size = {'z':self.info['shape'][0],
				'x':self.info['shape'][1],
				'y':self.info['shape'][2]}
		n = size[axis_name]
		Dk = 2*np.pi/self.info[axis_name + '_delta']
		return (np.arange(n) - n/2)*Dk/n
			
		
	@batch_handler.log_timing
	def gaussian_window(self, sigma, axes=None):
		l,m,n = self.info['shape']
		dz = self.info['z_delta']
		dx = self.info['x_delta']
		dy = self.info['y_delta']
		if axes==None:
			rad_3d = (self.z*self.z)[:,None,None] * np.ones((1,m,n))
			rad_3d += (self.x*self.x)[None,:,None] * np.ones((l,1,n))
			rad_3d += (self.y*self.y)[None,None,:] * np.ones((l,m,1))
			return dz*dx*dy*(sigma*np.sqrt(2*np.pi))**(-3) * np.exp(-rad_3d/(2.*sigma*sigma))
		else:
			z,x,y = axes
			rad_3d = (z*z)[:,None,None] * np.ones((1,m,n))
			rad_3d += (x*x)[None,:,None] * np.ones((l,1,n))
			rad_3d += (y*y)[None,None,:] * np.ones((l,m,1))
			return np.exp(-rad_3d/(2.*sigma*sigma))


	@batch_handler.log_timing	
	def partial_derivatives(self):
		l,m,n = self.info['shape']
		sigma = self.params['sigma_recon']
		x_0 = self.x[m/2]
		y_0 = self.y[n/2]
		w_x = np.zeros((l,m,n))
		w_y = np.zeros((l,m,n))
		for i in range(m):
			x_i = self.x[i] - x_0
			w_x[:,i,:] = -x_i/(sigma*sigma)*self.w[:,i,:]
		for j in range(n):
			y_j = self.y[j] - y_0
			w_y[:,:,j] = -y_j/(sigma*sigma)*self.w[:,:,j]
		return w_x,w_y

	
	@batch_handler.log_timing
	def quadratic_estimators(self, map, auto_rescale=True):
		(w_x,w_y) = self.partial_derivatives()
		delta_x = own_fftconvolve(map,w_x)
		delta_y = own_fftconvolve(map,w_y)
		map_log = own_fftconvolve(map,self.w)
		if map_log.min()<=(-1.+1.e-2):
			if auto_rescale:
				print "Min value of smoothed map is below -1 ! => %f"%map_log.min()
				print "Rescaling with epsilon = 1e-2"
				map_log = rescale(map_log,epsilon=1e-2)
			else:
				assert map_log.min()>-1., "Min value of map_log is below -1 ! => %f"%map_log.min()
		delta_x = delta_x/(1.+map_log)
		delta_y = delta_y/(1.+map_log)
		if self.params['verbose']: print "Logarithmic smoothing complete"
		t_xx = delta_x*delta_x
		t_yy = delta_y*delta_y
		t_xy = delta_x*delta_y
		if self.params['save_map_log']: self.map_log = map_log
		return t_xx,t_yy,t_xy
	
	@batch_handler.log_timing
	def fast_quad_est(self, auto_rescale=True):
		l,m,n = self.info['shape']
		sigma_k = self.params['sigma_recon']**(-1)
		kz = self.get_k_axis('z')
		kx = self.get_k_axis('x')
		ky = self.get_k_axis('y')
		w = self.gaussian_window(sigma_k, axes=[kz,kx,ky])
		
		d_x = 1.j*kx[None,:,None] * w
		d_y = 1.j*ky[None,None,:] * w
		d_x = np.fft.ifftshift(d_x)
		d_y = np.fft.ifftshift(d_y)
		f_map = np.fft.ifftshift(self.map)
		t1 = time.time()
		f_map = np.fft.fftn(f_map)
		t2 = time.time()
		print '[TIMING] fft:',t2-t1
		d_x = d_x * f_map
		d_y = d_y * f_map
		t1 = time.time()
		print '[TIMING] multiplication:',t1-t2
		d_x = np.fft.ifftn(d_x).real
		t2 = time.time()
		print '[TIMING] ifft:',t2-t1
		d_y = np.fft.ifftn(d_y).real
		t1 = time.time()
		print '[TIMING] ifft:',t1-t2
		d_x = np.fft.fftshift(d_x)
		d_y = np.fft.fftshift(d_y)
		
		if self.info['w_sum']>1:
			d_x = d_x/self.info['w_sum']
			d_y = d_y/self.info['w_sum']

		map_log = f_map*np.fft.ifftshift(w)
		t2 = time.time()
		map_log = np.fft.ifftn(map_log)
		map_log = np.fft.fftshift(map_log.real)
		t3 = time.time()
		print '[TIMING] ifft:',t3-t2
		if map_log.min()<=(-1.+1.e-2):
			if auto_rescale:
				print "Min value of smoothed map is below -1 ! => %f"%map_log.min()
				print "Rescaling with epsilon = 1e-2"
				map_log = cut(map_log,epsilon=1e-2)
			else:
				assert map_log.min()>-1., "Min value of map_log is below -1 ! => %f"%map_log.min()
		d_x = d_x/(1. + map_log)
		d_y = d_y/(1. + map_log)
		if self.params['save_map_log']: self.map_log = map_log
		return d_x*d_x,d_y*d_y,d_x*d_y
	
	@batch_handler.log_timing	
	def kernel(self):
		kz = self.get_k_axis('z')
		kx = self.get_k_axis('x')
		ky = self.get_k_axis('y')
		k_xx = (kx*kx)[:,None] * np.ones((1,len(ky)))
		k_yy = np.ones((len(kx),1)) * (ky*ky)[None,:]
		k_xy = kx[:,None] * ky[None,:]
		rad_2d = k_xx + k_yy #calculate kx^2 + ky^2
		k_xx = k_xx/rad_2d - 0.5
		k_yy = k_yy/rad_2d - 0.5
		k_xy = k_xy/rad_2d
		k_z_3d = (kz*kz)[:,None,None] * (1/rad_2d)[None,:,:]
		k_xx = (1 + k_z_3d)*k_xx[None,:,:]
		k_yy = (1 + k_z_3d)*k_yy[None,:,:]
		k_xy = (1 + k_z_3d)*k_xy[None,:,:]
		k_xx[:,len(kx)/2,len(ky)/2] = 0. #Temporary solution, ok only if x_centre,y_centre = 0
		k_yy[:,len(kx)/2,len(ky)/2] = 0.
		k_xy[:,len(kx)/2,len(ky)/2] = 0.
			#Shift so that the function ifftn work
		k_xx = np.fft.ifftshift(k_xx)
		k_yy = np.fft.ifftshift(k_yy)
		k_xy = np.fft.ifftshift(k_xy)
		return k_xx,k_yy,k_xy

	@batch_handler.log_timing
	def density_field(self, t_set):
		(t_xx,t_yy,t_xy) = t_set
			#Shift so that x=0 is the first element of the matrix
		t_xx = np.fft.ifftshift(t_xx)
		t_yy = np.fft.ifftshift(t_yy)
		t_xy = np.fft.ifftshift(t_xy)
			#Apply the fft
		t_xx = np.fft.fftn(t_xx)
		t_yy = np.fft.fftn(t_yy)
		t_xy = np.fft.fftn(t_xy)
			#Get kappa in k space
		kappa_xx = t_xx * self.k_set[0]
		kappa_yy = t_yy * self.k_set[1]
		kappa_xy = t_xy * self.k_set[2]
			#Merge contributions to kappa
		kappa = kappa_xx + kappa_yy + 2*kappa_xy
			#Get kappa in real space
		kappa = np.fft.ifftn(kappa)
			#Keep only the real part (imaginary part should be negligible)
		kappa = kappa.real
			#Shift so that x=0 is the centre of the matrix
		kappa = np.fft.fftshift(kappa)
		return kappa
	

	def modify_smoothing(self, new_sigma):
		self.params['sigma_smooth'] = new_sigma
		w_smooth = self.gaussian_window(new_sigma)
		self.map_smooth = own_fftconvolve(self.map,w_smooth)
		self.kappa_smooth = own_fftconvolve(self.clean,w_smooth) 

	
	@batch_handler.log_timing
	def modify_filters(self, new_k_hi, new_k_lo=None, reconstruct=False, kz_filt=False):
		if new_k_lo==None : new_k_lo = new_k_hi
		self.k_hi = new_k_hi
		self.k_lo = new_k_lo
		
		kz = self.get_k_axis('z')
		kx = self.get_k_axis('x')
		ky = self.get_k_axis('y')
		self.hi_map,self.lo_map = filters.spherical_separation_filter(self.map,[kz,kx,ky],new_k_hi,new_k_lo)
		if self.mask!=None:
			self.hi_map *= self.mask
			self.lo_map *= self.mask		

		if reconstruct:
			print "High map values between :",self.hi_map.min(),",",self.hi_map.max()
			if self.hi_map.min()<=-1:
				print "rescaling..."
				self.hi_map = rescale(self.hi_map, epsilon=1e-2)
			self.kappa_hi = self.density_field(self.quadratic_estimators(self.hi_map,log_smooth=True))
			if kz_filt:
				self.kappa_hi = apply_weight(kappa=self.kappa_hi,output=True)

		print "New filters applied"


	def noisy_modes(self):
		kz = self.get_k_axis('z')
		kx = self.get_k_axis('x')
		ky = self.get_k_axis('y')
		l = len(kz)
		m = len(kx)
		n = len(ky)
		rad_2d = (kx*kx)[:,None]*np.ones((1,n))
		rad_2d += (ky*ky)[None,:]*np.ones((m,1))
		rad_2d = rad_2d[None,:,:]*np.ones((l,1,1))
		self.weight = rad_2d/(rad_2d + (kz*kz)[:,None,None]*np.ones((1,m,n)))
		self.weight = np.fft.ifftshift(self.weight)		
		self.weight[0,0,0] = 0

		
	def apply_weight(self, kappa=None, output=False):
		if kappa==None: kappa=self.kappa
		try:
			clean = np.fft.ifftshift(kappa)
			clean = np.fft.fftn(clean)
			clean *= self.weight
			clean = np.fft.ifftn(clean)
			clean = np.fft.fftshift(clean.real)
			if output:
				return clean
			else:
				self.clean = clean
		except NameError:
			print "Error in apply_weight : weights are not defined!"
		
