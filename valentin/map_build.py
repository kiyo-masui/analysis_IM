import numpy as np
from valentin.my_utils import *
from valentin.filters import *
from utils import batch_handler
import time
import copy as cp

class MapBuilder:
	"""Class where every useful tools for density reconstruction is stored

		attributes :
			-map : initial map
			-w : gaussian window
			-map_smooth : map convolved with w
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
	
	def __init__(self, original_map, mask=None, sigma_gaussian=1.25, sigma_smooth=8., log_smooth=True, full=True, verbose=True, fast=False, auto_rescale=True):
#		if original_map.min()<=-1.:
#			print "Min value of input map is below -1 ! => %f"%original_map.min()
#			if auto_rescale:
#				print "Rescaling with epsilon = 1e-2"
#				original_map = rescale(original_map,epsilon=1e-2)
#			else:
#				return
		self.map = np.copy(original_map)
		self.info = cp.deepcopy(original_map.info)
		self.info['shape'] = original_map.shape
		self.z = original_map.get_axis('z') - self.info['z_centre']
		self.x = original_map.get_axis('x') - self.info['x_centre']
		self.y = original_map.get_axis('y') - self.info['y_centre']
		self.info['sigma_recon'] = sigma_gaussian #* self.info['x_delta']
		self.info['sigma_smooth'] = sigma_smooth
		
#		if mask!=None: 
		self.mask = mask
#		else : self.mask = np.ones_like(original_map)	
		
		sigma1 = self.info['sigma_recon']
		self.w = self.gaussian_window(sigma1)
#		if full==False: return
		test=np.sum(self.w)
		if test>1:
			self.w = self.w/test
		if full:
			sigma2 = self.info['sigma_smooth']
			w_smooth = self.gaussian_window(sigma2)
		if verbose: print "Gaussian windows generated"
		
		if full:
			self.map_smooth = own_fftconvolve(self.map,w_smooth)
			if verbose: print "Smoothing with gaussian window complete"
		
		if not fast: self.t_set = self.quadratic_estimators(self.map,sigma1,log_smooth)
		else: self.t_set = self.fast_quad_est()
		if verbose: print "Quadratic estimators complete"

		self.k_set = self.fast_kernel()
		if verbose: print "Kernel generated"

		self.kappa = self.density_field(self.t_set)
		if verbose: print "Density field reconstructed"

		if full:
			self.noisy_modes()
			self.final = np.fft.ifftshift(self.kappa)
			self.final = np.fft.fftn(self.final)
			self.final *= self.noise
			self.final = np.fft.ifftn(self.final)
			self.final = np.fft.fftshift(self.final.real)
#			self.final = kz_filter(self.kappa,[self.get_k_axis('z'),self.get_k_axis('x'),self.get_k_axis('y')])
			if verbose: print "Noisy modes filtered"
			
			self.kappa_smooth = own_fftconvolve(self.final, w_smooth)
			if verbose: print "Kappa smoothing complete"

#			self.k_hi = 0.2
#			self.k_lo = 0.2
#			self.hi_map,self.lo_map = spherical_separation_filter(self.map,[1,1,1],self.k_hi,self.k_lo)
#			self.kappa_hi = self.density_field(self.quadratic_estimators(self.hi_map,sigma))
#			self.kappa_lo = self.density_field(self.quadratic_estimators(self.lo_map,sigma))
#			print "High modes filtered map processed"
##		self.correlations()
##		print "Correlations computed"

		print "Done"
	
	def get_k_axis(self, axis_name):
		size = {'z':self.info['shape'][0],
				'x':self.info['shape'][1],
				'y':self.info['shape'][2]}
		n = size[axis_name]
		Dk = 2*np.pi/self.info[axis_name + '_delta']
		return (np.arange(n) - n/2)*Dk/n
			
#	@batch_handler.log_timing
#	def gaussian_window(self, sigma):
#		mesh = np.indices(self.info['shape'])
#		sph = np.zeros(mesh.shape)
#		sph[0] = self.z[mesh[0]] - self.info['z_centre']
#		sph[1] = self.x[mesh[1]] - self.info['x_centre']
#		sph[2] = self.y[mesh[2]] - self.info['y_centre']
#		rad = sph[0]**2 + sph[1]**2 + sph[2]**2
#		dz = self.info['z_delta']
#		dx = self.info['x_delta']
#		dy = self.info['y_delta'] 
#		return dz*dx*dy*(sigma*np.sqrt(2*np.pi))**(-3) * np.exp(-rad/(2.*sigma*sigma))
		
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
	def partial_derivatives(self, sigma=1.):
		l,m,n = self.info['shape']
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
	def quadratic_estimators(self, map, sigma=1., log_smooth=False, auto_rescale=True):
		t = time.time()
		(w_x,w_y) = self.partial_derivatives(sigma)
		t1 = time.time() - t
#		print 'Partial derivatives : %.3f' %t1
		delta_x = own_fftconvolve(map,w_x)
		delta_y = own_fftconvolve(map,w_y)
		map_log = own_fftconvolve(map,self.w)
		if log_smooth:
			if map_log.min()<=(-1.+1.e-2):
				if auto_rescale:
					print "Min value of smoothed map is below -1 ! => %f"%map_log.min()
					print "Rescaling with epsilon = 1e-2"
					map_log = rescale(map_log,epsilon=1e-2)
				else:
					assert map_log.min()>-1., "Min value of map_log is below -1 ! => %f"%map_log.min()
			print map_log.min()
#			map_log = rescale(map_log)
#			semi_log = (map_log<=0) + (map_log>0)*(1.+map_log)
#			delta_x = delta_x/semi_log
#			delta_y = delta_y/semi_log
			delta_x = delta_x/(1.+map_log)
			delta_y = delta_y/(1.+map_log)
			print "Logarithmic smoothing complete"
		t_xx = delta_x*delta_x
		t_yy = delta_y*delta_y
		t_xy = delta_x*delta_y
		t = time.time() - t
#		print t
		self.map_log = map_log
		return t_xx,t_yy,t_xy
	
	@batch_handler.log_timing
	def fast_quad_est(self):
		l,m,n = self.info['shape']
		sigma_k = self.info['sigma_recon']**(-1)
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

		t1 = time.time()
		map_log = f_map*np.fft.ifftshift(w)
		t2 = time.time()
		map_log = np.fft.ifftn(map_log)
		map_log = np.fft.fftshift(map_log.real)
		t3 = time.time()
		print '[TIMING] ifft:',t3-t2
#		if map_log.min()<=-1. : map_log = rescale(map_log)
		assert map_log.min()>-1., "Min value of map_log is below -1 ! => %f"%map_log.min()
		d_x = d_x/(1. + map_log)
		d_y = d_y/(1. + map_log)
		t2 = time.time()
		print '[TIMING] log_smooth:',t2-t1
		return d_x*d_x,d_y*d_y,d_x*d_y
	
#	def map_kernel(self):
#		t = time.time()
#		l,m,n = self.info['shape']
#		ker_xx = np.zeros((l,m,n))
#		ker_yy = np.zeros((l,m,n))
#		ker_xy = np.zeros((l,m,n))
#		
#		k_axes = []
#		k_axes.append(self.get_k_axis('z'))		
#		k_axes.append(self.get_k_axis('x'))
#		k_axes.append(self.get_k_axis('y'))
#		
#		for k in range(l):
#			for i in range(m):
#				for j in range(n):
#					k_z = k_axes[0][k]
#					k_x = k_axes[1][i]
#					k_y = k_axes[2][j]
#					cst2 = k_x*k_x + k_y*k_y + 0.0
#					if (cst2 != 0):
#						cst = (cst2 + k_z*k_z)/cst2
#						ker_xx[k,i,j] = cst*(k_x*k_x/cst2 - 0.5)
#						ker_yy[k,i,j] = cst*(k_y*k_y/cst2 - 0.5)
#						ker_xy[k,i,j] = cst*(k_x*k_y/cst2)
#						
#		#Shift so that the fonction ifftn work
#		ker_xx = np.fft.ifftshift(ker_xx)
#		ker_yy = np.fft.ifftshift(ker_yy)
#		ker_xy = np.fft.ifftshift(ker_xy)
#		
#		t = time.time() - t
#		print t
#		return ker_xx,ker_yy,ker_xy
	
	@batch_handler.log_timing	
	def fast_kernel(self):
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
		t = time.time()
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
#			#Shift so that x=0 is the centre of the matrix
		kappa = np.fft.fftshift(kappa)
		t = time.time() - t
#		print t
		return kappa
	
	@batch_handler.log_timing
	def modify_smoothing(self, new_sigma):
		self.info['sigma_smooth'] = new_sigma
		w_smooth = self.gaussian_window(new_sigma)
		self.map_smooth = own_fftconvolve(self.map,w_smooth)
		self.kappa_smooth = own_fftconvolve(self.final,w_smooth) 
	
	@batch_handler.log_timing
	def modify_filters(self, new_k_hi, new_k_lo=None, reconstruct=False, kz_filt=False):
		if new_k_lo==None : new_k_lo = new_k_hi
		sigma = self.info['sigma_recon']
		self.k_hi = new_k_hi
		self.k_lo = new_k_lo
		
		kz = self.get_k_axis('z')
		kx = self.get_k_axis('x')
		ky = self.get_k_axis('y')
#		dk = k[1]-k[0]
		self.hi_map,self.lo_map = spherical_separation_filter(self.map,[kz,kx,ky],new_k_hi,new_k_lo)
		if self.mask!=None:
			self.hi_map *= self.mask
			self.lo_map *= self.mask
		
		if reconstruct:
			print "High map values between :",self.hi_map.min(),",",self.hi_map.max()
			print "rescaling..."
			self.hi_map = rescale(self.hi_map, epsilon=1e-2)
			self.kappa_hi = self.density_field(self.quadratic_estimators(self.hi_map,sigma,log_smooth=True))
			if kz_filt:
#				self.kappa_hi = kz_filter(self.kappa_hi,[self.get_k_axis('z'),self.get_k_axis('x'),self.get_k_axis('y')])
				self.kappa_hi = np.fft.ifftshift(self.kappa_hi)
				self.kappa_hi = np.fft.fftn(self.kappa_hi)
				self.kappa_hi *= self.noise
				self.kappa_hi = np.fft.ifftn(self.kappa_hi)
				self.kappa_hi = np.fft.fftshift(self.kappa_hi.real)
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
		self.noise = rad_2d/(rad_2d + (kz*kz)[:,None,None]*np.ones((1,m,n)))
		self.noise = np.fft.ifftshift(self.noise)		
		self.noise[0,0,0] = 0

#	def correlations(self):
#		f_hi = np.fft.ifftshift(self.hi_map)
#		f_hi = np.fft.fftn(f_hi)
#		f_lo = np.fft.ifftshift(self.lo_map)
#		f_lo = np.fft.fftn(f_lo)
#		corr_filt = f_hi * f_lo.conjugate()
#		corr_filt = np.fft.ifftn(corr_filt)
#		corr_filt = corr_filt.real
#		corr_filt = np.fft.fftshift(corr_filt)
#
#		f_kap_hi = np.fft.ifftshift(self.kappa_hi)
#		f_kap_hi = np.fft.fftn(f_kap_hi)
#		f_lo = np.fft.ifftshift(self.lo_map)
#		f_lo = np.fft.fftn(f_lo)
#		corr_kap = f_kap_hi * f_lo.conjugate()
#		corr_kap = np.fft.ifftn(corr_kap)
#		corr_kap = corr_kap.real
#		corr_kap = np.fft.fftshift(corr_kap)
#		
##		self.corr_f = corr_filt
##		self.corr_k = corr_kap
#		return corr_filt,corr_kap
