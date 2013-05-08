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
def kernel(k,i,j)
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


            
#animation to show every slices z
fig, ax = plt.subplots()

im = ax.imshow(kappa[2,:,:].T)
fig.show()
for i in range (128):
	im.set_data(kappa[i,:,:].T)
	fig.canvas.draw()
	time.sleep(.2)
	


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

	


