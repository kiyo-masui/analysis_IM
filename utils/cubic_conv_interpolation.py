import numpy as np
import matplotlib.pyplot as plt
	

def cross(set_list):
	'''Given a list of sets, return the cross product.'''
	
	# By associativity of cross product, cross the first two sets together
	# then cross that with the rest. The big conditional in the list
	# comprehension is just to make sure that there are no nested lists
	# in the final answer.
	if len(set_list) == 1:
		ans = []
		for elem in set_list[0]:
			# In the 1D case, these elements are not lists.
			if type(elem) == list:
				ans.append(np.array(elem))
			else:
				ans.append(np.array([elem]))
		return ans
	else:
		A = set_list[0]
		B = set_list[1]
		cross_2 = [a+b if ((type(a) == list) and (type(b) == list)) else \
		(a+[b] if type(a) == list else ([a]+b if type(b) == list else \
		[a]+[b])) for a in A for b in B]
		remaining = set_list[2:]
		remaining.insert(0,cross_2)
		return cross(remaining)
		
		
def point_list(node):
	'''Return the list of lists containing the numbers, in order, of the
	indeces that must be looked at in the interpolation. node is a
	coordinate.'''
	
	ans = []
	for coord in node:
		ans.append([coord-1,coord,coord+1,coord+2])
	return ans
	
	
def left_interpolation_node(pnt, x0, step_sizes):
	'''Suppose pnt (list of coordinates in each dimension) is where the 
	sampled data is to be interpolated. Then pnt must be between two
	consecutive interpolation nodes (in each dimension), x_j and x_(j+1).
	Return x_j. The first point x0 and the step size in each dimension
	is needed for x_j to be found.
	Note: all parameters must have same dimensions.
	Note: returns the node in indeces, not actual units.'''
	
	pnt = np.array(pnt, dtype=float)
	left_node = np.floor((pnt - x0)/step_sizes)
	return np.array(left_node,dtype=int)
	
	
def get_value(arr,node):
	'''Get the value of n-dimensional array arr at n-coordinate index node.
	Note node is a 1D array in pixels, not actual units.'''
	
	# Copy node since we needed it later and shouldn't mess with it.
	node_copy = node.copy()
	# The max index possible for each dimension of arr.
	max_indeces = np.array(arr.shape) - 1
	# Initialize the bad_index. If it's stll -99 after looking through
	# all the coordinates of pnt, then pnt is in arr.
	bad_index = -99
	# True if one coordinate of pnt is under the min possible (0) and
	# false if one coordinate of pnt is over the max possible.
	under = False
	# Find any index out of bounds (lower).
	for i in range(len(node)):
		if node[i] < 0:
			bad_index = i
			under = True
		elif node[i] > max_indeces[i]:
			bad_index = i
			under = False
			
	if bad_index == -99:
		# No bad indeces so just return the value of array at that point.
		return arr[tuple(node)]
	elif under:
		# There is a bad index and it is -1. Use method similar to set of
		# equations on page 1158 of article to find value.
		node_copy[bad_index] += 1
		fst_node = node_copy.copy()
		node_copy[bad_index] += 1
		snd_node = node_copy.copy()
		node_copy[bad_index] += 1
		thd_node = node_copy.copy()
		return 3.*get_value(arr,fst_node) - 3.*get_value(arr,snd_node) + \
		get_value(arr,thd_node)
	else:
		# There is a bad index and it is one more than the max possible.
		# Note there cannot be any other cases.
		node_copy[bad_index] -= 1
		fst_node = node_copy.copy()
		node_copy[bad_index] -= 1
		snd_node = node_copy.copy()
		node_copy[bad_index] -= 1
		thd_node = node_copy.copy()
		return 3.*get_value(arr,fst_node) - 3.*get_value(arr,snd_node) + \
		get_value(arr,thd_node)

		
def u(s):
	'''Return the cubic convolution interpolation kernel given the pixel
	distance. See equation 15 from Robert G. Keys 1981 article.'''
	
	# Note: the s==0 and else cases are put in for correctness but
	# these should never be reached since this function does not need
	# to be called in those cases.
	s = abs(s)
	if 0 < s < 1:
		return 3./2.*s**3 - 5./2.*s**2 + 1.
	elif 1 < s < 2:
		return -1./2.*s**3 + 5./2.*s**2 - 4.*s + 2.
	elif s == 0:
		return 1.
	else:
		return 0.
		

def u_product(pnt, node, x0, step_sizes):
	'''Return the product of the "u" in each dimension. Want to get the 
	equivalent of the u part of the equation on pg 1157. pnt and node
	are arrays containing the point we wish to interpolate and the points
	with known values in the array, respectively. pnt is in actual units
	while node is in indeces. x0 is the initial point in n-dimensions.
	step_sizes is the distance between pixels in each dimension.'''
	
	step_sizes = np.array(step_sizes, dtype='float')
	node_pixel = x0 + node*step_sizes
	return np.product(map(u,(pnt-node_pixel)/step_sizes))
	

def reassign_weights(points, weights, max_inds):
	'''Return new points and weights such that the value and weight
	of any out of bound point is spread over in-bound points.
	max_ind is an array containing the max index allowed for each
	axis.'''
	# Keep track of bad points so that they can be assigned fake
	# inbound point indices afterwards.
	bad_points = []
	for i in range(len(points)):
		spread = find_weight_spread(points[i],points,max_inds)
		# If the point is bad, redistribute the weights.
		if not (np.sum(spread**2) == 1):
			bad_points.append(i)
			# Note that this can be done concurrently with finding the
			# bad points since the weights are only going to change for
			# in bound points and these points will never have their
			# weights spread over anything else.
			weights += spread*weights[i]
	# Assign fake indeces to bad points. Since the 0th index always exists
	# for an array, assign it that with a weight of zero.
	zero_index = np.zeros(len(max_inds),dtype='int')
	for j in bad_points:
		points[j] = zero_index
		weights[j] = 0.0
	return points,weights
	

def find_point_in_list(point, point_list):
	'''Return the index of a point in a point_list (or array).
	Note that point must occur once and only once in point_list.'''
	index = 0
	for i in range(len(point_list)):
		if tuple(point) == tuple(point_list[i]):
			index = i
	return index
	
	
def find_weight_spread(point, points, max_inds):
	'''Given a point "point" and the list of points needed
	for cubic interpolation, "points," return how much weight needs to
	be multiplied in each point in "points."
	max_ind is an array containing the max index allowed for each
	axis.'''
	point = np.array(point)
	min_inds = np.zeros(len(point))
	# Point is in bounds.
	if (False not in (point <= max_inds).tolist()) and \
	(False not in (point >= min_inds).tolist()):
		# Point should not spread its value around. So return
		# its location with a value of 1 and all other points
		# a value of 0.
		add_weights = np.zeros(len(points))
		one_index = find_point_in_list(point,points)
		add_weights[one_index] = 1.0
		return add_weights
	# Point is under bounds on at least one axis.
	elif (False in (point >= min_inds).tolist()):
		# Find the bad axis.
		for i in range(len(point)):
			if point[i] < 0:
				under_axis = i
		# Split its weight based on eqn 19. Note recursion
		# is to take care of points out of bounds on muliple axes.
		point_1 = np.copy(point)
		point_1[under_axis] += 1
		point_2 = np.copy(point)
		point_2[under_axis] += 2
		point_3 = np.copy(point)
		point_3[under_axis] += 3
		weights_1 = find_weight_spread(point_1, points, max_inds)
		weights_2 = find_weight_spread(point_2, points, max_inds)
		weights_3 = find_weight_spread(point_3, points, max_inds)
		return 3.0*weights_1 - 3.0*weights_2 + 1.0*weights_3
	# Point is over bounds on at least one axis.
	elif (False in (point <= max_inds).tolist()):
		# Find the bad axis.
		for i in range(len(point)):
			if point[i] > max_inds[i]:
				over_axis = i
		# Split its weight based on eqn 25. Note recursion
		# is to take care of points out of bounds on muliple axes.
		point_N = np.copy(point)
		point_N[over_axis] -= 1
		point_N1 = np.copy(point)
		point_N1[over_axis] -= 2
		point_N2 = np.copy(point)
		point_N2[over_axis] -= 3
		weights_N = find_weight_spread(point_N, points, max_inds)
		weights_N1 = find_weight_spread(point_N1, points, max_inds)
		weights_N2 = find_weight_spread(point_N2, points, max_inds)
		return 3.0*weights_N - 3.0*weights_N1 + 1.0*weights_N2
	else:
		print "You shouldn't be here..."
	

def	interpolate_weights(axes, pnt, x0, step_sizes, max_inds):
	'''Return the points and weights of all non-zero weighted points
	needed in the cubic interpolation. 'axes' is the list of axes to be
	interpolated over. Only the dimensions in axes will be looked at to
	get the appropriate nodes and all other dimensions will be ignored.
	See interpolate for description of other parameters.
	
	In the new interpolate, points that are out of the boundaries are
	split amongst points that are inboundaries so map maker doesn't
	crash.
	max_ind is an array containing the max index allowed for each
	axis. Now needed to be able to tell what an out of bound point is.'''
	
	# Get only the needed axes out of the info we have to be able to
	# call the code already written.
	new_pnt = []
	new_x0 = []
	new_step_sizes = []
	
	for i in range(len(axes)):
		new_pnt.append(pnt[i])
		new_x0.append(x0[i])
		new_step_sizes.append(step_sizes[i])
	
	new_pnt = np.array(new_pnt)
	new_x0 = np.array(new_x0)
	new_step_sizes = np.array(new_step_sizes)
	
	# Get the closest node to the "left".
	left_node = left_interpolation_node(new_pnt, new_x0, new_step_sizes)
	# Get the list of nodes that we need.
	needed_nodes = cross(point_list(left_node))
	
	# Get the corresponding list of weights.
	weights = []
	for node in needed_nodes:
		weights.append(u_product(new_pnt,node,new_x0,new_step_sizes))
	weights = np.array(weights,dtype=float)
	
	
	points = np.array(needed_nodes)
	points,weights = reassign_weights(points,weights,max_inds)
	return points,weights
	
	
def interpolate(axes, arr, pnt, x0, step_sizes):
	'''Return cubic convolution interpolation at m-coordinate point pnt.
	Note pnt is in the actual units of the axes. arr is the n-dimentional
	array that contains the known discrete values to be interpolated.
	x0 is a 1D m-element array that contains the first value of every
	dimension in actual units. step_sizes is a 1D m-element array that
	contains the step size between pixels in each dimension.
	axes is the list of axes (length m) to interpolate on. If m==n, then
	the estimated value of arr at point pnt is returned else, an interpolated
	array with shape of the axes not in 'axes' is returned.'''
	
	# Check if axes is a subset of all the axes.
	if len(axes) == arr.ndim:
		return interpolate_one_step(arr,pnt,x0,step_sizes)
	else:
		# Get the maximum possible index in each dimension in axes.
		max_inds = np.array(arr.shape) - 1
		max_needed_inds = []
		for ax in axes:
			max_needed_inds.append(max_inds[ax])
		max_inds = np.array(max_needed_inds)
		# Get the points we want to look at and their associated weights.
		points, weights = interpolate_weights(axes, pnt, x0, step_sizes, max_inds)
	
		all_slices = []
		# Since a point no longer represent one node, we have to make
		# the array it represents, assign the weights to each of these
		# arrays, then sum up these arrays for the final answer.
		for node in points:
			# Get lists of the indices of each axis not in axes.
			indices_lists = []
			sub_arr_shape = []
			for dim in range(arr.ndim):
				if dim not in axes:
					indices = range(arr.shape[dim])
					indices_lists.append(indices)
					sub_arr_shape.append(arr.shape[dim])
			# Now put in the indices of node in the proper place
			# based on its axis.
			for j in range(len(axes)):
				indices_lists.insert(axes[j], [node[j]])
			# Get the slice of the array we want to return. This contains
			# the indices we need to look at, not the values yet. This also
			# is flattened (1D) but can be made the right shape
			# with sub_arr_shape.
			slice = cross(indices_lists)
			# Get the values at each index. Note we must use the "c"
			# business since some of our node may be out of bounds.
			values_slice = []
			for index in slice:
				values_slice.append(get_value(arr,index))
			# Put out slice with the proper values into the proper shape.
			values_slice = np.array(values_slice)
			values_slice.shape = tuple(sub_arr_shape)
			all_slices.append(values_slice)
			
		# Now that we have all the slices (and their weights from before),
		# Multiply by their weights and add.
		interpol_sum = 0.
		for i in range(len(all_slices)):
			interpol_sum += all_slices[i] * weights[i]
		return interpol_sum
	
	
def interpolate_one_step(arr, pnt, x0, step_sizes):
	'''Return cubic convolution interpolation at n-coordinate point pnt.
	Note pnt is in the actual units of the axes. arr is the n-dimentional
	array that contains the known discrete values to be interpolated.
	x0 is a 1D n-element array that contains the first value of every
	dimension in actual units. step_sizes is a 1D n-element array that
	contains the step size between pixels in each dimension.'''
	
	# Get one of the nodes that surround the point we wish to interpolate.
	left_node = left_interpolation_node(pnt, x0, step_sizes)
	# Get the list of the 4 closest nodes (in each dimension) around the
	# interpolation point.
	needed_nodes = cross(point_list(left_node))
	# Do the n-dimensional equivalent to the formula on page 1157.
	total = 0.
	for node in needed_nodes:
		c_node = get_value(arr,node)
		u_node = u_product(pnt, node, x0, step_sizes)
		total += c_node * u_node
	return total

	
if __name__ == '__main__':
	
	# Make sure the kernel is right.
	x = np.arange(-3,3,0.001);
	ux = [];
	for i in x:
		ux.append(u(i));
	plt.plot(x,ux);
	plt.plot(x,np.zeros(len(x)),'k');
	plt.xlabel('s (pixel distance)'); plt.ylabel('kernel');
	plt.title('Interpolation Kernel. Figure 1.(d)');
	plt.show();

	## Some tests.
	print '\nx**2 with x from -5 to 5'
	print 'Point ..... Expected ..... Interpolation'
	
	#
	# 1D quadratic function.
	arr = np.array([25,16,9,4,1,0,1,4,9,16,25]);
	x0 = np.array([-5]);
	step_sizes = np.array([1]);
	
	# Some obvious points.
	pnt = np.array([0]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '0 ..... 0 ..... ' + str(ans)

	pnt = np.array([-1]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '-1 ..... 1 ..... ' + str(ans)

	pnt = np.array([3]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '3 ..... 9 ..... ' + str(ans)	
	
	# Some real tests.
	pnt = np.array([0.5]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '0.5 ..... 0.25 ..... ' + str(ans)	

	pnt = np.array([2.5]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '2.5 ..... 6.25 ..... ' + str(ans)	

	pnt = np.array([-3.456]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '-3.456 ..... 11.943936 ..... ' + str(ans)
	
	pnt = np.array([4.34]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '4.34 ..... 18.8356 ..... ' + str(ans)
	
	pnt = np.array([-4.78]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '-4.78 ..... 22.8484 ..... ' + str(ans)
	
	
	print '\nx**3 with x from -5 to 5'
	print 'Point ..... Expected ..... Interpolation'

	#
	# 1D cubic function.
	arr = np.array([125,64,27,16,1,0,1,16,27,64,125]);
	x0 = np.array([-5]);
	step_sizes = np.array([1]);
	
	# Some obvious points.
	pnt = np.array([0]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '0 ..... 0 ..... ' + str(ans)

	pnt = np.array([-3]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '-3 ..... 27 ..... ' + str(ans)	

	# Real tests. Watch how the approximation is not as good for non quadratics
	pnt = np.array([0.9]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '0.9 ..... 0.729 ..... ' + str(ans)	

	pnt = np.array([1.17]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '1.17 ..... 1.601613 ..... ' + str(ans)	

	pnt = np.array([3.65]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '3.65 ..... 48.627125 ..... ' + str(ans)	

	pnt = np.array([-4.5]);
	ans = interpolate([0], arr, pnt, x0, step_sizes);
	print '-4.5 ..... 91.125 ..... ' + str(ans)	

	
	def f(x,y):
		return x**2 + 4.*x*y + 6.*y**2

	print '\nx**2 + 4xy + 6y**2 with x from -5 to 5, y from -5 to 5'
	print 'Point ..... Expected ..... Interpolation'

	#
	# 2D quadratic function.
	# mesgrid does what I want backwards.
	y,x = np.meshgrid(np.arange(-5,5),np.arange(-5,5));
	arr = f(x,y);
	x0 = np.array([-5,-5]);
	step_sizes = np.array([1,1]);

	# Simple cases
	pnt = np.array([2,-1]);
	ans = interpolate([0,1], arr, pnt, x0, step_sizes);
	print '(2,-1) ..... 2 ..... ' + str(ans)	

	pnt = np.array([-3,4]);
	ans = interpolate([0,1], arr, pnt, x0, step_sizes);
	print '(-3,4) ..... 57 ..... ' + str(ans)	
	
	# Real tests.
	pnt = np.array([3.7,-2.3]);
	ans = interpolate([0,1], arr, pnt, x0, step_sizes);
	print '(3.7,-2.3) ..... 11.39 ..... ' + str(ans)	

	pnt = np.array([4.001,3.999]);
	ans = interpolate([0,1], arr, pnt, x0, step_sizes);
	print '(4.001,3.999) ..... 175.960003 ..... ' + str(ans)	

	pnt = np.array([4.75,-4.23]);
	ans = interpolate([0,1], arr, pnt, x0, step_sizes);
	print '(4.75,-4.23) ..... 49.5499 ..... ' + str(ans)	
	
	#
	# Interpolation of sub axes.
	print '\nSame setup as before but doing interpolation on one axes.'
	print 'Using first axis and point 3.7'
	print interpolate([0],arr,np.array([3.7]),np.array([-5]),np.array([1]));
	
	print '\nCheck by doing full interpolation on points (3.7,x)'
	print 'where x are all the values from the second axis.'
	for i in range(-5,5):
		print interpolate([0,1],arr,np.array([3.7,i]),x0,step_sizes);

	# x_large = np.arange(0,1,0.01);
	# f_large = x_large;

	# f = np.arange(0,1,0.05);
	# plt.plot(f_large,f_large,'rx');
	# plt.plot(f,f,'g+');
	# plt.show();
