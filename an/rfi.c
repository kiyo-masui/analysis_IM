

/* num_pieces = 32:  Used by the get_fit() function
    big_number = 50: Used by get_rms()
    threshold = 3: Used by get_rms()
    when_to_quit = 0.1: Used by get_rms()
    part=4:  Used by xi_square_fit() 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "rfi.h"
#include "Python.h"

double absl(double x)  {  return x < 0. ? -x : x;  }
double sq(double x)  {  return x*x;   }

int get_rms(double *array, int len, double *rms, double *mean)    {

	/*  This function calculates the mean and rms of  an array of data points.
	 array[] is the data, len is the length of the array, rms and mean are the output variables.  
	 For a 2048 long array of XX[], you would call get_rms(XX, 2048, &rms, &mean)     */

	/* Threshold = About how many sigma's is teh data treated as noise? */
	/* when_to_quit:   Quit when successive computations of the RMS differ by this amount. */
	   
	double lim,rms1,count,mu,*copy,
		   threshold = 3.,
		   when_to_quit = 0.1;

	int i,l,big_number = 50,status = 0;
	
	copy = (double *)(malloc(len*sizeof(double))); 
	for(i=0; i<len; copy[i]=array[i],i++);      /* Make a copy. */
	
	for(mu=0.,i=0; i<len; mu += (copy[i]/((double)(len))), i++);
	for(i=0; i<len; copy[i] -= mu, i++);    /* subtract the mean */
		
	/*   A naive computation of the RMS including noise spikes.    */		
	for(rms1=0.,i=0; i<len; rms1 += sq(copy[i]),i++);
	rms1 = sqrt(rms1/((double)(len)));
	lim = rms1; 
	
	for(l=0; l<big_number; l++)  {     /* never ending loop is too risky. */

		/*   Recompute rms and mean ignoring suspected noise. Set your threshold. */							
		count = 0; rms1 = 0.;		
		for(i=0; i<len; ++i)
			if(sq(copy[i]) < sq(threshold*lim))
			{  count += 1.; rms1 += sq(copy[i]);  }
		rms1 = sqrt(rms1/count);			
		
		if( absl((rms1-lim)/lim) < when_to_quit ) { status = 1; break; } /* When to exit. */
		lim = rms1;	
	}	

	*rms = rms1; *mean = mu;
	return status;
}
	
void xi_square_fit(double *freq, double *array, double *fit, int len)    {
	
	/*    This constructs a fit by making 'p' points by avearing over 'q' values. 
	 I have set 'p' = 4 (see assignment 'parts=4'). 'q' is then = len/4.
	 So,  we characterize the 128 points by means of 4 representative points. A straight line is fit through these 4 points by the least-squares
	 method, and stored in the array 'fit'.    
	 
	 A power law may be better, and I will try that later. For now, it's linear.
	 
	 freq is the frequency array.
	 array is the data array.
	 fit is the output array containing the fit.
	 len is the length..
	 */     
	
	double *T,*nu,count,rms,d_parts,d_pieces,mean,
	avg,x,y,x_bar,y_bar,sum_x_x,sum_x_y,a,b;
	
	int i,j,pieces,parts;
	
	for(i=0; i<len; fit[i]=0,++i);
	
	parts = 4;	
	d_parts = (double)(parts);
	pieces = len/parts;            
	d_pieces = (double)(pieces);
	
	T = (double *)(malloc(sizeof(double)*parts));
	nu = (double *)(malloc(sizeof(double)*parts));
	
	get_rms(array,len,&rms,&mean);
	
	for(i=0; i<parts; i++)   {
		
		nu[i] = 0.;
		T[i] = 0.;
		count = 0.;
		for(j=0; j<pieces; j++)  {
			
			if(absl(array[j + (i*pieces)]-mean) > 5.*rms) continue;			
			nu[i] +=  (freq[j + (i*pieces)]);
			T[i]  += (array[j + (i*pieces)]);
			count += 1.;
		}
		nu[i] /= count;
		T[i] /= count;
	}
	
	x_bar = 0.;  y_bar = 0.;  sum_x_x = 0.; sum_x_y = 0.;
	for(i=0; i<parts; i++)    {
		
		x = nu[i]; y = T[i];		
		x_bar += (x/d_parts);y_bar += (y/d_parts);		
		sum_x_x += (x*x); sum_x_y += (x*y);		
	}
	
	b = (sum_x_y - (d_parts*x_bar*y_bar)) / (sum_x_x - (d_parts*sq(x_bar)));
	a = y_bar - (b*x_bar);
		
	for(i=0; i<len; fit[i] = a + b*freq[i],i++);	
}


void get_fit(int len, double *array, double *f, double *fit_array)    {
	
	/*         This flattens the waveform stored in variable array[], by means of the fit. 
		  The fit is obtained and then removed from the data. */
	
	int num_pieces = 32, 
		i,j,block = len/num_pieces;
	
	double *pieces, *f_pieces, *fit,rms,mean;

	pieces   = (double *)(malloc(sizeof(double)*block));
	f_pieces = (double *)(malloc(sizeof(double)*block));
	fit = (double *)(malloc(sizeof(double)*block));
	
	for(i=0; i<num_pieces; i++)   {
		
		for(j=0; j<block; j++)    
			{  pieces[j] = array[j + (i*block)]; f_pieces[j] = f[j + (i*block)];  }
					
		xi_square_fit(f_pieces,pieces, fit, block);
		for(j=0; j<block; fit_array[j+(block*i)] = fit[j],j++);			
	}
}

void get_fit_py(int len1, double *array, int len2, double *f, int len3, double *fit_array) {
  // A swig compatible wrapper for get_fit
  
  // Check that len1, len2, and len3 are all equal.
  assert(len1 == len2);
  assert(len1 == len3);

  get_fit(len1, array, f, fit_array);
}


void flatten(double *array, double *fit_array, int len)    {

	int i;
	double rms,mean;
	
	for(i=0; i<len; array[i] -= fit_array[i],i++);
	get_rms(array,len,&rms,&mean);
	for(i=0; i<len; (array[i] -= mean), i++);		
}



void mask_array(int size, double *array, double sig, int *mask, int tol)    {

	/* Prepare a mask. 1 mean noise, 0 means clean. 
	 Flag all points above a threshold, and also nearby bins. */
	
	int i,k;
	for(i=0; i<size; i++)   		
		if(absl(array[i]) > sig)    {
			
			for(k=i; k>=i-tol; --k)  mask[k] = 1; 
			for(k=i; k<=i+tol; ++k)  mask[k] = 1; 			
		}
}


void rfi_find_dTdf(double *nu, double *arr, int *mask, int count, int lim)    {

	/*   Get the noisiest XX or YY frequencies, which may have been missed. */
	
	int i,k,points,*max_pos,*max_sign;
	double *diff,max;
	
	diff = (double *)(malloc(sizeof(double)*count));	
	max_pos = (int *)(malloc(sizeof(int)*lim));	
	max_sign = (int *)(malloc(sizeof(int)*lim));		
		
	diff[count-1] = 0.;	
	for(i=0; i<count-1; ++i)   		
		diff[i] = -(arr[i+1]-arr[i])/(nu[i+1]-nu[i]);  /* derivatives dT/df. */
				
	for(k=0; k<lim; k++) {
		
		max = 0.; max_pos[k] = 0;  max_sign[k] = 1.;
		for(i=0; i<count-1; i++)
			if(absl(diff[i]) > max)   { 
				
				max = absl(diff[i]); max_pos[k] = i; 
				max_sign[k] = (diff[i] > 0. ? 1 : -1); 
			}   				
		diff[max_pos[k]] = 0.;   /* Now look for the next biggest. */
	}
		
	for(i=0; i<lim; i++)   
		mask[i] = ( max_sign[i] > 0 ? max_pos[i]+1 : max_pos[i] );	
}

void hi_f_spikes(int len, double *array, double *f, int lim, int count, int tol)      {
	
	int i,j,
		*temp_mask = (int *)(malloc(sizeof(int)*len));
	
	for(i=0; i<lim; temp_mask[i]=0,i++);
	rfi_find_dTdf(f,array,temp_mask,count,lim);	
	for(i=0; i<lim; i++)
		for(j=temp_mask[i]-tol; j <= temp_mask[i]+tol; array[j++] = 0.);
}

void rfi_flag(int len,  int flat, int spike, int tol, double sig, double *fit, double *array, double *f, int *mask)  {

	int i;
	double rms,mean;

	if(flat)   flatten(array,fit,len);  /* Compensate for structure in the cross correlation  */				
	for(i=0; i<len; mask[i++]=0);				
	get_rms(array,len,&rms,&mean);		
	mask_array(len,array,sig*rms,mask,tol);	
}	
	

void clean(int len, double sig, int tol, int flat, int spike, int dTdf_limit, int dTdf_tol, double *fit, double *cross1, double *array1, double *f1, int *m)  { 

	int i,lim,ct,num_entries,	
		*ind = (int *)(malloc(sizeof(int)*len)),
		*mask = (int *)(malloc(sizeof(int)*len)),
		*temp_ind = (int *)(malloc(sizeof(int)*len));
	
	double 		
		*temp_f = (double *)(malloc(sizeof(double)*len)),
		*temp_arr = (double *)(malloc(sizeof(double)*len)),
		*cross = (double *)(malloc(sizeof(double)*len)),
		*array = (double *)(malloc(sizeof(double)*len)),	
		*f = (double *)(malloc(sizeof(double)*len));
	
	for(i=0; i<len; ind[i]=i, i++);
	for(i=0; i<len; (cross[i]=cross1[i],array[i]=array1[i],f[i]=f1[i]),i++); /* Make a copy. */
	
	rfi_flag(len,flat,spike,tol,sig,fit,cross,f,mask);		
	ct=0;
	for(i=0; i<len; i++)
		if(mask[i] == 0)            /* Collect unflagged data only. */
		{   temp_f[ct] = f[i]; temp_arr[ct] = array[i]; temp_ind[ct] = ind[i]; ++ct; }
	
	if(spike) /* Do you want a further check of the data? Is so, compute dT/df. */
		hi_f_spikes(len,temp_arr,temp_f,dTdf_limit,ct,dTdf_tol);	

	for(i=0; i<len; m[i]=1,i++);	/* Initialize: All points noisy. */
	for(i=0; i<len; (f[i]=0.,array[i]=0.,ind[i]=0),i++);			
	num_entries=0;
	for(i=0; i<ct; i++)
		if(temp_arr[i]!=0.) 
		{  f[num_entries] = temp_f[i]; array[num_entries] = temp_arr[i]; ind[num_entries]=temp_ind[i]; ++num_entries;   }
	
	for(i=0; i<num_entries; m[ind[i++]]=0);	 /* mark as clean. */	
}

void clean_py(double sig, int tol, int flat, int spike, int dTdf_limit, int dTdf_tol, 
              int len1, double *fit, int len2, double *cross1, int len3, double *array1, 
              int len4, double *f1, int len5, int *m) {
  // Swig compatible wrapper for clean.

  // Make sure all the lengths are the same
  assert(len1 == len2);
  assert(len1 == len3);
  assert(len1 == len4);
  assert(len1 == len5);

  clean(len1, sig, tol, flat, spike, dTdf_limit, dTdf_tol, fit, cross1, array1, f1, m);
}
  
