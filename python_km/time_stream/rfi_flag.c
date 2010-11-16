
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define POINTS 2048
#define big_number 100

double absl(double x)  {  return x < 0. ? -x : x;  }
double sq(double x)  {  return x*x;   }

void get_rms(double *array, int len, double *rms, double *mean)    {

	/*  This function calculates the mean and rms of  an array of data points.
	 array[] is the data, len is the length of the array, rms and mean are the output variables.  
	 For a 2048 long array of XX[], you would call get_rms(XX, 2048, &rms, &mean)     */
	
	double lim,rms1 = 0.,count,mu;
	int i,l;
	
/*   A naive computation of the mean including noise spikes.    */		
	mu = 0.;
	for(i=0; i<len; i++)
		mu += array[i];
	mu /= (double)(len);
/*   A similar computation of the rms.    */				
	for(i=0; i<len; ++i)
		rms1 += sq(array[i]-mu);
	rms1 = sqrt(rms1/len);
	lim = rms1; 
	
	for(l=0; l<big_number; l++)  {     /* never ending loop is too risky. */

		count = 0; rms1 = 0.;
		for(i=0; i<len; ++i)
		/*   Recompute rms and mean ignoring suspected noise. Set your threshold. */					
			if(sq(array[i]-mu) < sq(3.*lim))
			{  count += 1.; rms1 += sq(array[i]-mu);  }
		rms1 = sqrt(rms1/count);			
		
		count = 0; mu = 0.;
		for(i=0; i<len; ++i)
			if(sq(array[i]) < sq(3.*lim))
				{  count += 1.; mu += array[i];  }
		mu /= count;

		if( absl((rms1-lim)/lim) < 0.01 ) break;  /* When to exit. */
		lim = rms1;	
	}	

	*rms = rms1;	
	*mean = mu;
}
	
void xi_square_fit(double *freq, double *array, double *fit, int len)    {

#define PARTS 4    
	
	/*    This constructs a fit by making 'p' points by avearing over 'q' values. 
	 freq is the frequency array.
	 array is the data array.
	 fit is the output array containing the fit.
	 */     
	
	double *T,*nu,count,d_parts,d_pieces,min,offset,rms,mean,
		   x_bar,y_bar,sum_x_x,sum_x_y,a,b;	
	int i,j,pieces,parts;
	
	for(i=0; i<len; fit[i]=0,++i);
	
	parts = (int)(PARTS); d_parts = (double)(PARTS);
	pieces = len/parts; d_pieces = (double)(pieces);

	T = (double *)(malloc(sizeof(double)*parts));
	nu = (double *)(malloc(sizeof(double)*parts));
	
	get_rms(array,len,&rms,&mean);
	for(i=0; i<parts; i++)   {
		
		nu[i] = 0.; T[i] = 0.; count = 0.;
		for(j=0; j<pieces; j++)  {
			
			if(absl(array[j + (i*pieces)]-mean) > 5.*rms) continue;	 /* your noise threshold here. */
			nu[i] +=  (freq[j + (i*pieces)]);
			T[i]  += (array[j + (i*pieces)]);
			count += 1.;
		}
		nu[i] /= count; T[i] /= count;
	}
	
	min = T[0];
	for(i=1; i<parts; i++)
		if(T[i] < min) min = T[i];
	offset = ( min < 0. ? 1.-min : 0. ); 
	for(i=0; i<parts; T[i]+=offset,i++); /* Don't take log of a negative number! */
	
	x_bar = 0.;  y_bar = 0.;  sum_x_x = 0.; sum_x_y = 0.;
	for(i=0; i<parts; i++)    {     /* Fit a straight line to a log-log plot */		
		
		x_bar += (log(nu[i])/d_parts); y_bar += (log(T[i])/d_parts);		
		sum_x_x += (log(nu[i])*log(nu[i])); sum_x_y += (log(nu[i])*log(T[i]));		
	}
	
	b = (sum_x_y - (d_parts*x_bar*y_bar)) / (sum_x_x - (d_parts*sq(x_bar)));
	a = y_bar - (b*x_bar);
		
	for(i=0; i<len; (fit[i] = exp(a + b*log(freq[i]))), i++);
}

void get_fit(double *array, double *f, double *fit_array)    {
	
	/*         This flattens the waveform stored in variable array[], by means of the fit. 
		  The fit is obtained and then removed from the data. */
	
	#define  NUM_PIECES 32
	#define  BLOCK 64
	
	int num_pieces = (int)(NUM_PIECES), block = (int)(BLOCK),i,j;
	double pieces[NUM_PIECES][BLOCK], f_pieces[NUM_PIECES][BLOCK], fit[BLOCK],rms,mean;
	
	for(i=0; i<num_pieces; i++)
		for(j=0; j<block; j++)    {
			
			pieces[i][j] = array[j + (i*block)];
			f_pieces[i][j] = f[j + (i*block)];
		}	
	
	for(i=0; i<num_pieces; i++)   {
		
		xi_square_fit(f_pieces[i], pieces[i], fit, block);
		for(j=0; j<block; fit_array[j+(block*i)] = fit[j],j++);			
	}
}


void flatten(double *array, double *fit_array)    {

	int i;
	double rms,mean;
	
	for(i=0; i<POINTS; array[i] -= fit_array[i],i++);
	get_rms(array,POINTS,&rms,&mean);
	for(i=0; i<POINTS; (array[i] -= mean), i++);		
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

void hi_f_spikes(double *array, double *f, int lim, int count, int tol)      {
	
	int i,j,temp_mask[POINTS];
	
	for(i=0; i<lim; temp_mask[i]=0,i++);
	rfi_find_dTdf(f,array,temp_mask,count,lim);	
	for(i=0; i<lim; i++)
		for(j=temp_mask[i]-tol; j <= temp_mask[i]+tol; array[j++] = 0.);
}

void rfi_flag(int flat, int spike, int tol, double sig, double *fit, double *array, double *f, int *mask)  {

	int i;
	double rms,mean;

	if(flat)   flatten(array,fit);  /* Compensate for structure in the cross correlation  */				
	for(i=0; i<POINTS; mask[i++]=0);				
	get_rms(array,POINTS,&rms,&mean);		
	mask_array(POINTS,array,sig*rms,mask,tol);	
}	
	

void clean(double sig, int tol, int flat, int spike, int dTdf_limit, int dTdf_tol, double *fit, double *cross1, double *array1, double *f1, int *m)  { 

	int i,ind[POINTS],lim,ct,mask[POINTS],num_entries,temp_ind[POINTS];
	double temp_f[POINTS], temp_arr[POINTS],
	       cross[POINTS], array[POINTS], f[POINTS];
		
	for(i=0; i<POINTS; ind[i]=i, i++);
	for(i=0; i<POINTS; (cross[i]=cross1[i],array[i]=array1[i],f[i]=f1[i]),i++); /* Make a copy. */
	
	rfi_flag(flat,spike,tol,sig,fit,cross,f,mask);		
	ct=0;
	for(i=0; i<POINTS; i++)
		if(mask[i] == 0)            /* Collect unflagged data only. */
		{   temp_f[ct] = f[i]; temp_arr[ct] = array[i]; temp_ind[ct] = ind[i]; ++ct; }
	
	if(spike) /* Do you want a further check of the data? Is so, compute dT/df. */
		hi_f_spikes(temp_arr,temp_f,dTdf_limit,ct,dTdf_tol);	

	for(i=0; i<POINTS;m[i]=1,i++);	/* Initialize: All points noisy. */
	for(i=0; i<POINTS; (f[i]=0.,array[i]=0.,ind[i]=0),i++);			
	num_entries=0;
	for(i=0; i<ct; i++)
		if(temp_arr[i]!=0.) 
		{  f[num_entries] = temp_f[i]; array[num_entries] = temp_arr[i]; ind[num_entries]=temp_ind[i]; ++num_entries;   }
	
	for(i=0; i<num_entries; m[ind[i++]]=0);	 /* mark as clean. */	
}

