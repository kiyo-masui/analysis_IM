#include "cMakePower.h"

int fillingf(FillConf *conf){
	
	double disc = ((double)conf->nbox_x)/conf->boxshape[0];
	//printf("\tdisc = %lg\n", disc);

	for (int i=conf->boxinf0[0]; i<conf->boxinf1[0]; i++)
		for (int j=conf->boxinf0[1]; j<conf->boxinf1[1]; j++)
			for (int k=conf->boxinf0[2]; k<conf->boxinf1[2]; k++){
				int box_idx = (int)(
					((int)(i/disc))*conf->boxshape[1]*conf->boxshape[2]+
					((int)(j/disc))*conf->boxshape[2]+
					((int)(k/disc)));
				//printf("(%d %lg)\n ",box_idx, conf->box[box_idx]);

				double box_r = sqrt(
					conf->box_x[i]*conf->box_x[i] + 
					conf->box_y[j]*conf->box_y[j] +
					conf->box_z[k]*conf->box_z[k]);
				double alpha = atan2(conf->box_y[j], conf->box_x[i]);
				double delta = asin(conf->box_z[k]/box_r);

				double ra_min = conf->ra0;
				double ra_max = conf->ra0 + conf->dra*conf->mapshape[1];
				double dec_min = conf->dec0;
				double dec_max = conf->dec0 + conf->ddec*conf->mapshape[2];

				if (alpha<ra_min || alpha>=ra_max || 
					delta<dec_min || delta>=dec_max)
					continue;

				int ra_idx = (int)((alpha-ra_min)/conf->dra);
				int dec_idx = (int)((delta-dec_min)/conf->ddec);
				int r_idx = 0;
				for (int ii=0; ii<conf->rn; ii++){
					if (box_r < conf->r[ii]){
						r_idx --;
						break;
					}
					r_idx ++;
				}
				if ((r_idx==-1)||(r_idx==conf->rn))
					continue;

				int map_idx = (int)(
					r_idx*conf->mapshape[1]*conf->mapshape[2]+
					ra_idx*conf->mapshape[2]+
					dec_idx);

				double value = conf->map[map_idx];
				conf->box[box_idx] += value*(1./pow(disc, 3.));
			}
	

	return 0;

}

int fillingf2(FillConf2 *conf){

//	for (int i=0; i<conf->ran; i++)
//		printf("%5.3f\t", conf->ra[i]);
//	printf("\n");
//	printf("%lg \n", conf->mapinf[3]);
//	for(int i=0; i<conf->mapshape[2]; i++)
//		printf("%5.3f\t", conf->map[i]);
//	printf("\n\n");

	double V = conf->boxinf[3]*conf->boxinf[3]*conf->boxinf[3];
	double *complete = (double*)malloc(
		conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]*sizeof(double));
	double Veff;
	
	for(int i=0; i<conf->rn-1; i++)
		for(int j=0; j<conf->decn; j++){
			int z = (int)
				((conf->r[i]*cos(0.5*PI-conf->dec[j])
				-conf->boxinf[2])/conf->boxinf[3]);
			double dr = conf->r[i+1] - conf->r[i];
			double v = conf->r[i]*conf->r[i]*sin(0.5*PI-conf->dec[j]);
			double alpha = v*dr*conf->mapinf[1]*conf->mapinf[2]/V;
//			printf("%5.4e\t", alpha);
			for(int k=0; k<conf->ran; k++){
				int indx = (int)(
					((int)(i/conf->mapinf[3]))*conf->mapshape[1]*conf->mapshape[2]+
					((int)(k/conf->mapinf[3]))*conf->mapshape[2]+
					((int)(j/conf->mapinf[3])));
				double value = conf->map[indx];
				double value2 = conf->map2[indx];
				if(value==0 && value2==0) continue;
				Veff = Veff + alpha*V;
				int x = (int)
					((conf->r[i]*sin(0.5*PI-conf->dec[j])*cos(conf->ra[k])
					-conf->boxinf[0])/conf->boxinf[3]);
				int y = (int)
					((conf->r[i]*sin(0.5*PI-conf->dec[j])*sin(conf->ra[k])
					-conf->boxinf[1])/conf->boxinf[3]);
				indx = x*conf->boxshape[1]*conf->boxshape[2]+y*conf->boxshape[2]+z;
				conf->box[indx] += value*alpha;
				conf->box2[indx] += value2*alpha;
				complete[indx] += alpha;

			}
		}
	
//	double MAX = -1.e10; double MIN = 1.e10;
//	double MAX2 = -1.e10; double MIN2 = 1.e10;
//	for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
//		if(conf->box[i]>MAX) MAX = conf->box[i];
//		if(conf->box[i]<MIN) MIN = conf->box[i];
//		if(conf->box2[i]>MAX2) MAX2 = conf->box2[i];
//		if(conf->box2[i]<MIN2) MIN2 = conf->box2[i];
//	}
//	printf("MAX = %lg , MIN = %lg\n", MAX, MIN);
//	printf("MAX2 = %lg , MIN2 = %lg\n", MAX2, MIN2);
//	printf("Veff = %lg\n", Veff);
//	int N0 = 0, N1 = 0, N2 = 0, N3 = 0;
//	for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
//		if(complete[i]>1) N0++;
//		if(complete[i]==1) N1++;
//		if(complete[i]<1 && complete[i]>0){
//			N2++;
//			//conf->box[i] += (1.-complete[1])*1.e19;
//		}
//		if(complete[i]<=0) N3++;
//	}
//	printf("\n");
//	printf("    More: %d/%d\n", N0, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);
//	printf("Complete: %d/%d\n", N1, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);
//	printf("    Part: %d/%d\n", N2, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);
//	printf("   Empty: %d/%d\n", N3, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);

	free(complete);

	return 0;

}

int nfillingf(FillConf2 *conf){

//	for (int i=0; i<conf->ran; i++)
//		printf("%5.3f\t", conf->ra[i]);
//	printf("\n");
//	printf("%lg \n", conf->mapinf[3]);
//	for(int i=0; i<conf->mapshape[2]; i++)
//		printf("%5.3f\t", conf->map[i]);
//	printf("\n\n");

	double V = conf->boxinf[3]*conf->boxinf[3]*conf->boxinf[3];
	double *complete = (double*)malloc(
		conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]*sizeof(double));
	for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
		complete[i] = 0.;
	}
	
	for(int i=0; i<conf->rn-1; i++)
		for(int j=0; j<conf->decn; j++){
			int z = (int)
				((conf->r[i]*cos(0.5*PI-conf->dec[j])
				-conf->boxinf[2])/conf->boxinf[3]);
			double dr = conf->r[i+1] - conf->r[i];
			double v = conf->r[i]*conf->r[i]*sin(0.5*PI-conf->dec[j]);
			double alpha = v*dr*conf->mapinf[1]*conf->mapinf[2]/V;
//			printf("%5.4e\t", alpha);
			for(int k=0; k<conf->ran; k++){
				int indx = (int)(
					((int)(i/conf->mapinf[3]))*conf->mapshape[1]*conf->mapshape[2]+
					((int)(k/conf->mapinf[3]))*conf->mapshape[2]+
					((int)(j/conf->mapinf[3])));
				double value = conf->map[indx];
				double value2 = conf->map2[indx];
				if(value==0 && value2==0) continue;
				int x = (int)
					((conf->r[i]*sin(0.5*PI-conf->dec[j])*cos(conf->ra[k])
					-conf->boxinf[0])/conf->boxinf[3]);
				int y = (int)
					((conf->r[i]*sin(0.5*PI-conf->dec[j])*sin(conf->ra[k])
					-conf->boxinf[1])/conf->boxinf[3]);
				indx = x*conf->boxshape[1]*conf->boxshape[2]+y*conf->boxshape[2]+z;
				conf->box[indx] += value*alpha*alpha;
				conf->box2[indx] += value2*alpha*alpha;
				complete[indx] += alpha;

			}
		}
	
	double MAX = -1.e10; double MIN = 1.e10;
	double MAX2 = -1.e10; double MIN2 = 1.e10;
	for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
		if(conf->box[i]>MAX) MAX = conf->box[i];
		if(conf->box[i]<MIN) MIN = conf->box[i];
		if(conf->box2[i]>MAX2) MAX2 = conf->box2[i];
		if(conf->box2[i]<MIN2) MIN2 = conf->box2[i];
	}
	//printf("MAX = %lg , MIN = %lg\n", MAX, MIN);
	//printf("MAX2 = %lg , MIN2 = %lg\n", MAX2, MIN2);
	int N1 = 0, N2 = 0, N3 = 0;
	for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
		if(complete[i]>=1) N1++;
		if(complete[i]<1 && complete[i]>0){
			N2++;
			conf->box[i] += (1.-complete[1])*1.e19;
		}
		if(complete[i]<=0) N3++;
	}
	//printf("Complete: %d/%d\n", N1, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);
	//printf("Part    : %d/%d\n", N2, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);
	//printf("Empty   : %d/%d\n", N3, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);

	free(complete);
	//return 0;

}

int makepk(FFT *fft, PK *pk){
//	printf("%lg\n\n", fft->data[1]);
//	for(int i=0; i<pk->N; i++)
//		printf("%5.4f\t",pk->val[i]);
//	printf("\n");
	int Nx = fft->sizex;
	int Ny = fft->sizey;
	int Nz = fft->sizez;
	double kunitx = 1./Nx;
	double kunity = 1./Ny;
	double kunitz = 1./Nz;

//	double kmax = sqrt(Nx*Nx*kunitx*kunitx 
//		+Ny*Ny*kunity*kunity+Nz*Nz*kunitz*kunitz);
//	double kmin = 1.*kunitx;
	double kmin = pk->k[0];
	double kmax = pk->k[pk->N-1];
	printf("%lg %lg\n", kmin, kmax);
	double dk = pow(10, log10(kmax/kmin)/pk->N);
	#ifdef Linearkbin
		printf("\t::Using Linear k bin\n");
		dk = (kmax-kmin)/pk->N;
	#endif
	double *kn = (double *)malloc(pk->N*sizeof(double));
	double dkp = pow(10, log10(kmax/kmin)/pk->Np);
	double dkv = pow(10, log10(kmax/kmin)/pk->Nv);
	double **kn2 = (double **)malloc(pk->Np*sizeof(double*));
	for(int i=0; i<pk->N; i++){
		kn[i] = 0;
	}
	for(int i=0; i<pk->Np; i++){
		kn2[i] = (double *)malloc(pk->Nv*sizeof(double));
		for(int j=0; j<pk->Nv; j++)
			kn2[i][j] = 0.;
	}
	
	
	for(int i=1; i<Nx*Ny*Nz; i++){
		int x = (int)(i/(Ny*Nz));
		int y = (int)((i-x*Ny*Nz)/Nz);
		int z = (int)(i-x*Ny*Nz-y*Nz);
		int idx, idxp, idxv;
		double result0, result1;
		result1 = fft->data[i];

		if(x<0.5*Nx && y<0.5*Ny && z<0.5*Nz){
			result0 = sqrt(x*x*kunitx*kunitx
				+y*y*kunity*kunity+z*z*kunitz*kunitz);
			if((result0>=kmin)&&(result0<kmax)){
				idx = (int)(log10(result0/kmin)/log10(dk));
				#ifdef Linearkbin
					idx = (int)((result0-kmin)/dk);
				#endif
				pk->val[idx] = pk->val[idx] + result1;
				kn[idx] = kn[idx] + 1.;
//				result0 = fabs(x*kunitx);
//				if(result0!=0){
//					idxp = (int)(log10(result0/kmin)/log10(dkp));
//					result0 = sqrt(y*y*kunity*kunity+z*z*kunitz*kunitz);
//					if(result0!=0){
//						idxv = (int)(log10(result0/kmin)/log10(dkv));
//						pk->val2[idxp*pk->Nv+idxv] = 
//							pk->val2[idxp*pk->Nv+idxv] + result1;
//						kn2[idxp][idxv] = kn2[idxp][idxv] + 1.;
//						//printf("%d %d \t", idxp, idxv);
//					}
//				}
			}
		}

		if(x<0.5*Nx && y>0.5*Ny && z<0.5*Nz){
			result0 = sqrt(x*x*kunitx*kunitx
				+(y-Ny)*(y-Ny)*kunity*kunity+z*z*kunitz*kunitz);
			if((result0>=kmin)&&(result0<kmax)){
				idx = (int)(log10(result0/kmin)/log10(dk));
				#ifdef Linearkbin
					idx = (int)((result0-kmin)/dk);
				#endif
				pk->val[idx] = pk->val[idx] + result1;
				kn[idx] = kn[idx] + 1.;
//				result0 = fabs(x*kunitx);
//				if(result0!=0){
//					idxp = (int)(log10(result0/kmin)/log10(dkp));
//					result0 = sqrt((y-Ny)*(y-Ny)*kunity*kunity+z*z*kunitz*kunitz);
//					if(result0!=0){
//						idxv = (int)(log10(result0/kmin)/log10(dkv));
//						pk->val2[idxp*pk->Nv+idxv] = 
//							pk->val2[idxp*pk->Nv+idxv] + result1;
//						kn2[idxp][idxv] = kn2[idxp][idxv] + 1.;
//					}
//				}
			}
		}

		if(x>0.5*Nx && y<0.5*Ny && z<0.5*Nz){
			result0 = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx
				+y*y*kunity*kunity+z*z*kunitz*kunitz);
			if((result0>=kmin)&&(result0<kmax)){
				idx = (int)(log10(result0/kmin)/log10(dk));
				#ifdef Linearkbin
					idx = (int)((result0-kmin)/dk);
				#endif
				pk->val[idx] = pk->val[idx] + result1;
				kn[idx] = kn[idx] + 1.;
//				result0 = fabs(x-Nx);
//				if(result0!=0){
//					idxp = (int)(log10(result0/kmin)/log10(dkp));
//					result0 = sqrt(y*y+z*z);
//					if(result0!=0){
//						idxv = (int)(log10(result0/kmin)/log10(dkv));
//						pk->val2[idxp*pk->Nv+idxv] = pk->val2[idxp*pk->Nv+idxv] + result1;
//						kn2[idxp][idxv] = kn2[idxp][idxv] + 1.;
//					}
//				}
			}
		}

		if(x>0.5*Nx && y>0.5*Ny && z<0.5*Nz){
			result0 = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx
				+(y-Ny)*(y-Ny)*kunity*kunity+z*z*kunitz*kunitz);
			if((result0>=kmin)&&(result0<kmax)){
				idx = (int)(log10(result0/kmin)/log10(dk));
				#ifdef Linearkbin
					idx = (int)((result0-kmin)/dk);
				#endif
				pk->val[idx] = pk->val[idx] + result1;
				kn[idx] = kn[idx] + 1.;
//				result0 = fabs(x-Nx);
//				if(result0!=0){
//					idxp = (int)(log10(result0/kmin)/log10(dkp));
//					result0 = sqrt((y-Ny)*(y-Ny)+z*z);
//					if(result0!=0){
//						idxv = (int)(log10(result0/kmin)/log10(dkv));
//						pk->val2[idxp*pk->Nv+idxv] = pk->val2[idxp*pk->Nv+idxv] + result1;
//						kn2[idxp][idxv] = kn2[idxp][idxv] + 1.;
//					}
//				}
			}
		}

		if(x<0.5*Nx && y<0.5*Ny && z>0.5*Nz){
			result0 = sqrt(x*x*kunitx*kunitx
				+y*y*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
			if((result0>=kmin)&&(result0<kmax)){
				idx = (int)(log10(result0/kmin)/log10(dk));
				#ifdef Linearkbin
					idx = (int)((result0-kmin)/dk);
				#endif
				pk->val[idx] = pk->val[idx] + result1;
				kn[idx] = kn[idx] + 1.;
//				result0 = fabs(x*kunitx);
//				if(result0!=0){
//					idxp = (int)(log10(result0/kmin)/log10(dkp));
//					result0 = sqrt(y*y*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
//					if(result0!=0){
//						idxv = (int)(log10(result0/kmin)/log10(dkv));
//						pk->val2[idxp*pk->Nv+idxv] = 
//							pk->val2[idxp*pk->Nv+idxv] + result1;
//						kn2[idxp][idxv] = kn2[idxp][idxv] + 1.;
//					}
//				}
			}
		}

		if(x<0.5*Nx && y>0.5*Ny && z>0.5*Nz){
			result0 = sqrt(x*x*kunitx*kunitx+
				(y-Ny)*(y-Ny)*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
			if((result0>=kmin)&&(result0<kmax)){
				idx = (int)(log10(result0/kmin)/log10(dk));
				#ifdef Linearkbin
					idx = (int)((result0-kmin)/dk);
				#endif
				pk->val[idx] = pk->val[idx] + result1;
				kn[idx] = kn[idx] + 1.;
//				result0 = fabs(x*kunitx);
//				if(result0!=0){
//					idxp = (int)(log10(result0/kmin)/log10(dkp));
//					result0 = 
//						sqrt((y-Ny)*(y-Ny)*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
//					if(result0!=0){
//						idxv = (int)(log10(result0/kmin)/log10(dkv));
//						pk->val2[idxp*pk->Nv+idxv] =
//							pk->val2[idxp*pk->Nv+idxv] + result1;
//						kn2[idxp][idxv] = kn2[idxp][idxv] + 1.;
//					}
//				}
			}
		}

		if(x>0.5*Nx && y<0.5*Ny && z>0.5*Nz){
			result0 = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx+
				y*y*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
			if((result0>=kmin)&&(result0<kmax)){
				idx = (int)(log10(result0/kmin)/log10(dk));
				#ifdef Linearkbin
					idx = (int)((result0-kmin)/dk);
				#endif
				pk->val[idx] = pk->val[idx] + result1;
				kn[idx] = kn[idx] + 1.;
//				result0 = fabs(x-Nx);
//				if(result0!=0){
//					idxp = (int)(log10(result0/kmin)/log10(dkp));
//					result0 = sqrt(y*y+(z-Nz)*(z-Nz));
//					if(result0!=0){
//						idxv = (int)(log10(result0/kmin)/log10(dkv));
//						pk->val2[idxp*pk->Nv+idxv] = pk->val2[idxp*pk->Nv+idxv] + result1;
//						kn2[idxp][idxv] = kn2[idxp][idxv] + 1.;
//					}
//				}
			}
		}

		if(x>0.5*Nx && y>0.5*Ny && z>0.5*Nz){
			result0 = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx+
				(y-Ny)*(y-Ny)*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
			if((result0>=kmin)&&(result0<kmax)){
				idx = (int)(log10(result0/kmin)/log10(dk));
				#ifdef Linearkbin
					idx = (int)((result0-kmin)/dk);
				#endif
				pk->val[idx] = pk->val[idx] + result1;
				kn[idx] = kn[idx] + 1.;
//				result0 = fabs(x-Nx);
//				if(result0!=0){
//					idxp = (int)(log10(result0/kmin)/log10(dkp));
//					result0 = sqrt((y-Ny)*(y-Ny)+(z-Nz)*(z-Nz));
//					if(result0!=0){
//						idxv = (int)(log10(result0/kmin)/log10(dkv));
//						pk->val2[idxp*pk->Nv+idxv] = pk->val2[idxp*pk->Nv+idxv] + result1;
//						kn2[idxp][idxv] = kn2[idxp][idxv] + 1.;
//					}
//				}
			}
		}
	}

	for(int i=0; i<pk->N; i++){
		if(kn[i]!=0) 
			pk->val[i] /= kn[i];
//		pk->k[i] = kmin*pow(dk, i);
//		#ifdef Linearkbin
//			pk->k[i] = kmin+dk*i;
//		#endif
	}
	for(int i=0; i<pk->Np; i++){
		pk->k2[i] = kmin*pow(dkv, i);
		pk->k2[i+pk->Np] = kmin*pow(dkp, i);
		for(int j=0; j<pk->Nv; j++)
			if(kn2[i][j]!=0)
				pk->val2[i*pk->Nv+j] /= kn2[i][j];
	}

	free(kn);
	free(kn2);
	return 0;
}
