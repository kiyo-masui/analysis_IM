#include "cMakePower.h"

int fillingf(FillConf *conf){
    
    double disc = ((double)conf->nbox_x)/conf->boxshape[0];
    //printf("\tdisc = %lg\n", disc);
    //printf("\tbox j = %d\n", conf->boxinf1[2]);
    //int pixl_num = 
    //    (conf->boxinf1[0] - conf->boxinf0[0]) *
    //    (conf->boxinf1[1] - conf->boxinf0[1]) *
    //    (conf->boxinf1[2] - conf->boxinf0[2]) ;
    //double *box_cont = (double*)malloc( pixl_num * sizeof(double));
    //for (int i=0; i<pixl_num; i++)
    //    box_cont[i] = 0.;

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
                //conf->box[box_idx] += value;
                //if (value != 0) box_cont[box_idx] += 1.;
            }
    
    //for (int i=0; i<pixl_num; i++){
    //    if (box_cont[i] != 0.){
    //        conf->box[i] /= box_cont[i];
    //    }
    //}
    return 0;

}

int fillingf2(FillConf2 *conf){

//  for (int i=0; i<conf->ran; i++)
//      printf("%5.3f\t", conf->ra[i]);
//  printf("\n");
//  printf("%lg \n", conf->mapinf[3]);
//  for(int i=0; i<conf->mapshape[2]; i++)
//      printf("%5.3f\t", conf->map[i]);
//  printf("\n\n");

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
//          printf("%5.4e\t", alpha);
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
    
//  double MAX = -1.e10; double MIN = 1.e10;
//  double MAX2 = -1.e10; double MIN2 = 1.e10;
//  for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
//      if(conf->box[i]>MAX) MAX = conf->box[i];
//      if(conf->box[i]<MIN) MIN = conf->box[i];
//      if(conf->box2[i]>MAX2) MAX2 = conf->box2[i];
//      if(conf->box2[i]<MIN2) MIN2 = conf->box2[i];
//  }
//  printf("MAX = %lg , MIN = %lg\n", MAX, MIN);
//  printf("MAX2 = %lg , MIN2 = %lg\n", MAX2, MIN2);
//  printf("Veff = %lg\n", Veff);
//  int N0 = 0, N1 = 0, N2 = 0, N3 = 0;
//  for(int i=0; i<conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]; i++){
//      if(complete[i]>1) N0++;
//      if(complete[i]==1) N1++;
//      if(complete[i]<1 && complete[i]>0){
//          N2++;
//          //conf->box[i] += (1.-complete[1])*1.e19;
//      }
//      if(complete[i]<=0) N3++;
//  }
//  printf("\n");
//  printf("    More: %d/%d\n", N0, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);
//  printf("Complete: %d/%d\n", N1, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);
//  printf("    Part: %d/%d\n", N2, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);
//  printf("   Empty: %d/%d\n", N3, conf->boxshape[0]*conf->boxshape[1]*conf->boxshape[2]);

    free(complete);

    return 0;

}

int nfillingf(FillConf2 *conf){

//  for (int i=0; i<conf->ran; i++)
//      printf("%5.3f\t", conf->ra[i]);
//  printf("\n");
//  printf("%lg \n", conf->mapinf[3]);
//  for(int i=0; i<conf->mapshape[2]; i++)
//      printf("%5.3f\t", conf->map[i]);
//  printf("\n\n");

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
//          printf("%5.4e\t", alpha);
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

    // Get the global fft information
    int Nx = fft->sizex;
    int Ny = fft->sizey;
    int Nz = fft->sizez;
    double kunitx = 1./Nx * pk->kunit;
    double kunity = 1./Ny * pk->kunit;
    double kunitz = 1./Nz * pk->kunit;

    // Initialize the 3D power spectrum
    double kmin = pk->k[0];
    double kmax = pk->k[pk->N];
    double dk = pow(10, log10(kmax/kmin)/pk->N);
    #ifdef Linearkbin
        printf("\t::Using Linear k bin\n");
        dk = (kmax-kmin)/pk->N;
    #endif
    double *kn = (double *)malloc(pk->N*sizeof(double));
    for(int i=0; i<pk->N; i++){
        kn[i] = 0.;
        pk->val[i] = 0.;
    }

    // Initialize the 2D power spectrum
    double dkp = pow(10, log10(kmax/kmin)/pk->Np);
    double dkv = pow(10, log10(kmax/kmin)/pk->Nv);
    double **kn2 = (double **)malloc(pk->Np*sizeof(double*));
    for(int i=0; i<pk->Np; i++){
        kn2[i] = (double *)malloc(pk->Nv*sizeof(double));
        for(int j=0; j<pk->Nv; j++){
            kn2[i][j] = 0.;
            pk->val2[i*pk->Nv+j] = 0.;
        }
    }

    // Function for 3D power 
    int pkplus(double k, double val, double *p, double *pn){
        if((k>=kmin)&&(k<kmax)){
            int idx = (int)(log10(k/kmin)/log10(dk));
            #ifdef Linearkbin
                idx = (int)((k-kmin)/dk);
            #endif
            p[idx] = p[idx] + val*pow(k, 3.)/2./3.1415926/3.1415926;
            pn[idx] = pn[idx] + 1.;
        }
        return 0;
    }

    // Function for 2D power 
    int pk2plus(double kp, double kv, double val, double *p, double **pn){
        if((kp>=kmin)&&(kp<kmax)&&(kv>=kmin)&&(kv<kmax)){
            double k = sqrt(kp*kp + kv*kv);
            int idxp = (int)(log10(kp/kmin)/log10(dkp));
            int idxv = (int)(log10(kv/kmin)/log10(dkv));
            int idx  = idxp*pk->Nv + idxv;
            p[idx] = p[idx] + val*pow(k, 3.)/2./3.1415926/3.1415926;
            pn[idxp][idxv] = pn[idxp][idxv] + 1.;
        }
        return 0;
    }

    double nyquist = 0.5;

    for(int i=1; i<Nx*Ny*Nz; i++){
        int xx = (int)(i/(Ny*Nz));
        int yy = (int)((i-xx*Ny*Nz)/Nz);
        int zz = (int)(i-xx*Ny*Nz-yy*Nz);
        double result0, result0p, result0v, result1;

        result1 = fft->data[i];

        double x = xx;
        double y = yy;
        double z = zz;

        // Quadrant 1
        if(x<nyquist*Nx && y<nyquist*Ny && z<nyquist*Nz){
            result0 = sqrt(x*x*kunitx*kunitx
                +y*y*kunity*kunity+z*z*kunitz*kunitz);
            pkplus(result0, result1, pk->val, kn);

            result0p = sqrt(x*x*kunitx*kunitx);
            result0v = sqrt(y*y*kunity*kunity+z*z*kunitz*kunitz);
            pk2plus(result0p, result0v, result1, pk->val2, kn2);
        }

        // Quadrant 2 
        if(x<nyquist*Nx && y>=(1.-nyquist)*Ny && z<nyquist*Nz){
            result0 = sqrt(x*x*kunitx*kunitx
                +(y-Ny)*(y-Ny)*kunity*kunity+z*z*kunitz*kunitz);
            pkplus(result0, result1, pk->val, kn);

            result0p = sqrt(x*x*kunitx*kunitx);
            result0v = sqrt((y-Ny)*(y-Ny)*kunity*kunity+z*z*kunitz*kunitz);
            pk2plus(result0p, result0v, result1, pk->val2, kn2);
        }

        // Quadrant 3 
        if(x>=(1.-nyquist)*Nx && y<nyquist*Ny && z<nyquist*Nz){
            result0 = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx
                +y*y*kunity*kunity+z*z*kunitz*kunitz);
            pkplus(result0, result1, pk->val, kn);

            result0p = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx);
            result0v = sqrt(y*y*kunity*kunity+z*z*kunitz*kunitz);
            pk2plus(result0p, result0v, result1, pk->val2, kn2);
        }

        // Quadrant 4 
        if(x>=(1.-nyquist)*Nx && y>=(1.-nyquist)*Ny && z<nyquist*Nz){
            result0 = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx
                +(y-Ny)*(y-Ny)*kunity*kunity+z*z*kunitz*kunitz);
            pkplus(result0, result1, pk->val, kn);

            result0p = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx);
            result0v = sqrt((y-Ny)*(y-Ny)*kunity*kunity+z*z*kunitz*kunitz);
            pk2plus(result0p, result0v, result1, pk->val2, kn2);
        }

        // Quadrant 5
        if(x<nyquist*Nx && y<nyquist*Ny && z>=(1.-nyquist)*Nz){
            result0 = sqrt(x*x*kunitx*kunitx
                +y*y*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
            pkplus(result0, result1, pk->val, kn);

            result0p = sqrt(x*x*kunitx*kunitx);
            result0v = sqrt(y*y*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
            pk2plus(result0p, result0v, result1, pk->val2, kn2);
        }

        // Quadrant 6
        if(x<nyquist*Nx && y>=(1.-nyquist)*Ny && z>=(1.-nyquist)*Nz){
            result0 = sqrt(x*x*kunitx*kunitx+
                (y-Ny)*(y-Ny)*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
            pkplus(result0, result1, pk->val, kn);

            result0p = sqrt(x*x*kunitx*kunitx);
            result0v = sqrt((y-Ny)*(y-Ny)*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
            pk2plus(result0p, result0v, result1, pk->val2, kn2);
        }

        // Quadrant 7
        if(x>=(1.-nyquist)*Nx && y<nyquist*Ny && z>=(1.-nyquist)*Nz){
            result0 = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx+
                y*y*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
            pkplus(result0, result1, pk->val, kn);

            result0p = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx);
            result0v = sqrt(y*y*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
            pk2plus(result0p, result0v, result1, pk->val2, kn2);
        }

        // Quadrant 8
        if(x>=(1.-nyquist)*Nx && y>=(1.-nyquist)*Ny && z>=(1.-nyquist)*Nz){
            result0 = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx+
                (y-Ny)*(y-Ny)*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
            pkplus(result0, result1, pk->val, kn);

            result0p = sqrt((x-Nx)*(x-Nx)*kunitx*kunitx);
            result0v = sqrt((y-Ny)*(y-Ny)*kunity*kunity+(z-Nz)*(z-Nz)*kunitz*kunitz);
            pk2plus(result0p, result0v, result1, pk->val2, kn2);
        }
    }

    for(int i=0; i<pk->N; i++){
        if(kn[i]!=0) 
            pk->val[i] /= kn[i];
            pk->kn[i] = kn[i];
    }
    for(int i=0; i<pk->Np; i++){
        for(int j=0; j<pk->Nv; j++){
            if(kn2[i][j]!=0){
                pk->val2[i*pk->Nv+j] /= kn2[i][j];
                pk->kn2[i*pk->Nv+j] = kn2[i][j];
            }
        }
    }

    free(kn);
    free(kn2);
    return 0;
}
