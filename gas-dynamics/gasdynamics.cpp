#include "gasdynamics.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <omp.h>
#include <ctime>
#include <vector>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <filesystem>

////////////////////////////////////////////////
//// na podstawie programu prof. T. Chwieja ////
////////////////////////////////////////////////

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void GasDynamics::inicjalizuj() {
    generator.seed(ziarno);
    suma_czasu=0.0;
    
    nxy=nx*ny;
    delta_x=(xmax-xmin)/nx;
    delta_y=(ymax-ymin)/ny;
    
    ntot=0;
    for(int i=0;i<n_mieszanina;i++) ntot+=nc[i];
        
    czastki.resize(ntot);
    indx0.resize(nx*ny,std::vector<int>(10));
    indx.resize(ntot);
    temp_komorki.resize(nx,std::vector<double>(ny,temp));
    gestosc_komorki.resize(nx,std::vector<double>(ny,0.0));
    vx_komorki.resize(nx,std::vector<double>(ny));
    vy_komorki.resize(nx,std::vector<double>(ny));
    predkosc_x.resize(nx);
    predkosc_y.resize(ny);
    cisnienie_x.resize(nx);
    cisnienie_y.resize(ny);
    cisnienie_komorki.resize(nx,std::vector<double>(ny));
    tensor_cisnienia.resize(nx,std::vector<std::vector<double>>(ny,std::vector<double>(4)));
    
    rozklad_poczatkowy();
    sortuj();
        
    krawedz_zewn[0][0]=xmin;
    krawedz_zewn[0][1]=ymin;
    krawedz_zewn[1][0]=xmin;
    krawedz_zewn[1][1]=ymax;
    krawedz_zewn[2][0]=xmax;
    krawedz_zewn[2][1]=ymax;
    krawedz_zewn[3][0]=xmax;
    krawedz_zewn[3][1]=ymin;
    krawedz_zewn[4][0]=xmin;
    krawedz_zewn[4][1]=ymin;
}

void GasDynamics::ewolucja(double czas_max, int nit) {
    int i, it=0;
    std::string nazwa, katalog;
    
    katalog="wyniki";
    if(!std::filesystem::exists(katalog)){
        std::filesystem::create_directory(katalog);
    }
    
    FILE *fp;
    fp=fopen("xy_1.dat","w");
        
    hist_predkosci_wszystkie("hist_0.dat",5.0,50);
    zapisz_polozenie_predkosc("rv_0.dat");

    while(suma_czasu<czas_max || it<nit){
        krok();
        
        if((it%100)==0){   
            nazwa=katalog+std::string("/nptv_")+std::to_string(it)+std::string(".dat");
        }
        oblicz_srednia_droge_swobodna();
        oblicz_autokorelacje();
        wsp_dyfuzji+=wsp_cv*dt/3.0;
        
        fprintf(fp," %15.7E   %15.7E  \n",czastki[0].x,czastki[0].y);
        fflush(fp);
        it++;
    }

    fclose(fp);
}

inline void GasDynamics::sortuj() {
    int i,j,ix,iy,nsum,numer,m;
    double x,y;
    
    for(i=0;i<nx*ny;i++){
        for(j=0;j<4;j++){
            indx0[i][j]=0; 
        }
    }
    
    for(i=0;i<ntot;i++){
        x=czastki[i].x;
        y=czastki[i].y;
        if(x<xmin) czastki[i].x=xmin+delta_x*1.0E-3;
        if(x>xmax) czastki[i].x=xmax-delta_x*1.0E-3;
        if(y<ymin) czastki[i].y=ymin+delta_y*1.0E-3;
        if(y>ymax) czastki[i].y=ymax-delta_y*1.0E-3;
    }
    
    for(i=0;i<ntot;i++) indx[i]=-1;
    
    for(i=0;i<ntot;i++){
        komorka_ix_iy(czastki[i].x,czastki[i].y,&ix,&iy);
        czastki[i].ix=ix;
        czastki[i].iy=iy;
        m=globalny_indeks_komorki(ix,iy);
        indx0[m][0]++; 
    }
    
    nsum=0;
    for(m=0;m<nx*ny;m++){
        indx0[m][1]=nsum;      
        nsum+=indx0[m][0];
        indx0[m][2]=nsum-1; 
    }
    
    for(i=0;i<ntot;i++){
        ix=czastki[i].ix;
        iy=czastki[i].iy;
        if(ix>=0 && ix<nx && iy>=0 && iy<ny){
            m=globalny_indeks_komorki(ix,iy);
            indx0[m][3]++; 
            numer=indx0[m][1]+indx0[m][3]-1; 
            indx[numer]=i; 
        }
    }
}

inline void GasDynamics::komorka_ix_iy(double x, double y, int *ix, int *iy) {
    *ix=(int)((x-xmin)/delta_x);
    *iy=(int)((y-ymin)/delta_y);
}

inline int GasDynamics::globalny_indeks_komorki(int ix, int iy) {
    return ix+iy*nx;
}

inline void GasDynamics::krok() {
    int i,j,k,ix,iy,m,im,ichunk;

    double vmax=0.0;
    for(i=0;i<ntot;i++) vmax=std::max(vmax,czastki[i].v);
    dt=std::min(delta_x,delta_y)*0.999/vmax; 
    
    for(i=0;i<ntot;i++){
        czastki[i].irun=0;
    }
    
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<4;k++) tensor_cisnienia[i][j][k]=0.0;
        }
    }
    
    omp_set_num_threads(liczba_watkow);
    ichunk=3;    
    #pragma omp parallel default(shared) private(ix,iy,m,im,i)
    {
    #pragma omp for schedule(dynamic,ichunk) ordered
            for(iy=0;iy<ny;iy++){
                for(ix=0;ix<nx;ix++){
                m=globalny_indeks_komorki(ix,iy);
                for(im=indx0[m][1];im<=indx0[m][2];im++){
                    i=indx[im];
            
                    if((icol==0 || icol==1) && czastki[i].irun==0){
                        zderzenia_czastek(i,icol);
                    }
                }
            }
        }
    }  
    
    sortuj();
    oblicz_temperature_komorki();
    oblicz_cisnienie();
    oblicz_gestosc_komorki();
    oblicz_predkosc_komorki();
    suma_czasu+=dt;
}

inline int GasDynamics::zderzenia_czastek(int i, int itype_scat) {
    int i2,ix,iy,ix2,iy2,im,im2,m,m2;
    double x1,y1,x2,y2,vx1,vy1,vx2,vy2,vx_rel,vy_rel;
    double x_rel,y_rel;
    double rc1,rc2;
    double a,b,c,delta,dt0,dt1,dt2,dtmin;
    int i2_min;
    double rr_eff,rr_rel,vv_rel,rv_rel;
    int icount_1,icount_2;
    
    CZASTKA p1, p2;
    
    if(itype_scat==0 && czastki[i].irun==0){
        double dti=dt;
        int icount; 
        
        do{
            icount=odbicie_od_bariery(czastki[i],&dti); 
        }while(icount==1);
        
        if(dti>0.0) swobodny_lot(czastki[i],&dti);
        czastki[i].irun=1;
        return 0;
    }
    
    if(itype_scat>0 && czastki[i].irun==0){
        dtmin=dt*1.0; 
        i2_min=-1; 
        ix=czastki[i].ix;
        iy=czastki[i].iy;
        for(ix2=std::max(0,ix-1);ix2<=std::min(nx-1,ix+1);ix2++){
            for(iy2=std::max(0,iy-1);iy2<=std::min(ny-1,iy+1);iy2++){
                m2=ix2+iy2*nx;
                for(im2=indx0[m2][1];im2<=indx0[m2][2];im2++){
                i2=indx[im2];
                
					if(czastki[i2].irun==0 && i2!=i){
						x1=czastki[i].x;
						y1=czastki[i].y;
						vx1=czastki[i].vx;
						vy1=czastki[i].vy;
						rc1=czastki[i].rc;
								
						x2=czastki[i2].x;
						y2=czastki[i2].y;
						vx2=czastki[i2].vx;
						vy2=czastki[i2].vy;
						rc2=czastki[i2].rc;
								
						vx_rel=vx1-vx2;
						vy_rel=vy1-vy2;
								
						x_rel=x1-x2;
						y_rel=y1-y2;

						rv_rel=x_rel*vx_rel+y_rel*vy_rel;
						rr_rel=std::pow(x_rel,2)+std::pow(y_rel,2);
						vv_rel=std::pow(vx_rel,2)+std::pow(vy_rel,2);
						rr_eff=std::pow(rc1+rc2,2);
						
						a=1.0;  
						b=2.0*rv_rel/vv_rel;
						c=(rr_rel-rr_eff)/vv_rel;
						delta=b*b-4.0*a*c;
						
						if(delta>=0.0){
							dt1=(-b-std::sqrt(delta))/(2.0*a);  
							dt2=(-b+std::sqrt(delta))/(2.0*a); 
							
							if(dt1>0.0 && dt2>0.0) dt0=std::min(dt1,dt2);
							else if(dt1>0.0 && dt2<0.0) dt0=dt1;
							else if(dt1<0.0 && dt2>0.0) dt0=dt2;
							else dt0=-1.0; 
								
							if(dt0>0.0 && dt0<dtmin){
								dtmin=dt0;
								i2_min=i2;                  
							}
						}
					}
                }
            }
        }
                    
        if(i2_min>=0 && i2_min<ntot){
            p1=czastki[i];
            p2=czastki[i2_min];     
                        
            dt1=dtmin;
            dt2=dtmin;
            icount_1=odbicie_od_bariery(p1,&dt1); 
            icount_2=odbicie_od_bariery(p2,&dt2); 
                        
            if(icount_1==0 && icount_2==0){
                dt1=dtmin;
                dt2=dtmin;
                swobodny_lot(czastki[i],&dt1);
                swobodny_lot(czastki[i2_min],&dt2);      
                
                p1=czastki[i];
                p2=czastki[i2_min];
                
                rozpraszanie_czastek_kule(czastki[i],czastki[i2_min],itype_scat);                
                oblicz_tensor_cisnienia_oddzialywanie(p1,p2,czastki[i],czastki[i2_min]);
                
                dt1=dt-dtmin;         
                dt2=dt-dtmin;      
            } else {
                dt1=dt;  
                dt2=dt;             
            }
                        
            if(dt1>0.0){    
                    do{
                        icount_1=odbicie_od_bariery(czastki[i],&dt1); 
                    }while(icount_1==1);
                        
                    do{
                        icount_2=odbicie_od_bariery(czastki[i2_min],&dt2);                        
                    }while(icount_2==1);                                          
            }                                                                           
                        
            if(dt1>0.0) swobodny_lot(czastki[i],&dt1);
            if(dt2>0.0) swobodny_lot(czastki[i2_min],&dt2);     
                        
            czastki[i].irun=1;
            czastki[i2_min].irun=1;
                        
        } else if(i2_min<0) {
            dt1=dt;
            do{
                icount_1=odbicie_od_bariery(czastki[i],&dt1); 
            }while(icount_1==1);
            
            if(dt1>0.0) swobodny_lot(czastki[i],&dt1);
            czastki[i].irun=1;
        }
    }
    return 0;
}

inline void GasDynamics::rozpraszanie_czastek_kule(CZASTKA & p1, CZASTKA & p2, int itype_scat) {
    double x1,y1,x2,y2,vx1,vy1,vx2,vy2;
    double mc1,mc2;
    double t5,t6,t7,t8,t9,t10,t11,t13,t14,t15,t17,t18,t22,t35,t29,t23,t25;
    
    x1=p1.x;
    y1=p1.y;
    vx1=p1.vx;
    vy1=p1.vy;
    mc1=p1.mc;
    
    x2=p2.x;
    y2=p2.y;
    vx2=p2.vx;
    vy2=p2.vy;
    mc2=p2.mc;
    
    if(itype_scat==2){
        t5=mc1*vx1;
        t6=-x1+x2;
        t7=std::abs(t6);
        t8=t7*t7;
        t9=y1-y2;
        t10=std::abs(t9);
        t11=t10*t10;
        t13=std::sqrt(t8+t11);
        t14=1.0/t13;
        t15=-t6*t14;
        t18=t9*t14;
        t22=-t6*t14*(t15*t5+t18*mc1*vy1);
        t35=(2.0*(t22-t6*t14*(t15*mc2*vx2+t18*mc2*vy2))/(mc1+mc2)*mc1-2.0*t22+t5)/mc1;
        p1.vx=t35;
    
        t6=-x1+x2;
        t7=std::abs(t6);
        t8=t7*t7;
        t9=-y1+y2;
        t10=std::abs(t9);
        t11=t10*t10;
        t13=std::sqrt(t8+t11);
        t14=1.0/t13;
        t15=-t6*t14;
        t17=mc1*vy1;
        t18=-t9*t14;
        t22=-t9*t14*(t15*mc1*vx1+t18*t17);
        t35=(2.0*(t22-t9*t14*(t15*mc2*vx2+t18*mc2*vy2))/(mc1+mc2)*mc1-2.0*t22+t17)/mc1;
        p1.vy=t35;
    
        t6=-x1+x2;
        t7=std::abs(t6);
        t8=t7*t7;
        t9=-y1+y2;
        t10=std::abs(t9);
        t11=t10*t10;
        t13=std::sqrt(t8+t11);
        t14=1.0/t13;
        t15=-t6*t14;
        t18=-t9*t14;
        t23=mc2*vx2;
        t29=-t6*t14*(t15*t23+t18*mc2*vy2);
        t35=(2.0*(-t6*t14*(t15*mc1*vx1+t18*mc1*vy1)+t29)/(mc1+mc2)*mc2-2.0*t29+t23)/mc2;
        p2.vx=t35;
    
        t6=x1-x2;
        t7=std::abs(t6);
        t8=t7*t7;
        t9=-y1+y2;
        t10=std::abs(t9);
        t11=t10*t10;
        t13=std::sqrt(t8+t11);
        t14=1.0/t13;
        t15=t6*t14;
        t18=-t9*t14;
        t25=mc2*vy2;
        t29=-t9*t14*(t15*mc2*vx2+t18*t25);
        t35=(2.0*(-t9*t14*(t15*mc1*vx1+t18*mc1*vy1)+t29)/(mc1+mc2)*mc2-2.0*t29+t25)/mc2;
        p2.vy=t35;
    } else if(itype_scat==1) {
        double vcm_x,vcm_y,teta,vx10,vy10,vx11,vy11,vx21,vy21,ss,cc;    
        vcm_x=(mc1*vx1+mc2*vx2)/(mc1+mc2);
        vcm_y=(mc1*vy1+mc2*vy2)/(mc1+mc2);
        
        vx10=vx1-vcm_x;
        vy10=vy1-vcm_y;
        teta=losuj_U(generator)*2.0*M_PI;
        cc=std::cos(teta);
        ss=std::sin(teta);
        vx11=cc*vx10-ss*vy10;
        vy11=ss*vx10+cc*vy10;
        vx21=-vx11*mc1/mc2;
        vy21=-vy11*mc1/mc2;
        
        p1.vx=vx11+vcm_x;
        p1.vy=vy11+vcm_y;
        p2.vx=vx21+vcm_x;
        p2.vy=vy21+vcm_y;
    }
    
    p1.v=std::sqrt(std::pow(p1.vx,2)+std::pow(p1.vy,2));                
    p2.v=std::sqrt(std::pow(p2.vx,2)+std::pow(p2.vy,2));
    p1.ncol++;
    p2.ncol++;
}

inline void GasDynamics::oblicz_cisnienie() {
    int ix,iy;
    oblicz_tensor_cisnienia_kinetyczny();
    
    for(ix=0;ix<nx;ix++){
        for(iy=0;iy<ny;iy++){
            cisnienie_komorki[ix][iy]=(tensor_cisnienia[ix][iy][0]+tensor_cisnienia[ix][iy][3])/2.0/delta_x/delta_y;
        }
    }
    
    for(ix=0;ix<nx;ix++){
        cisnienie_x[ix]=0.0;
        for(iy=0;iy<ny;iy++){
            cisnienie_x[ix]+=cisnienie_komorki[ix][iy]/ny;
        }
    }
    
    for(iy=0;iy<ny;iy++){
        cisnienie_y[iy]=0.0;
        for(ix=0;ix<nx;ix++){   
            cisnienie_y[iy]+=cisnienie_komorki[ix][iy]/nx;
        }
    }
}

inline void GasDynamics::oblicz_tensor_cisnienia_kinetyczny() {
    double vx_av,vy_av,dvx,dvy,mc;
    int i,ix,iy,im,m;
    
    for(ix=0;ix<nx;ix++){
        for(iy=0;iy<ny;iy++){
            m=globalny_indeks_komorki(ix,iy);
            vx_av=0.0;
            vy_av=0.0;
            for(im=indx0[m][1];im<=indx0[m][2];im++){
                i=indx[im];
                vx_av+=czastki[i].vx/indx0[m][0];
                vy_av+=czastki[i].vy/indx0[m][0];
            }
            
            for(im=indx0[m][1];im<=indx0[m][2];im++){
                i=indx[im];
                dvx=czastki[i].vx-vx_av;
                dvy=czastki[i].vy-vy_av;
                mc=czastki[i].mc;
                tensor_cisnienia[ix][iy][0]+=mc*dvx*dvx;
                tensor_cisnienia[ix][iy][1]+=mc*dvx*dvy;
                tensor_cisnienia[ix][iy][2]+=mc*dvy*dvx;
                tensor_cisnienia[ix][iy][3]+=mc*dvy*dvy;
            }
        }
    }
}

inline void GasDynamics::oblicz_tensor_cisnienia_oddzialywanie(CZASTKA & p1_old, CZASTKA & p2_old, CZASTKA & p1_new, CZASTKA & p2_new) {
    int ix1,iy1,ix2,iy2;
    double x1,y1,x2,y2;
    double vx1,vy1,vx2,vy2;
    double mc1,mc2;
    double x12,y12,x21,y21;
    double px1,py1,px2,py2;
    
    x1=p1_old.x;
    y1=p1_old.y;
    vx1=p1_old.vx;
    vy1=p1_old.vy;
    mc1=p1_old.mc;
    
    x2=p2_old.x;
    y2=p2_old.y;
    vx2=p2_old.vx;
    vy2=p2_old.vy;
    mc2=p2_old.mc;
    
    komorka_ix_iy(x1,y1,&ix1,&iy1);
    komorka_ix_iy(x2,y2,&ix2,&iy2);
    
    x12=x1-x2;
    y12=y1-y2;
    x21=-x12;
    y21=-y12;
    
    px1=mc1*(p1_new.vx-vx1);
    py1=mc1*(p1_new.vy-vy1);
    px2=mc2*(p2_new.vx-vx2);
    py2=mc2*(p2_new.vy-vy2);
    
    tensor_cisnienia[ix1][iy1][0]+=x12*px1/2.0/dt;
    tensor_cisnienia[ix1][iy1][1]+=x12*py1/2.0/dt;
    tensor_cisnienia[ix1][iy1][2]+=y12*px1/2.0/dt;
    tensor_cisnienia[ix1][iy1][3]+=y12*py1/2.0/dt;
    
    tensor_cisnienia[ix2][iy2][0]+=x21*px2/2.0/dt;
    tensor_cisnienia[ix2][iy2][1]+=x21*py2/2.0/dt;
    tensor_cisnienia[ix2][iy2][2]+=y21*px2/2.0/dt;
    tensor_cisnienia[ix2][iy2][3]+=y21*py2/2.0/dt;
}

inline void GasDynamics::swobodny_lot(CZASTKA & par, double *dti) {
    par.x=par.x+par.vx*(*dti);
    par.y=par.y+par.vy*(*dti);
    par.path=par.path+par.v*(*dti);
    par.vxt+=par.vx*(*dti);
    par.vyt+=par.vy*(*dti);
    par.cv+=(*dti)*(par.vx0*par.vx+par.vy0*par.vy);
    *dti=0.0;
}
            
inline int GasDynamics::odbicie_od_bariery(CZASTKA & par, double *dti) {
    double x1,y1,x2,y2,vx,vy,v;
    int i,inter;
    double u1x,u1y,u2x,u2y;
    double wx,wy,zx,zy;
    double dti2,xx,yy,r;
    double cross_z;
        
    if(par.ix==0 || par.ix==(nx-1) || par.iy==0 || par.iy==(ny-1)){
        for(i=0;i<wezly_zewn;i++){
            x1=par.x;
            y1=par.y;
            vx=par.vx;
            vy=par.vy;
        
            x2=x1+vx*(*dti);
            y2=y1+vy*(*dti);
            
            u1x=krawedz_zewn[i][0];
            u1y=krawedz_zewn[i][1];
            u2x=krawedz_zewn[i+1][0];
            u2y=krawedz_zewn[i+1][1];
                
            cross_z=vx*(u2y-u1y)-vy*(u2x-u1x);
            if(cross_z>=0.0) continue;
                
            przeciecie(u1x,u1y,u2x,u2y,x1,y1,x2,y2,&wx,&wy,&zx,&zy,&inter);   
                
            if(inter==1){
                xx=wx-x1;
                yy=wy-y1;
                r=std::sqrt(xx*xx+yy*yy);
                par.path=par.path+r;
                v=par.v;
                dti2=r/v;                   
                swobodny_lot(par,&dti2); 
                *dti=*dti-dti2; 
                    
                if(tempi[i]<=1.0E-10){
                    v=std::sqrt(vx*vx+vy*vy);
                    par.vx=zx*v;
                    par.vy=zy*v;
                } else {
                    double vnew,vnew2,ux,uy,u_norm;
                    double vx_parallel,vy_parallel;
                    double vx_perpendicular,vy_perpendicular;
                                        
                    ux=u2x-u1x;
                    uy=u2y-u1y;
                    u_norm=std::sqrt(ux*ux+uy*uy);
                    ux=ux/u_norm;
                    uy=uy/u_norm;
                    
                    rozklad_maxwella(&vnew,&vnew2,tempi[i],par.mc);
                        
                    vnew=std::sqrt(vnew*vnew+vnew2*vnew2);           
                    vx_perpendicular=std::abs(vnew)*uy;
                    vy_perpendicular=-std::abs(vnew)*ux;
                        
                    rozklad_maxwella(&vnew,&vnew2,tempi[i],par.mc);
                        
                    vx_parallel=vnew*ux;
                    vy_parallel=vnew*uy;
                        
                    par.vx=vx_parallel+vx_perpendicular;
                    par.vy=vy_parallel+vy_perpendicular;
                    par.v=std::sqrt(std::pow(par.vx,2)+std::pow(par.vy,2));
                }
                par.nbound_col++; 
                return 1;
            }
        }
    }
    return 0; 
}

inline int GasDynamics::przeciecie(double u1x, double u1y, double u2x, double u2y, double v1x, double v1y, double v2x, double v2y, double *wx, double* wy, double *zx, double *zy, int *inter) {
    double u,cc,ss;
    double ux,uy,px,py;
    double v1x_p,v1y_p,v2x_p,v2y_p;
    double a,b,a_num,a_denom;

    ux=u2x-u1x;
    uy=u2y-u1y;
    u=std::sqrt(ux*ux+uy*uy);
    
    cc=ux/u; 
    ss=uy/u; 
    
    v1x_p=ss*(v1x-u1x)-cc*(v1y-u1y);
    v1y_p=cc*(v1x-u1x)+ss*(v1y-u1y);
    v2x_p=ss*(v2x-u1x)-cc*(v2y-u1y);
    v2y_p=cc*(v2x-u1x)+ss*(v2y-u1y);
    
    a_num=(v2y_p-v1y_p);
    a_denom=(v2x_p-v1x_p);
    a=a_num/a_denom; 
    b=v2y_p-a*v2x_p;
    
    if(std::abs(a_num)>std::abs(1.0E+10*a_denom)){
        *inter=0;
        *wx=0.0;
        *wy=0.0;
        return 0;
    } else if(b>=(-u*1.0E-12) && b<=u*(1.0+1.0E-12) && v1x_p*v2x_p<=0.0){
        if(b<0.0) b=0.0;
        if(b>u) b=u;
        
        *inter=1;
        *wx=cc*b+u1x;
        *wy=ss*b+u1y;
        
        ux=ss*(-1.0)*v1x_p+cc*v1y_p;  
        uy=-cc*(-1.0)*v1x_p+ss*v1y_p;
        px=ss*(-1.0)*v2x_p+cc*v2y_p;  
        py=-cc*(-1.0)*v2x_p+ss*v2y_p;
    
        ux=px-ux;
        uy=py-uy;
        
        u=std::sqrt(ux*ux+uy*uy);
        *zx=ux/u;
        *zy=uy/u;
    } else {
        *inter=0;       
        *wx=0.0;
        *wy=0.0;
    }
    return 0;
}

void GasDynamics::zapisz_polozenie_predkosc(const char* text) {
    FILE *fp;
    fp=fopen(text,"w");
    for(int i=0;i<ntot;i++){
        fprintf(fp," %15.7E  %15.7E  %15.7E  %15.7E  %15.7E  %15.7E \n  ",\
        czastki[i].x,czastki[i].y,czastki[i].vx,czastki[i].vy,czastki[i].rc,czastki[i].mc);
    }
    fclose(fp);
}

void GasDynamics::zapisz_granice_komorek(const char* text) {
    FILE *fp=fopen(text,"w");
    for(int i=0;i<=nx;i++){
        fprintf(fp,"%12.3g  %12.3g\n",delta_x*i,ymin);
        fprintf(fp,"%12.3g  %12.3g\n\n\n",delta_x*i,ymax);
    }
    for(int i=0;i<=ny;i++){
        fprintf(fp,"%12.3g  %12.3g\n",xmin,delta_y*i);
        fprintf(fp,"%12.3g  %12.3g\n\n\n",xmax,delta_y*i);
    }
    fclose(fp);
}

void GasDynamics::hist_predkosci_wszystkie(const char* text, double v_multi, int nhist) {
    double sigma,dist,vmax,dv,vx,vy,v;
    int i,j,ic;
    double **hist_v_num; 
    double **hist_v_teo;       
    
    hist_v_num=new double*[n_mieszanina];
    hist_v_teo=new double*[n_mieszanina];
    for(int i=0;i<n_mieszanina;i++){
        hist_v_num[i]=new double[nhist];
        hist_v_teo[i]=new double[nhist];
    }
    
    sigma=std::sqrt(kb*temp/mc[0]);
    vmax=v_multi*sigma;
    dv=vmax/nhist;
    
    for(ic=0;ic<n_mieszanina;ic++){
        for(i=0;i<nhist;i++){
            hist_v_num[ic][i]=0.0;
        }
    }
    
    for(i=0;i<ntot;i++){
        vx=czastki[i].vx;   
        vy=czastki[i].vy;
        v=std::sqrt(vx*vx+vy*vy);
        j=(int)(v/dv);
        ic=czastki[i].ic;
        if(j<nhist){
            hist_v_num[ic][j]+=1.0/nc[ic]/dv;
        }else{
            hist_v_num[ic][nhist-1]+=1.0/nc[ic]/dv; 
        }
    }
        
    for(ic=0;ic<n_mieszanina;ic++){
        sigma=std::sqrt(kb*temp/mc[ic]);
        for(i=0;i<nhist;i++){
            v=dv*(i+0.5);
            dist=1.0/sigma/sigma*v*std::exp(-v*v/2.0/sigma/sigma);
            hist_v_teo[ic][i]=dist;
        }
    }
    
    FILE *fp=fopen(text,"w");
    
    for(i=0;i<nhist;i++){
        v=dv*(i+0.5);
        fprintf(fp,"%12.4E ",v);
        for(ic=0;ic<n_mieszanina;ic++) fprintf(fp,"%12.4E ",hist_v_num[ic][i]);
        for(ic=0;ic<n_mieszanina;ic++) fprintf(fp,"%12.4E ",hist_v_teo[ic][i]);
        fprintf(fp,"\n");
    }
    fclose(fp);
    
    for(int i=0;i<n_mieszanina;i++) delete [] hist_v_num[i];
    delete [] hist_v_num;
    
    for(int i=0;i<n_mieszanina;i++) delete [] hist_v_teo[i];
    delete [] hist_v_teo;
}

inline void GasDynamics::zapisz_nptv(const char * filename, int msr0) {
    double temp_av,dens_av,p_av,v_av,vx_av,vy_av,jx_av,nR_av;
    int msr,ksr;
    FILE *fp;
    fp=fopen(filename,"w");
        
    msr=std::min(msr0,nx);
    
    for(int ii=0;ii<nx;ii+=msr){
        temp_av=0.0;
        dens_av=0.0;
        p_av=0.0;
        vx_av=0.0;
        vy_av=0.0;
        jx_av=0.0;
        nR_av=0.0;
        for(int ix=ii;ix<std::min(ii+msr,nx);ix++){
            ksr=std::min(ii+msr,nx)-ii;
            for(int iy=0;iy<ny;iy++){
                temp_av+=temp_komorki[ix][iy]/ksr/ny;
                dens_av+=gestosc_komorki[ix][iy]/ksr/ny;
                p_av+=cisnienie_komorki[ix][iy]/ksr/ny;
                vx_av+=vx_komorki[ix][iy]/ksr/ny;
                vy_av+=vy_komorki[ix][iy]/ksr/ny;
                jx_av+=vx_komorki[ix][iy]*gestosc_komorki[ix][iy]/ksr/ny;
            }
        }
        v_av=std::sqrt(std::pow(vx_av,2)+std::pow(vy_av,2)); 
        nR_av=p_av*delta_x*delta_y/temp_av;
        fprintf(fp,"  %12.5E",(ii+0.5*msr)*delta_x);
        fprintf(fp,"  %12.5E",dens_av);
        fprintf(fp,"  %12.5E",p_av);
        fprintf(fp,"  %12.5E",temp_av);
        fprintf(fp,"  %12.5E",v_av);
        fprintf(fp,"  %12.5E",jx_av);
        fprintf(fp,"  %12.5E",nR_av);
        fprintf(fp,"\n");           
    }
    fclose(fp);
}

inline void GasDynamics::oblicz_autokorelacje() {
    int i,k;
    wsp_cv=0.0;
    
    k=0;
    for(i=0;i<ntot;i++){
        if(czastki[i].nbound_col==0){
            wsp_cv+=czastki[i].cv;
            k++;    
        }
    }
    if(k>0) wsp_cv=wsp_cv/k/suma_czasu; 
}

inline void GasDynamics::oblicz_srednia_droge_swobodna() {
    double sigma,ro2d;
    int k=0;
    
    droga_swobodna_num=0.0;
    droga_swobodna_teo=0.0;
    
    for(int i=0;i<ntot;i++){
        if(czastki[i].ncol>0){  
            droga_swobodna_num+=czastki[i].path/czastki[i].ncol;
            k++;
        }
    }
    
    droga_swobodna_num=droga_swobodna_num/k;
    sigma=4.0*czastki[0].rc;
    ro2d=ntot/(xmax-xmin)/(ymax-ymin);  
    droga_swobodna_teo=1.0/sigma/ro2d/std::sqrt(2.0);   
}

inline void GasDynamics::rozklad_maxwella(double *vx, double *vy, double temp_i, double mci) {
    double sigma=std::sqrt(kb*temp_i/mci);
    *vx=losuj_N(generator)*sigma;
    *vy=losuj_N(generator)*sigma;
}

inline double GasDynamics::temperatura_efektywna() {
    double ekin=0.0;
    for(int i=0;i<ntot;i++){
        ekin+=czastki[i].mc*std::pow(czastki[i].v,2)/2.0;
    }
    return ekin/ntot/kb;
}

inline void GasDynamics::oblicz_temperature_komorki() {
    for(int ix=0;ix<nx;ix++){
        for(int iy=0;iy<ny;iy++){
            temp_komorki[ix][iy]=0.0;
            int m=globalny_indeks_komorki(ix,iy);
            for(int im=indx0[m][1];im<=indx0[m][2];im++){
                int i=indx[im];
                temp_komorki[ix][iy]+=czastki[i].mc*std::pow(czastki[i].v,2)/2.0;
            }
            if(indx0[m][0]>0) temp_komorki[ix][iy]=temp_komorki[ix][iy]/kb/indx0[m][0];
        }
    }
}

inline void GasDynamics::oblicz_gestosc_komorki() {
    for(int ix=0;ix<nx;ix++){
        for(int iy=0;iy<ny;iy++){
            int m=globalny_indeks_komorki(ix,iy);
            gestosc_komorki[ix][iy]=indx0[m][0]/delta_x/delta_y;
        }
    }
}

inline void GasDynamics::oblicz_predkosc_komorki() {
    for(int ix=0;ix<nx;ix++){
        for(int iy=0;iy<ny;iy++){
            vx_komorki[ix][iy]=0.0;
            vy_komorki[ix][iy]=0.0;
            int m=globalny_indeks_komorki(ix,iy);
            int ile=indx0[m][0];
            for(int im=indx0[m][1];im<=indx0[m][2];im++){
                int i=indx[im];
                vx_komorki[ix][iy]+=czastki[i].vx/ile;
                vy_komorki[ix][iy]+=czastki[i].vy/ile;
            }
        }
    }
}

inline void GasDynamics::test() {
    int n=10000;
    double mci=czastki[0].mc;
    double temp_i=tempi[0];
    double etot=0.0;
    double vx, vy;
    for(int i=0;i<n;i++){
        rozklad_maxwella(&vx,&vy,temp_i,mci);   
        etot+=mci*(vx*vx+vy*vy)/2.0;
    }
    etot=etot/n/kb;
}

void GasDynamics::rozklad_poczatkowy() {
    if(init_dist==0){
        double x,y,vx,vy,rci,mci;
        int k=0,ieof;
        
        FILE *fp_in=fopen("pos_vel_start.dat","r");
        if(fp_in==NULL){   
            printf("Brak pliku z polozeniami czastek -> KONIEC\n"); 
            exit(-1); 
        } else {
            for(int ic=0;ic<n_mieszanina;ic++){
                for(int i=0;i<nc[ic];i++){
                    ieof=fscanf(fp_in,"%lf %lf %lf  %lf %lf  %lf",&x,&y,&vx,&vy,&rci,&mci);  
                    if(ieof!=EOF){
                        czastki[k].x=x; 
                        czastki[k].y=y; 
                        czastki[k].vx=vx;   
                        czastki[k].vy=vy;
                        czastki[k].x0=x;    
                        czastki[k].y0=y;    
                        czastki[k].vx0=vx;  
                        czastki[k].vy0=vy;
                        czastki[k].v=std::sqrt(vx*vx+vy*vy); 
                        czastki[k].rc=rci;  
                        czastki[k].mc=mci;  
                        czastki[k].ic=ic;
                        czastki[k].irun=0;
                        czastki[k].ix=-1;
                        czastki[k].iy=-1;
                        k++;
                    }
                }
            }
            if(k!=ntot){
                printf("za malo danych: ntot= %d, wczytano k=%d\nEXIT\n",ntot,k);
                exit(0);
            }
        }
    } else { 
        int k=0;
        double x,y,v,vx,vy,vx0,vy0,sigma;
        
        for(int ic=0;ic<n_mieszanina;ic++){
            for(int i=0;i<nc[ic];i++){
                if(init_dist==1 || init_dist==2){
                    x=losuj_U(generator)*(xmax-xmin)+xmin;
                    y=losuj_U(generator)*(ymax-ymin)+ymin;
                }else if(init_dist==3 || init_dist==4){
                    x=losuj_U(generator)*delta_x+xmin;
                    y=losuj_U(generator)*delta_y+ymin;
                }   
                sigma=std::sqrt(kb*temp/mc[ic]);
                vx=losuj_N(generator)*sigma;
                vy=losuj_N(generator)*sigma;
                if(init_dist==1 || init_dist==3){
                    v=std::sqrt(2.0)*sigma;
                    vx0=v*vx/std::sqrt(vx*vx+vy*vy);
                    vy0=v*vy/std::sqrt(vx*vx+vy*vy);
                    vx=vx0;
                    vy=vy0;
                }
                
                czastki[k].x=x; 
                czastki[k].y=y; 
                czastki[k].vx=vx;   
                czastki[k].vy=vy;
                czastki[k].x0=x;    
                czastki[k].y0=y;    
                czastki[k].vx0=vx;  
                czastki[k].vy0=vy;
                czastki[k].v=std::sqrt(vx*vx+vy*vy); 
                czastki[k].rc=rc[ic];   
                czastki[k].mc=mc[ic];   
                czastki[k].ic=ic;
                czastki[k].irun=0;
                czastki[k].ix=-1;
                czastki[k].iy=-1;
                czastki[k].indx=k; 
                czastki[k].ncol=0;
                k++;
            }
        }
        
        double temp1=temperatura_efektywna();
        for(int i=0;i<ntot;i++){
            vx=czastki[i].vx*std::sqrt(temp/temp1);
            vy=czastki[i].vy*std::sqrt(temp/temp1);
            v=std::sqrt(vx*vx+vy*vy);
            czastki[i].vx=vx;
            czastki[i].vy=vy;
            czastki[i].v=v;
        }
    }   
}

void GasDynamics::czytaj(const char * text) {
    FILE *fp_in=fopen(text,"r");
    
    if(fp_in==NULL){   
        printf("Brak pliku -> KONIEC\n"); 
        exit(-1); 
    } else {
        fscanf(fp_in,"%lf %lf %lf  %lf",&xmin,&xmax,&ymin,&ymax);    
        fscanf(fp_in,"%d %d",&nx,&ny);         
        fscanf(fp_in,"%lf ",&kb);
        fscanf(fp_in,"%lf %lf %lf %lf %lf ",&temp,&tempi[0],&tempi[1],&tempi[2],&tempi[3]);            
        fscanf(fp_in,"%d ",&init_dist); 
        fscanf(fp_in,"%d ",&n_mieszanina); 
        for(int i=0;i<n_mieszanina;i++) 
            fscanf(fp_in,"%d  %lf %lf ",&nc[i],&mc[i],&rc[i]);

        fscanf(fp_in,"%d ",&wezly); 
        for(int i=0;i<wezly;i++) 
            fscanf(fp_in,"%lf %lf ",&krawedz[i][0],&krawedz[i][1]);
            
        krawedz[wezly][0]=krawedz[0][0];
        krawedz[wezly][1]=krawedz[0][1];
    }
    fclose(fp_in);
}