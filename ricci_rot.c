/*
    Ricci_rot: Visualization of Ricci Flow of Surfaces of Revolution
    Copyright (C) 2004  Robert M. Sinclair
    This emscripten/THREE.js adaptation Copyright (C) 2018 Anthony Carapetis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <emscripten.h>

enum N { N = 801 };
double c5 = 0, 
       c3 = 0,      // Shape parameters for initial data
       h[N], m[N],  // Metric components; i.e. h,m from the paper
       sm[N],       // Square root of m
       x[N], y[N],  // Coordinates of curve to be rotated
       oh[N], om[N], osm[N], // Old h, m, sm
       nh[N], nm[N],         // New h, m
       tmp[N], sh[N], 
       tm = 0.0;

int still_ok = 1;
int metric = 0;
int do_deriv_adj = 1;
int n_fft = 50;

double para_i[16], para_l[16], para_io[16], para_lo[16];

EMSCRIPTEN_KEEPALIVE
int resolution() { return N; }

EMSCRIPTEN_KEEPALIVE
double* coord_x() { return x; }

EMSCRIPTEN_KEEPALIVE
double* coord_y() { return y; }

void init() {
    double t, v;
    int i;

    for (i=0; i<N; i++) {
        t=((i-1.0)*M_PI)/(N-3.0);
        h[i]=1.0;
        v=(sin(t)+c5*sin(5.0*t)+c3*sin(3.0*t))
            /(1.0+3.0*c3+5.0*c5);
        if (i>20 && i<N-21 && fabs(v)<0.05) m[i]=-1.0;
        else if (v<0.0) m[i]=-1.0;
        else m[i]=pow(v,2.0);
    }
    m[1]=m[N-2]=0.0;
    m[0]=m[2];
    m[N-1]=m[N-3];

    for (i=0; i<N; i++) sm[i]=sqrt(m[i]);

    for (i=0; i<16; i++) {
        para_i[i]=(i*N+N/2)/16;
        para_l[i]=-1.0;
    }
}

int make_xy();

EMSCRIPTEN_KEEPALIVE
int reset_shape(double cc3, double cc5) {
    double oc3 = c3;
    double oc5 = c5;
    c3 = cc3;
    c5 = cc5;
    init();
    int ret = make_xy();
    if (!ret) {
        c3 = oc3;
        c5 = oc5;
        init();
        make_xy();
    }
    return ret;
}

void filter() {
    double ftt, deriv, av, c, s, dc, ds, tc;
    int i, k;

    for (i=0; i<N; i++) tmp[i]=0.0;

    av=0.0;
    for (k=0; k<n_fft; k+=2) {
        ftt=0.5*(h[1]*cos(k*((1-1.0)*M_PI)/(N-3.0))
                +h[(N-1)/2]*cos(k*(((N-1)/2-1.0)*M_PI)/(N-3.0)));
        c=cos(k*((2-1.0)*M_PI)/(N-3.0));
        s=sin(k*((2-1.0)*M_PI)/(N-3.0));
        dc=cos(k*((1.0)*M_PI)/(N-3.0));
        ds=sin(k*((1.0)*M_PI)/(N-3.0));
        for (i=2; i<(N-1)/2; i++) {
            ftt+=h[i]*c;
            tc=+dc*c-ds*s;
            s=+ds*c+dc*s;
            c=tc;
        }
        if (k==0) ftt/=0.5*(N-3.0);
        else ftt*=4.0/(N-3.0);
        av+=fabs(ftt*k);
        if (k>n_fft/2) {
            if (ftt*k/av/h[1]>0.1 && k>n_fft-20 && tm>0.001) still_ok=0;
            ftt*=1.0-sqrt((k-n_fft/2)/(0.5*n_fft+0.1));
            if (ftt>0.001) ftt*=0.1;
            if (ftt>0.01) ftt*=0.1;
        }
        for (i=0; i<=(N-1)/2; i++) tmp[i]+=ftt*cos(k*((i-1.0)*M_PI)/(N-3.0));
        for (i=(N+1)/2; i<N; i++) tmp[i]=tmp[N-1-i];
    }
    for (i=0; i<N; i++) { 
        if (tmp[i]>0.000001) h[i]=tmp[i]; 
        else h[i]=0.000001; 
    }

    deriv=0.0;
    for (i=0; i<N; i++) tmp[i]=0.0;
    for (k=1; k<n_fft; k+=2) {
        ftt=0.5*(sm[1]*sin(k*((1-1.0)*M_PI)/(N-3.0))
                +sm[(N-1)/2]*sin(k*(((N-1)/2-1.0)*M_PI)/(N-3.0)));
        c=cos(k*((2-1.0)*M_PI)/(N-3.0));
        s=sin(k*((2-1.0)*M_PI)/(N-3.0));
        dc=cos(k*((1.0)*M_PI)/(N-3.0));
        ds=sin(k*((1.0)*M_PI)/(N-3.0));
        for (i=2; i<(N-1)/2; i++) {
            ftt+=sm[i]*s;
            tc=+dc*c-ds*s;
            s=+ds*c+dc*s;
            c=tc;
        }
        ftt*=4.0/(N-3.0); if (k>n_fft/2)
        {
            ftt*=1.0-(k-n_fft/2)/(0.5*n_fft+0.1);
            if (ftt>0.001) ftt*=0.1;
            if (ftt>0.01) ftt*=0.1;
        }
        deriv+=ftt*k;
        for (i=0; i<=(N-1)/2; i++) tmp[i]+=ftt*sin(k*((i-1.0)*M_PI)/(N-3.0));
        for (i=(N+1)/2; i<N; i++) tmp[i]=tmp[N-1-i];
    }

    if (do_deriv_adj==1) ftt=sqrt(h[1])/deriv; else ftt=1.0;
    for (i=1; i<=(N-1)/2; i++) {
        sm[i]=fabs(tmp[i])*(ftt+i*(i/500.0)/n_fft)/(1.0+i*(i/500.0)/n_fft);
        m[i]=sm[i]*sm[i];
    }
    sm[0]=sm[2]; m[0]=m[2];
    for (i=(N+1)/2; i<N; i++) { sm[i]=sm[N-1-i]; m[i]=m[N-1-i]; }
}

void rescale();

EMSCRIPTEN_KEEPALIVE
void undo_step() {
    for (int i=0; i<N; i++) { h[i]=oh[i]; m[i]=om[i]; sm[i]=osm[i]; }
}

EMSCRIPTEN_KEEPALIVE
void step(double dt) {
    for (int i=0; i<N; i++) { oh[i]=h[i]; om[i]=m[i]; osm[i]=sm[i]; }
    double dm, dh, dsm, ddh, l, dm_m, hder;
    double dm1, dm2, dm3, dm4, dh1, dh2, dh3, dh4, dsm1, dsm2, dsm3, dsm4;
    int i, j, k, n;

    l=(N-3.0)/M_PI;

    for (n=0; n<1; n++) {
        for (k=0; k<4; k++) {
            for (j=0; j<1; j++) {
                ddh=(h[1+1]-2.0*h[1]+h[1-1])*l*l;
                dm_m=(2.0*m[3]-8.0*m[2])*l*l/(2.0*m[2]);
                nh[1]=h[1]-2.0*dt*(ddh/(2.0*h[1])-0.25*dm_m);
                nm[1]=0.0;

                ddh=(h[N-2+1]-2.0*h[N-2]+h[N-2-1])*l*l;
                dm_m=(2.0*m[N-4]-8.0*m[N-3])*l*l/(2.0*m[N-3]);
                nh[N-2]=h[N-2]-2.0*dt*(ddh/(2.0*h[N-2])-0.25*dm_m);
                nm[N-2]=0.0;

                for (i=2; i<N-2; i++)
                {
                    if (i==2 || i==N-3)
                    {
                        dh=(h[i+1]-h[i-1])*l*0.5;
                        dm=(m[i+1]-m[i-1])*l*0.5;
                        dsm=(sm[i+1]-2.0*sm[i]+sm[i-1])*l*l;
                    }
                    else
                    {
                        dh1=(h[i+1]-h[i-1])*l*0.5-(h[i+2]-2.0*h[i+1]+2.0*h[i-1]-h[i-2])*l*0.5/6.0;
                        dh2=(h[i+2]-h[i-2])*l*0.25-(h[i+2]-2.0*h[i+1]+2.0*h[i-1]-h[i-2])*l*0.5/6.0;
                        dh3=(h[i+1]-h[i-1])*l*0.5;
                        dh4=(h[i+2]-h[i-2])*l*0.25;
                        if (fabs(dh1)<fabs(dh2)) dh=dh1; else dh=dh2;
                        if (fabs(dh3)<fabs(dh)) dh=dh3;
                        if (fabs(dh4)<fabs(dh)) dh=dh4;
                        dm1=(m[i+1]-m[i-1])*l*0.5-(m[i+2]-2.0*m[i+1]+2.0*m[i-1]-m[i-2])*l*0.5/6.0;
                        dm2=(m[i+2]-m[i-2])*l*0.25-(m[i+2]-2.0*m[i+1]+2.0*m[i-1]-m[i-2])*l*0.5/6.0;
                        dm3=(m[i+1]-m[i-1])*l*0.5;
                        dm4=(m[i+2]-m[i-2])*l*0.25;
                        if (fabs(dm1)<fabs(dm2)) dm=dm1; else dm=dm2;
                        if (fabs(dm3)<fabs(dm)) dm=dm3;
                        if (fabs(dm4)<fabs(dm)) dm=dm4;
                        dsm1=(sm[i+1]-2.0*sm[i]+sm[i-1])*l*l-(sm[i-2]-4.0*sm[i-1]+6.0*sm[i]-4.0*sm[i+1]+sm[i+2])*l*l/8.0;
                        dsm2=(sm[i+2]-2.0*sm[i]+sm[i-2])*l*l*0.25-(sm[i-2]-4.0*sm[i-1]+6.0*sm[i]-4.0*sm[i+1]+sm[i+2])*l*l/8.0;
                        dsm3=(sm[i+1]-2.0*sm[i]+sm[i-1])*l*l;
                        dsm4=(sm[i+2]-2.0*sm[i]+sm[i-2])*l*l*0.25;
                        if (fabs(dsm1)<fabs(dsm2)) dsm=dsm1; else dsm=dsm2;
                        if (fabs(dsm3)<fabs(dsm)) dsm=dsm3;
                        if (fabs(dsm4)<fabs(dsm)) dsm=dsm4;
                    }

                    dm_m=dm/m[i];
                    hder=-dsm/sm[i]+0.25*dm_m*dh/h[i];

                    nh[i]=h[i]-2.0*dt*hder;
                    nm[i]=m[i]-2.0*dt*hder*m[i]/h[i];
                }

                nh[0]=nh[2];     nm[0]=nm[2];
                nh[N-1]=nh[N-3]; nm[N-1]=nm[N-3];

                for (i=1; i<N-1; i++)
                {
                    h[i]=nh[i];
                    m[i]=nm[i];
                }

                m[0]=m[2];     h[0]=h[2];
                m[N-1]=m[N-3]; h[N-1]=h[N-3];

                for (i=0; i<N; i++) sm[i]=sqrt(m[i]);
            }
            filter(); 
        }
        rescale(); 
    }
}

EMSCRIPTEN_KEEPALIVE
int make_xy()
{
    double l, dm, hh, dx1, dx2, dx3;
    int i;

    l=(N-3.0)/M_PI;

    x[1]=0.0;
    y[1]=0.0;

    for (i=0; i<N; i++) if (m[i]!=m[i] || h[i]!=h[i] || sm[i]!=sm[i]) return 0;

    for (i=2; i<N-1; i++)       
    {
        if (m[i]<-0.1) return 0;
        if (m[i]<0.0) y[i]=0.0;
        else y[i]=sm[i];

        dm=(sm[i]-sm[i-1])*l;
        hh=h[i-1];
        if (hh-dm*dm<0.0) dx1=0.0;
        else dx1=sqrt(hh-dm*dm)/l;

        dm=(sm[i]-sm[i-1])*l;
        hh=0.5*(h[i-1]+h[i]);
        if (hh-dm*dm<0.0) dx2=0.0;
        else dx2=sqrt(hh-dm*dm)/l;

        dm=(sm[i+1]-sm[i-1])*l*0.5;
        hh=h[i];
        if (hh-dm*dm<0.0) dx3=0.0;
        else dx3=sqrt(hh-dm*dm)/l;

        if (dx1==0.0 && dx2==0.0 && dx3==0.0) return 0;
        x[i]=x[i-1]+0.25*dx1+0.5*dx2+0.25*dx3;
    }

    for (i=1; i<N-2; i++) x[i]-=0.5*x[N-2];
    x[N-2]*=0.5;

    x[0]=x[1];     y[0]=y[1];
    x[N-1]=x[N-2]; y[N-1]=y[N-2];

    if (still_ok==0) return 0;  
    return 1;
}

void rescale()
{
    double len, a, b, c, l, s1;
    int i, ii, j;

    for (i=0; i<=(N-1)/2; i++) sh[i]=sqrt(h[i]);
    for (i=(N+1)/2; i<N; i++) sh[i]=sh[N-1-i];

    len=0.5*(sh[1]+sh[(N-1)/2]);
    for (i=2; i<(N-1)/2; i++) len+=sh[i];
    len*=2.0*M_PI/(N-3.0);

    for (i=0; i<16; i++) para_l[i]=-1.0;

    l=0.0;
    for (i=1; i<N; i++)
    {
        a=0.5*(sh[i+1]-sh[i])*M_PI/(N-3.0);
        b=sh[i]*M_PI/(N-3.0);

        for (j=0; j<16; j++)
            if (i>para_i[j]-1 && para_l[j]==-1.0)
                para_l[j]=l+a*pow(para_i[j]-i,2.0)+b*(para_i[j]-i);

        l+=a+b;
    }

    for (j=0; j<16; j++) para_i[j]=1.0+(N-3.0)*para_l[j]/len;

    tmp[1]=0.0;
    ii=2;
    l=0.0;
    for (i=1; i<N; i++)
    {
        a=0.5*(sh[i+1]-sh[i])*M_PI/(N-3.0);
        b=sh[i]*M_PI/(N-3.0);
another:
        if (ii>N-3) break;
        c=l-len*(ii-1.0)/(N-3.0);
        if (fabs(2.0*a*c/(b*b))<1e-8) s1=-c/b;
        else s1=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);

        if (s1<=1.0)
        {
            tmp[ii++]=sm[i]+s1*(sm[i+1]-sm[i]);
            goto another;
        }
        else if (i==N-3 && s1<=1.00001)
            tmp[ii++]=sm[i+1];

        l+=a+b;
    }

    tmp[N-2]=0.0;
    tmp[0]=tmp[2]; tmp[N-1]=tmp[N-3];

    for (i=0; i<N; i++) { sm[i]=tmp[i]; m[i]=tmp[i]*tmp[i]; }

    a=pow(len/M_PI,2.0);
    for (i=0; i<N; i++) h[i]=a;
}

