/*
  Ricci_rot: Visualization of Ricci Flow of Surfaces of Revolution
  Copyright (C) 2004  Robert M. Sinclair

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>

#define N 801
double c5, c3, h[N], m[N], sm[N], nh[N], nm[N], x[N], y[N], oh[N], om[N], osm[N],
       tmp[N], sh[N], dt, tm;
int still_ok;
int metric;
int do_deriv_adj;
int n_fft;

double para_i[16], para_l[16], para_io[16], para_lo[16];

void init()
{
  double t, v;
  int i;

  for (i=0; i<N; i++)
  {
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

  for (i=0; i<16; i++)
  {
    para_i[i]=(i*N+N/2)/16;
    para_l[i]=-1.0;
  }
}

void filter()
{
  double ftt, deriv, av, c, s, dc, ds, tc;
  int i, k;
  
  for (i=0; i<N; i++) tmp[i]=0.0;
  
  av=0.0;
  for (k=0; k<n_fft; k+=2)
  {
    ftt=0.5*(h[1]*cos(k*((1-1.0)*M_PI)/(N-3.0))
      +h[(N-1)/2]*cos(k*(((N-1)/2-1.0)*M_PI)/(N-3.0)));
    c=cos(k*((2-1.0)*M_PI)/(N-3.0));
    s=sin(k*((2-1.0)*M_PI)/(N-3.0));
    dc=cos(k*((1.0)*M_PI)/(N-3.0));
    ds=sin(k*((1.0)*M_PI)/(N-3.0));
    for (i=2; i<(N-1)/2; i++)
    {
      ftt+=h[i]*c;
      tc=+dc*c-ds*s;
       s=+ds*c+dc*s;
       c=tc;
    }
    if (k==0) ftt/=0.5*(N-3.0);
    else ftt*=4.0/(N-3.0);
    av+=fabs(ftt*k);
    if (k>n_fft/2)
    {
      if (ftt*k/av/h[1]>0.1 && k>n_fft-20 && tm>0.001) still_ok=0;
      ftt*=1.0-sqrt((k-n_fft/2)/(0.5*n_fft+0.1));
      if (ftt>0.001) ftt*=0.1;
      if (ftt>0.01) ftt*=0.1;
    }
    for (i=0; i<=(N-1)/2; i++) tmp[i]+=ftt*cos(k*((i-1.0)*M_PI)/(N-3.0));
    for (i=(N+1)/2; i<N; i++) tmp[i]=tmp[N-1-i];
  }
  for (i=0; i<N; i++) { if (tmp[i]>0.000001) h[i]=tmp[i]; else h[i]=0.000001; }

  deriv=0.0;
  for (i=0; i<N; i++) tmp[i]=0.0;
  for (k=1; k<n_fft; k+=2)
  {
    ftt=0.5*(sm[1]*sin(k*((1-1.0)*M_PI)/(N-3.0))
            +sm[(N-1)/2]*sin(k*(((N-1)/2-1.0)*M_PI)/(N-3.0)));
    c=cos(k*((2-1.0)*M_PI)/(N-3.0));
    s=sin(k*((2-1.0)*M_PI)/(N-3.0));
    dc=cos(k*((1.0)*M_PI)/(N-3.0));
    ds=sin(k*((1.0)*M_PI)/(N-3.0));
    for (i=2; i<(N-1)/2; i++)
    {
      ftt+=sm[i]*s;
      tc=+dc*c-ds*s;
       s=+ds*c+dc*s;
       c=tc;
    }
    ftt*=4.0/(N-3.0);
    if (k>n_fft/2)
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
  for (i=1; i<=(N-1)/2; i++)
  {
    sm[i]=fabs(tmp[i])*(ftt+i*(i/500.0)/n_fft)/(1.0+i*(i/500.0)/n_fft);
    m[i]=sm[i]*sm[i];
  }
  sm[0]=sm[2]; m[0]=m[2];
  for (i=(N+1)/2; i<N; i++) { sm[i]=sm[N-1-i]; m[i]=m[N-1-i]; }
}

void rescale();

void step()
{
  double dm, dh, dsm, ddh, l, dm_m, hder;
  double dm1, dm2, dm3, dm4, dh1, dh2, dh3, dh4, dsm1, dsm2, dsm3, dsm4;
  int i, j, k, n;

  l=(N-3.0)/M_PI;

  for (n=0; n<1; n++) {
  for (k=0; k<4; k++) {
  for (j=0; j<1; j++)
  {
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
  filter(); }
  rescale(); }
}

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

void make_vertex(double theta, int n, int no)
{
  GLfloat xx, yy, zz, r, c, s, n1, n2, inl;

  c=cos(theta);
  s=sin(theta);
  r=y[n];
  xx=c*r;   
  yy=s*r;
  zz=x[n];
  n1=y[n+1]-y[n-1];
  n2=x[n-1]-x[n+1];
  inl=1.0/sqrt(n1*n1+n2*n2);
  if (n2<0.0) { n2=-n2; n1=-n1; }
  if (no==1)
  {
    glNormal3f(c*n2*inl,n1*inl,s*n2*inl);  
  }
  else
  {
    xx+=0.0015*c*n2*inl;
    yy+=0.0015*s*n2*inl;
    zz+=0.0015*n1*inl;
  }
  glVertex3f((GLfloat)7.0*xx,(GLfloat)7.0*zz,(GLfloat)7.0*yy);
}

void make_vertex_d(double theta, double dn, int no)
{
  GLfloat xx, yy, zz, r, c, s, n1, n2, inl, t;
  int n;

  n=(int)floor(dn);
  c=cos(theta);
  s=sin(theta);
  r=y[n];
  xx=c*r;   
  yy=s*r;
  t=dn-n;
  zz=t*x[n+1]+(1.0-t)*x[n];
  n1=y[n+1]-y[n-1];
  n2=x[n-1]-x[n+1];
  inl=1.0/sqrt(n1*n1+n2*n2);
  if (n2<0.0) { n2=-n2; n1=-n1; }
  if (no==1)
  {
    glNormal3f(c*n2*inl,n1*inl,s*n2*inl);  
  }
  else
  {
    xx+=0.0015*c*n2*inl;
    yy+=0.0015*s*n2*inl;
    zz+=0.0015*n1*inl;
  }
  glVertex3f((GLfloat)7.0*xx,(GLfloat)7.0*zz,(GLfloat)7.0*yy);
}

int choose_surface=1;

void draw_full_surface()
{
  int i, j, N_theta=50, skip=10;
  GLfloat mdiff[]={0.8,0.8,0.8,1.0};
  GLfloat memm[]={0.05,0.05,0.05,1.0};
  GLfloat shininess[]={200.0};
  
  for (j=skip; j<N-2-skip; j+=skip)
  {
    glBegin(GL_TRIANGLE_STRIP);
      glMaterialfv(GL_FRONT,GL_DIFFUSE,mdiff);
      glMaterialfv(GL_FRONT,GL_EMISSION,memm);
      glMaterialfv(GL_FRONT,GL_SHININESS,shininess);
      
      for (i=0; i<=N_theta; i++)
      {
        make_vertex(i*2.0*M_PI/N_theta,j,1);
        make_vertex(i*2.0*M_PI/N_theta,j+skip,1);
      }
    glEnd();
  }
  
  glBegin(GL_TRIANGLE_FAN);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,mdiff);
    glMaterialfv(GL_FRONT,GL_EMISSION,memm);
    glMaterialfv(GL_FRONT,GL_SHININESS,shininess);
    
    glNormal3f((GLfloat)0.0,(GLfloat)1.0,(GLfloat)0.0);
    make_vertex(0.0,N-2,0);
    for (i=0; i<=N_theta; i++)
      make_vertex(i*2.0*M_PI/N_theta,j-skip,1);
  glEnd();
  
  glBegin(GL_TRIANGLE_FAN);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,mdiff);
    glMaterialfv(GL_FRONT,GL_EMISSION,memm);
    glMaterialfv(GL_FRONT,GL_SHININESS,shininess);
    
    glNormal3f((GLfloat)0.0,(GLfloat)-1.0,(GLfloat)0.0);
    make_vertex(0.0,1,0);
    for (i=0; i<=N_theta; i++)
      make_vertex(i*2.0*M_PI/N_theta,skip,1);
  glEnd();
  
  glDisable(GL_LIGHTING);

  for (i=0; i<=32; i++)
  {
    glBegin(GL_LINE_STRIP);
    if (choose_surface==1) glColor3f((GLfloat)0.0,(GLfloat)0.0,(GLfloat)1.0);
    else if (still_ok==1) glColor3f((GLfloat)0.0,(GLfloat)0.2,(GLfloat)0.0);
    else glColor3f((GLfloat)0.0,(GLfloat)0.0,(GLfloat)0.0);
    for (j=1; j<N-2; j+=skip)
    {
      make_vertex(i*2.0*M_PI/32,j,0);
    }
    make_vertex(i*2.0*M_PI/32,N-2,0);
    glEnd();
  }
  
  for (j=0; j<16; j++)
  {
    glBegin(GL_LINE_STRIP);
    if (choose_surface==1) glColor3f((GLfloat)0.0,(GLfloat)0.0,(GLfloat)1.0);
    else if (still_ok==1) glColor3f((GLfloat)0.0,(GLfloat)0.2,(GLfloat)0.0);
    else glColor3f((GLfloat)0.0,(GLfloat)0.0,(GLfloat)0.0);
    for (i=0; i<=N_theta; i++)
      make_vertex_d(i*2.0*M_PI/N_theta,para_i[j],0);
    glEnd();
  }
  glEnable(GL_LIGHTING);
}

void draw_surface()
{
  int i, j;
  
  if (metric==1)
  {
    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);
      glColor3f((GLfloat)0.0,(GLfloat)0.0,(GLfloat)1.0);
      for (i=0; i<N; i++)
        glVertex3f((GLfloat)(20.0*(i-0.5*N)/N),(GLfloat)(-8.0+10.0*m[i]),(GLfloat)(0.0));
      glColor3f((GLfloat)0.0,(GLfloat)1.0,(GLfloat)0.0);
      for (i=0; i<N; i++)
        glVertex3f((GLfloat)(20.0*(i-0.5*N)/N),(GLfloat)(-8.0+10.0*h[i]),(GLfloat)(0.0));

      glColor3f((GLfloat)1.0,(GLfloat)1.0,(GLfloat)1.0);
      for (i=1; i<N-1; i++)
      {
        glVertex3f((GLfloat)( 5.0*y[i]),(GLfloat)(5.0*x[i]),(GLfloat)(0.0));
        glVertex3f((GLfloat)(-5.0*y[i]),(GLfloat)(5.0*x[i]),(GLfloat)(0.0));
      }

    glEnd();
    for (j=0; j<16; j++) {
    glBegin(GL_LINE_STRIP);
      glColor3f((GLfloat)1.0,(GLfloat)1.0,(GLfloat)1.0);
      glVertex3f((GLfloat)(
          5.0*(y[(int)floor(para_i[j])]
         +(para_i[j]-floor(para_i[j]))*(y[(int)floor(para_i[j])+1]-y[(int)floor(para_i[j])]))),
                 (GLfloat)(
          5.0*(x[(int)floor(para_i[j])]
         +(para_i[j]-floor(para_i[j]))*(x[(int)floor(para_i[j])+1]-x[(int)floor(para_i[j])]))),
                 (GLfloat)(0.0));
      glVertex3f((GLfloat)(
         -5.0*(y[(int)floor(para_i[j])]
         +(para_i[j]-floor(para_i[j]))*(y[(int)floor(para_i[j])+1]-y[(int)floor(para_i[j])]))),
                 (GLfloat)(
          5.0*(x[(int)floor(para_i[j])]
         +(para_i[j]-floor(para_i[j]))*(x[(int)floor(para_i[j])+1]-x[(int)floor(para_i[j])]))),
                 (GLfloat)(0.0));
    glEnd(); }
    glEnable(GL_LIGHTING);
  }
  else
    draw_full_surface();
}

int rotate=0, movelight=0, origx, origy;
int rot_angle=0;
int spinxlight=0;
int spinylight=0;
int flow=0;
GLdouble rotation_matrix[16], rot_cos=1.0, rot_sin=0.0;
int xpos, ypos, oxpos, oypos;

struct Vector_3
{
  double x, y, z;
};

void new_Vector_3(struct Vector_3 *a, double xx, double yy, double zz)
{
  a->x=xx;
  a->y=yy;
  a->z=zz;
}

double Vector_3_dot(struct Vector_3 *a, struct Vector_3 *b)
{
  return a->x*b->x+a->y*b->y+a->z*b->z;
}

void Vector_3_normalize(struct Vector_3 *a)
{
  double il;

  il=1.0/sqrt(a->x*a->x+a->y*a->y+a->z*a->z);
  (a->x)*=il;
  (a->y)*=il;
  (a->z)*=il;
}

struct Matrix_3_3
{
  double m[3][3];
} rmat, refmat;

void Matrix_3_3_copy(struct Matrix_3_3 *a, struct Matrix_3_3 *b)
{
  int i, j;
  
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      b->m[i][j]=a->m[i][j];
}

void new_Matrix_3_3(struct Matrix_3_3 *a,
                    double m00, double m01, double m02,
                    double m10, double m11, double m12,
                    double m20, double m21, double m22)
{
  a->m[0][0]=m00; a->m[0][1]=m01; a->m[0][2]=m02;
  a->m[1][0]=m10; a->m[1][1]=m11; a->m[1][2]=m12;
  a->m[2][0]=m20; a->m[2][1]=m21; a->m[2][2]=m22;
}

void Matrix_3_3_multiply(struct Matrix_3_3 *a, struct Matrix_3_3 *b,
                         struct Matrix_3_3 *c)
{
  int i, j, k;

  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      c->m[i][j]=0.0;
      for (k=0; k<3; k++)
        c->m[i][j]+=a->m[i][k]*b->m[k][j];
    }
  }
}
  
void motion(int x, int y)
{
  struct Matrix_3_3 t1, t2, t3, t4;
  struct Vector_3 axis0, axis1, axis2;
  double c, s, angle, oc5, oc3;
  int i;
  
  if (rotate && metric!=1)
  {
    xpos=x; ypos=-y;  

    if (xpos==oxpos && ypos==oypos)
    {
      Matrix_3_3_copy(&refmat,&rmat);
      glutPostRedisplay();
      return;
    }

    angle=sqrt((double)((oypos-ypos)*(oypos-ypos))
              +(double)((oxpos-xpos)*(oxpos-xpos)))/300.0;
    c=cos(angle);
    s=sin(angle);

    new_Vector_3(&axis0,(double)( rot_cos*(oypos-ypos)-rot_sin*(xpos-oxpos)),
                        (double)(+rot_sin*(oypos-ypos)+rot_cos*(xpos-oxpos)),0.0);
    new_Vector_3(&axis1,(double)( rot_cos*(xpos-oxpos)-rot_sin*(ypos-oypos)),
                        (double)(+rot_sin*(xpos-oxpos)+rot_cos*(ypos-oypos)),0.0);
    new_Vector_3(&axis2,0.0,0.0,1.0);
    Vector_3_normalize(&axis0);
    Vector_3_normalize(&axis1);  

    new_Matrix_3_3(&t1,
        axis0.x,axis1.x,axis2.x,
        axis0.y,axis1.y,axis2.y,
        axis0.z,axis1.z,axis2.z);

    new_Matrix_3_3(&t2,1.0,0.0,0.0,
                       0.0,  c,  s,
                       0.0, -s,  c);

    new_Matrix_3_3(&t3,
        axis0.x,axis0.y,axis0.z,
        axis1.x,axis1.y,axis1.z,
        axis2.x,axis2.y,axis2.z);
    
    Matrix_3_3_multiply(&t3,&refmat,&t4);
    Matrix_3_3_multiply(&t2,&t4,&t3);
    Matrix_3_3_multiply(&t1,&t3,&rmat);
    
    glutPostRedisplay();
  }
  if (movelight)
  {
    spinylight=(spinylight+(x-origx))%720;
    spinxlight=(spinxlight+(y-origy))%720;
    origx=x;
    origy=y;
    glutPostRedisplay();
  }
  if (choose_surface)
  {
    still_ok=1;
    oc5=c5;
    oc3=c3;
    for (i=0; i<N; i++) { oh[i]=h[i]; om[i]=m[i]; osm[i]=sm[i]; }
    for (i=0; i<16; i++) { para_io[i]=para_i[i]; para_lo[i]=para_l[i]; }
    c3=(c3+0.001*(x-origx));
    c5=(c5+0.001*(y-origy));
    init();
    if (make_xy()==1) glutPostRedisplay();
    else
    {
      for (i=0; i<N; i++) { h[i]=oh[i]; m[i]=om[i]; sm[i]=osm[i]; }
      for (i=0; i<16; i++) { para_i[i]=para_io[i]; para_l[i]=para_lo[i]; }
      c5=oc5;
      c3=oc3;
      make_xy();
      glutPostRedisplay();
    }
    origx=x;
    origy=y;
  }
}

void mouse(int button, int state, int x, int y)
{
  switch(button)
  {
    case GLUT_LEFT_BUTTON:     
      if (state==GLUT_DOWN && choose_surface==0)
      {
        oxpos=x;
        oypos=-y;
        Matrix_3_3_copy(&rmat,&refmat);
        
        origx=x;
        origy=y;
        rotate=1;
      }
      else if (state==GLUT_DOWN && choose_surface==1)
      {
        origx=x;
        origy=y;
        rotate=0;
      }
      else
        rotate=0;
      break;

    case GLUT_MIDDLE_BUTTON:
      if (state==GLUT_DOWN)
      {
        origx=x;
        origy=y;
        movelight=1;
      }
      else
        movelight=0;
  }
}

int spotlight=1;

void main_menu_select(int value)
{
  switch (value)
  {
    case 21: /* Choose surface */
      choose_surface=1;
      flow=0;
      movelight=0;
      rotate=0;
      glutPostRedisplay();
      return;
    case 7: /* Flow */
      choose_surface=0;
      flow=1;
      movelight=0;
      rotate=0;
      glutPostRedisplay();
      return;
    case 999:
      exit(0);
  }
}

void spotlight_select(int value)
{
  switch (value)
  {
    case 0:
      spotlight=0;
      glutPostRedisplay();
      return;
    case 1:
      spotlight=1;
      glutPostRedisplay();
      return;
    case 999:
      exit(0);
  }
}

void make_menus(void)
{
  int mode_menu;

  mode_menu=glutCreateMenu(spotlight_select);  
  glutAddMenuEntry("Turn off", 0);
  glutAddMenuEntry("Turn on", 1);
  glutCreateMenu(main_menu_select);
  glutAddMenuEntry("New shape", 21);
  glutAddMenuEntry("Flow", 7);
  glutAddSubMenu("Spotlight", mode_menu);
  glutAddMenuEntry("Quit", 999);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void display(void)
{
  GLfloat position[]={100.0,0.0,0.0,0.0};
  GLfloat position1[]={0.0,-70.0,0.0,0.0};
  GLfloat spec[]={0.0,0.95,0.0,1.0};
  GLfloat dff1[]={0.8,0.4,0.4,1.0};

  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  glPointSize((GLfloat)1.0);

  glPushMatrix();
  
  if (metric==0)
  {
    glLightfv(GL_LIGHT1,GL_POSITION,position1);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,dff1);
    
    rotation_matrix[ 0]= rot_cos*rmat.m[0][0]+rot_sin*rmat.m[1][0];
    rotation_matrix[ 1]=-rot_sin*rmat.m[0][0]+rot_cos*rmat.m[1][0];
    rotation_matrix[ 2]= rmat.m[2][0];
    rotation_matrix[ 3]= 0.0;
    rotation_matrix[ 4]= rot_cos*rmat.m[0][1]+rot_sin*rmat.m[1][1];
    rotation_matrix[ 5]=-rot_sin*rmat.m[0][1]+rot_cos*rmat.m[1][1];
    rotation_matrix[ 6]= rmat.m[2][1];
    rotation_matrix[ 7]= 0.0;
    rotation_matrix[ 8]= rot_cos*rmat.m[0][2]+rot_sin*rmat.m[1][2];
    rotation_matrix[ 9]=-rot_sin*rmat.m[0][2]+rot_cos*rmat.m[1][2];
    rotation_matrix[10]= rmat.m[2][2];
    rotation_matrix[11]= 0.0;
    rotation_matrix[12]= 0.0;
    rotation_matrix[13]= 0.0;
    rotation_matrix[14]= 0.0;
    rotation_matrix[15]= 1.0;
    
    glLoadIdentity();
    glTranslatef((GLfloat)0.0,(GLfloat)0.0,(GLfloat)-30.0);
    glMultMatrixd(rotation_matrix);
  }
  else
  {
    glLoadIdentity();
    glTranslatef((GLfloat)0.0,(GLfloat)0.0,(GLfloat)-30.0);
  }
  
    if (spotlight==1)
    {
      glEnable(GL_LIGHT0);
      glPushMatrix();
        glRotatef((GLfloat)(0.5*spinxlight),(GLfloat)1.0,(GLfloat)0.0,(GLfloat)0.0);
        glRotatef((GLfloat)(0.5*spinylight),(GLfloat)0.0,(GLfloat)1.0,(GLfloat)0.0);
        glLightfv(GL_LIGHT0,GL_POSITION,position);
        glLightfv(GL_LIGHT0,GL_SPECULAR,spec);
      glPopMatrix();
    }
    else
    {
      glDisable(GL_LIGHT0);
    }
  
  draw_surface();
    
  glPopMatrix();

  glutSwapBuffers();
}

void special(int key, int x, int y)
{
  int i;
  
  switch(key)
  {
    case GLUT_KEY_LEFT:
      rot_angle--;
      rot_cos=cos(0.01*rot_angle);
      rot_sin=sin(0.01*rot_angle);
      glutPostRedisplay();
      break;
    case GLUT_KEY_RIGHT:
      rot_angle++;
      rot_cos=cos(0.01*rot_angle);
      rot_sin=sin(0.01*rot_angle);
      glutPostRedisplay();
      break;
    case GLUT_KEY_UP:
      choose_surface=0;
      flow=1;
      for (i=0; i<N; i++) { oh[i]=h[i]; om[i]=m[i]; osm[i]=sm[i]; }
      for (i=0; i<16; i++) { para_io[i]=para_i[i]; para_lo[i]=para_l[i]; }
      dt=fabs(dt);
      step();
      if (make_xy()==1)
      {
        still_ok=1;
        glutPostRedisplay();
        tm+=dt;
      }
      else
      {
        still_ok=0;
        for (i=0; i<N; i++) { h[i]=oh[i]; m[i]=om[i]; sm[i]=osm[i]; }
        for (i=0; i<16; i++) { para_i[i]=para_io[i]; para_l[i]=para_lo[i]; }
        rescale(); filter();
        step();
        if (make_xy()==1) { tm+=dt; still_ok=1; }
        else
        {
          for (i=0; i<N; i++) { h[i]=oh[i]; m[i]=om[i]; sm[i]=osm[i]; }
          for (i=0; i<16; i++) { para_i[i]=para_io[i]; para_l[i]=para_lo[i]; }
          make_xy();
          still_ok=0;
        }
        glutPostRedisplay();
      }
      break;
    case GLUT_KEY_DOWN:
      choose_surface=0;
      flow=1;
      for (i=0; i<N; i++) { oh[i]=h[i]; om[i]=m[i]; osm[i]=sm[i]; }
      for (i=0; i<16; i++) { para_io[i]=para_i[i]; para_lo[i]=para_l[i]; }
      dt=-fabs(dt);
      step();
      if (make_xy()==1)
      {
        still_ok=1;
        glutPostRedisplay();
      }
      else
      {
        still_ok=0;
        for (i=0; i<N; i++) { h[i]=oh[i]; m[i]=om[i]; sm[i]=osm[i]; }
        for (i=0; i<16; i++) { para_i[i]=para_io[i]; para_l[i]=para_lo[i]; }
        rescale(); filter();
        step();
        if (make_xy()==1) still_ok=1;
        else
        {
          for (i=0; i<N; i++) { h[i]=oh[i]; m[i]=om[i]; sm[i]=osm[i]; }
          for (i=0; i<16; i++) { para_i[i]=para_io[i]; para_l[i]=para_lo[i]; }
          make_xy();
          still_ok=0;
        }
        glutPostRedisplay();
      }
      break;
  }
}

void myReshape(int w, int h)
{
  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0,(GLfloat)w/(GLfloat)h,1.0,100.0);
  glMatrixMode(GL_MODELVIEW);
  glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'n':
      choose_surface=1;
      flow=0;
      movelight=0;
      rotate=0;
      init();
      make_xy();
      glutPostRedisplay();
      break;
    case 'f':
      choose_surface=0;
      flow=1;
      movelight=0;
      rotate=0;
      glutPostRedisplay();
      break;
    case 'm':
      metric=1;
      glutPostRedisplay();
      break;
    case 's':
      metric=0;
      glutPostRedisplay();
      break;
  }
}

int main(int argc, char** argv)
{
  int i, j;
  
  n_fft=50;
  still_ok=1;
  metric=0;
  do_deriv_adj=1;
  dt=0.0001;
  c5=c3=0.0;
  tm=0.0;
  init();
  make_xy();
  
  for (i=0; i<3; i++) for (j=0; j<3; j++) rmat.m[i][j]=refmat.m[i][j]=0.0;
  for (i=0; i<3; i++) rmat.m[i][i]=refmat.m[i][i]=1.0;
  
  glutInitWindowSize(800,800);
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_MULTISAMPLE);
  glutCreateWindow("Ricci_Rot");
  glEnable(GL_LIGHTING);       
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);    
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutSpecialFunc(special);
  glutKeyboardFunc(keyboard);
  glutReshapeFunc(myReshape);
  glutDisplayFunc(display);
  make_menus();
  glutMainLoop() ;
  return 0;
}
