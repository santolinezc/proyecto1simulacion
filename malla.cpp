//malla
#include <fstream>
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

using namespace std;

//------CONSTANTS------------------------------------------------------------
const int Nx=100,Ny=3,Nz=1,N=Nx*Ny;
const double L=3,K=20,K1=1e3,Kcundall=10;
const double Gamma=30,MU=0.4;
const double Gamma1=1.5;
const double gamma1=0.25;
const double g=1;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;

//-------DECLARATIONS---------------------------------------------------------
//classes--------------------
class mass;
class collider;

//functions------------------------
void movedt(mass *Masses,collider Hooke, double dt);
void initanimation2D(double h);
void startdrawing2D(void);
void initanimation3D(void);
void startdrawing3D(void);

void enddrawing(void);

//--------CLASSES-----------------------------------------------------------

class mass{
private:
  double m, R,theta,I, theta_r;
  vector3D F,r,v, w,tau;
public:
  void init(double m0,double theta0,vector3D r0,vector3D v0,vector3D w0,double R0);
  void eraseforce(void);
  void addforce(vector3D F0);
  void addtorque(vector3D tau0);
  void move_r(double dt,double cte);
  void move_v(double dt,double cte);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
  double Getvx(void){return v.x();};
  double Gettheta_r(void){return theta_r;};
  void drawmass2D(void);
  void drawmass3D(void);
  friend class mesh;
  friend class collider;
};



class collider{
private:
  vector3D interp[N+1][N+1]; bool interacting[N+1][N+1];
public:
  void iniciate(void);
  void allforces(mass *Masses,double dt);
  void interactionforce(mass & Mass1,mass & Mass2,vector3D & inerp,bool & interacting,double dt);
  void springforcebetween(mass & Mass1,mass & Mass2);
  void periodicconditions(mass & Mass1,mass & Mass2);
};

//--------MAIN-----------------------------------------------------------------

int main (void)
{
  mass Springmass[N+1];
  collider Hooke;
    
  int i,j,k, Ndraw;
  vector3D v0,r0 ,w0,l;
  
  double m0=2,R0=0.5;
  double t,tdraw, dt=1e-1, tmax=75;

  double mmass=m0/2., Rmass=20*R0, h=(Ny-1)*L+Rmass+R0;

  initanimation2D(h);Ndraw=15000;
  
  r0.cargue(0,0,0);
  v0.cargue(0,0,0);
  w0.cargue(0,0,0);
 
  //iniciates mesh
  Hooke.iniciate();
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx;j++){
       l.cargue(j*L,i*L,0);
       Springmass[i*Nx+j].init(m0,0,r0+l,v0,w0,R0);
    }
  }
  
  //iniciates mass
  
  r0.cargue((Nx-1)*L,h,0);
  w0.cargue(0,0,2);
  Springmass[N].init(mmass,0,r0,v0,w0,Rmass);

  for(t=tdraw=0; t < tmax; t+=dt,tdraw+=dt){
    
    if(tdraw > tmax/Ndraw){
      startdrawing2D();
      for(i=0;i<N+1;i++)Springmass[i].drawmass2D();
      enddrawing();
      tdraw=0;
      }
    //cout<<t<<"\t ";for(i=0;i<N;i++)cout<<Springmass[i].Getx()<<" "<<Springmass[i].Gety()<</*" "<<Springmass[i].Getz()<<*/"\t";cout<<endl;
    //cout<<Springmass[N].Getvx()<<" "<<Springmass[N].Gettheta_r()<<endl;
     
    movedt(Springmass,Hooke,dt);
    
   
      
    }
  
  return 0;
}


//--------FUNCTIONS-----------------------------------------

//---------Functions mass-------------

void mass::init(double m0,double theta0,vector3D r0,vector3D v0,vector3D w0,double R0)
{
  m=m0; v=v0;r=r0; R=R0;theta=theta0;w=w0;
  F.cargue(0,0,0); tau.cargue(0,0,0); I=2*m*R*R/5;
  theta_r=0;
}

void mass::eraseforce(void)
{
  F.cargue(0,0,0); tau.cargue(0,0,0);
}

void mass::addforce(vector3D F0)
{
  F+=F0;
}

void mass::addtorque(vector3D tau0)
{
  tau+=tau0;
}

void mass::move_r(double dt,double cte)
{
  vector3D Id;
  Id.cargue(0,0,1);
  r+=v*(dt*cte); theta+=Id*w*(cte*dt);
}

void mass::move_v(double dt,double cte)
{
  v+=F*(dt*cte/m); w+=tau*(cte*dt)/I;
}

void mass::drawmass3D(void)
{
  cout << ", " << r.x() << "+" << R << "*cos(u)*sin(v)," << r.y() << "+" << R << "*sin(u)*sin(v),"<< r.z() << "+" << R << "*cos(v) ";
  //cout << r.x()<< "+"<<R*cos(theta)/7. <<"*t" <<","<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

void mass::drawmass2D(void)
{
  cout << ", " << r.x() << "+" << R << "*cos(t)," << r.y() << "+" << R << "*sin(t),";
  cout << r.x()<< "+"<<R*cos(theta)/7. <<"*t" <<","<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

//----------------Functions collider--------------

void collider::iniciate(void)
{
  int i,j;
  for(i=0;i<N;i++){
    for(j=i+1;j<N+1;j++){
      interp[i][j].cargue(0,0,0);
      interacting[i][j]=false;
    }
  }
}

void collider::allforces(mass *Masses,double dt)
{
  int i,j,k,q;
  vector3D x0,y0,x,y,L0,L1,F,ru,ru1,gvec;
  //erase forces
  for(i=0;i<N+1;i++)Masses[i].eraseforce();

  //gravity
  gvec.cargue(0,-g,0);
  Masses[N].addforce(gvec*Masses[N].m);
  //Wall boundry conditions
  for(i=0;i<Ny;i++){
    L0.cargue(-L,i*L,0); L1.cargue((Nx)*L,i*L,0);
    x=Masses[i*Nx].r-L0; y= Masses[(i+1)*Nx-1].r-L1;
    ru=x/norma(x);       ru1=y/norma(y);
    x0=L*ru;             y0=L*ru1;
    F=-K*(x-x0);         Masses[i*Nx].addforce(F);
    F=-K*(y-y0);         Masses[(i+1)*Nx-1].addforce(F);
    }

  // horizontal springs in x
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx-1;j++){
      springforcebetween(Masses[i*Nx+j],Masses[i*Nx+j+1]);
    }
  }
  
  //horizontal springs in y
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx;j++){
      if(i<Ny-1) springforcebetween(Masses[i*Nx+j],Masses[(i+1)*Nx+j]);
    }
  }
  
  // Damping
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx;j++){
      Masses[i*Nx+j].addforce(-Masses[i*Nx+j].m*gamma1*Masses[i*Nx+j].v);
    }
  }
 
  
  //Periodic boundry conditions
  //for(i=0;i<Ny;i++)periodicconditions(Masses[(i+1)*Nx-1],Masses[i*Nx]);
  
  //interaction force between mass and mesh
  for(i=0;i<N;i++){interactionforce(Masses[i],Masses[N],interp[i][N],interacting[i][N],dt);}

}

void collider::springforcebetween(mass & Mass1,mass &Mass2)
{
  vector3D x0,x,F,F1, ru,gvec;  
  x=Mass2.r-Mass1.r;
  ru=x/norma(x);
  x0=L*ru;
  F=-K*(x-x0);
  Mass2.addforce(F);
  Mass1.addforce((-1)*F);
}

void collider::periodicconditions(mass & Mass1,mass & Mass2)
{
  vector3D x0,x,F,F1, ru,gvec,nxvec;
  nxvec.cargue((Nx)*L,0,0);
  x=Mass2.r-Mass1.r+nxvec;
  ru=x/norma(x);
  x0=L*ru;
  F=-K*(x-x0);
  Mass2.addforce(F);
  Mass1.addforce((-1)*F);
}

void collider::interactionforce(mass & Mass1,mass & Mass2,vector3D &interp,bool & interacting, double dt)
{
  int i,j;
  vector3D r21,n,vc,vcn,vct,Fn,Ft,F,t,F1,F2,tau,Id;
  double R1,R2,d,s,m1,m2,m21,compvcn,compFn,normvct,Ftmax,normFt;
  double ERFF=1e-8;
  //Id.cargue(0,-g,0);
  R1=Mass1.R; R2=Mass2.R;
  
  r21=Mass2.r-Mass1.r; d=norma(r21); n=r21/d;
  s=R1+R2-d;

  if(s>0){
    m1=Mass1.m; m2=Mass2.m; m21=m1*m2/(m1+m2);

    // contact velocity
    vc=Mass2.v-Mass1.v -(Mass2.w^n)*R2-(Mass1.w^n)*R1;
    compvcn=vc*n; vcn=compvcn*n; vct=vc-vcn;normvct=norma(vct);
    if(normvct<ERFF) t.cargue(0,0,0); else t=vct/normvct;

    //Hertz force
    compFn=K1*pow(s,1.5);

    //plastic disipation
    compFn-=m21*sqrt(s)*Gamma*compvcn;
    if(compFn<0) compFn=0;

    Fn=n*compFn;

    //tangential forces

    interp+=(vct*dt);
    Ft=interp*(-Kcundall);
    Ftmax=MU*compFn; normFt=norma(Ft);
    if(normFt>Ftmax) Ft=interp*(-Ftmax/norma(interp));
    
    
    
    F=Fn+Ft;
    F2=F-(Mass2.w^n)*(R2*Gamma1*m2);
    F1=(-1)*F-(Mass1.w^n)*(R1*Gamma1*m1);
    
    //Mass2.theta_r=(sqrt((F*F)*R2*R2*(R2*R2*(F2*F2)-tau.z()*tau.z())));
    Mass2.addforce(F); Mass2.addtorque((n*(-R2))^(F-Fn));
    Mass1.addforce((-1)*F); Mass1.addtorque((n*R1)^(F1+Fn));

    interacting=true;
  }
  else if(interacting==true){interp.cargue(0,0,0);interacting==false;}
}

//----------General Functions---------------------------

void movedt(mass *Masses,collider Hooke, double dt)
{
  //move with Omelyan
  int i,q=N+1;
   
  for(i=0;i<q;i++) Masses[i].move_r(dt,Zeta);
   Hooke.allforces(Masses,dt);
  for(i=0;i<q;i++)Masses[i].move_v(dt,(1-2*Lambda)/2);
  
  for(i=0;i<q;i++)Masses[i].move_r(dt,Xi);
  Hooke.allforces(Masses,dt);
  for(i=0;i<q;i++) Masses[i].move_v(dt,Lambda);

  for(i=0;i<q;i++)Masses[i].move_r(dt,1-2*(Xi+Zeta));
  Hooke.allforces(Masses,dt);
  for(i=0;i<q;i++)Masses[i].move_v(dt,Lambda);

  for(i=0;i<q;i++)Masses[i].move_r(dt,Xi);
  Hooke.allforces(Masses,dt);
  for(i=0;i<q;i++) Masses[i].move_v(dt,(1-2*Lambda)/2);
  
  for(i=0;i<q;i++)Masses[i].move_r(dt,Zeta);
}

//--------------------ANIMACION------------------------------

void initanimation2D(double h)
{
  cout << "set term gif animate" << endl;
  cout << "set output 'gify.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [" << -L << ":"<<(Nx)*L << "]" << endl;
  cout << "set yrange [" << -L-h << ":"<<(Ny)*L+3*h/2 << "]" << endl;
  cout << "set size ratio -1" << endl;
  cout << "set parametric" << endl;
  cout << "set trange[0:7]" << endl;
  cout << "set isosamples 12" << endl;
}

void initanimation3D(void)
{
  cout << "set term gif animate" << endl;
  cout << "set output 'gify.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [" << -L << ":"<<(Nx+1)*L << "]" << endl;
  cout << "set yrange [" << -3*L << ":"<<(Ny+1)*L << "]" << endl;
  cout << "set zrange [" << -L << ":"<<(Nz+1)*L << "]" << endl;
  cout << "set size ratio -1" << endl;
  cout << "set parametric" << endl;
  cout << "set urange[0:7]" << endl;
  cout << "set vrange[0:3.14]" << endl;
  cout << "set isosamples 30" << endl;
}

void startdrawing2D(void)
{
  cout << " plot 0,0 ";
}

void startdrawing3D(void)
{
  cout << " splot 0,0,0 ";
}


void enddrawing(void)
{
  cout << endl;
}


