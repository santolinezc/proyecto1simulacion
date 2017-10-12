//malla
#include <fstream>
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

using namespace std;

//------CONSTANTS------------------------------------------------------------
const int Nx=10,Ny=1,Nz=1,N=Nx*Ny;
const double L=10,K=10,K1=1e7,Kcundall=10;
const double Gamma=50,Mu=0.4;
const double gamma1=0.25;
const double g=1,G=1;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;

//-------DECLARATIONS---------------------------------------------------------
//classes--------------------
class mass;
class collider;

//functions------------------------
void movedt(mass *Masses,collider Hooke, double dt);
void initanimation2D(void);
void startdrawing2D(void);
void initanimation3D(void);
void startdrawing3D(void);

void enddrawing(void);

//--------CLASSES-----------------------------------------------------------

class mass{
private:
  double m, R,theta,I;
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
  void drawmass2D(void);
  void drawmass3D(void);
  friend class mesh;
  friend class collider;
};



class collider{
private:
  vector3D interp[N]; bool interacting[N];
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
  mass Springmass[N];
  collider Hooke;

  int i,j,k, Ndraw;
  vector3D v0,r0 ,w0,l;
  
  double m0=1,R0=1.;
  double t,tdraw, dt=1e-1, tmax=200;
  
  r0.cargue(0,0,0);
  v0.cargue(0,0,0);
  w0.cargue(0,0,0);

  initanimation2D();Ndraw=15000;
   
  //r0.cargue(Nx*L/2,L*Ny+10*R0,0);
  
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx;j++){
      //if(j==0) l.cargue((j)*L,j*L/2,0);
       l.cargue(j*L,i*L,0);
       Springmass[i*Nx+j].init(m0,0,r0+l,v0,w0,R0);
    }
  }
  r0.cargue(L/2,0,0);
  Springmass[0].init(m0,0,r0,v0,w0,R0);
  //r0.cargue((Nx-1)*L,3*L/2,0);
  //Springmass[Nx+Nx-1].init(m0,0,r0,v0,w0,R0);
  for(t=tdraw=0; t < tmax; t+=dt,tdraw+=dt){
     
    if(tdraw > tmax/Ndraw){
      startdrawing2D();
      for(i=0;i<N;i++)Springmass[i].drawmass2D();
      enddrawing();
      tdraw=0;
      }
    //cout<<t<<"\t ";for(i=0;i<N;i++)cout<<Springmass[i].Getx()<<" "<<Springmass[i].Gety()<</*" "<<Springmass[i].Getz()<<*/"\t";cout<<endl;
 
     
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
  cout << ", " << r.x() << "+" << R << "*cos(t)," << r.y() << "+" << R << "*sin(t)";
  //cout << r.x()<< "+"<<R*cos(theta)/7. <<"*t" <<","<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

//----------------Functions collider--------------

void collider::iniciate(void)
{
  int i,j;
  for(i=0;i<N;i++){
    interp[i].cargue(0,0,0);
    interacting[i]=false;
  }
}

void collider::allforces(mass *Masses,double dt)
{
  int i,j,k,q;
  vector3D x0,y0,x,y,L0,L1,F,ru,ru1;
 
  for(i=0;i<N;i++)Masses[i].eraseforce();

  /*/Wall boundry conditions
  for(i=0;i<Ny;i++){
    L0.cargue(-L,i*L,0); L1.cargue((Nx)*L,i*L,0);
    x=Masses[i*Nx].r-L0; y= Masses[(i+1)*Nx-1].r-L1;
    ru=x/norma(x);       ru1=y/norma(y);
    x0=L*ru;             y0=L*ru1;
    F=-K*(x-x0);         Masses[i*Nx].addforce(F);
    F=-K*(y-y0);         Masses[(i+1)*Nx-1].addforce(F);
    }*/

  // horizontal springs
  for(i=0;i<Ny;i++){
    for(j=0;j<Nx-1;j++){
      springforcebetween(Masses[i*Nx+j],Masses[i*Nx+j+1]);
    }
  }
  
  //Vertical springs
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
    for(i=0;i<Ny;i++)periodicconditions(Masses[(i+1)*Nx-1],Masses[i*Nx]);
  
  //for(i=0;i<N;i++)interactionforce(Masses[i],Masses[N],interp[i],interacting[i],dt);
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
  nxvec.cargue((Nx+1)*L,0,0);
  x=Mass2.r-Mass1.r+nxvec;
  ru=x/norma(x);
  x0=L*ru;
  F=-K*(x-x0);
  Mass2.addforce(F);
  Mass1.addforce((-1)*F);
}

void collider::interactionforce(mass & Mass1,mass & Mass2,vector3D &interp,bool & interacting, double dt)
{
  /* vector3D F,Ft,Fn,r21,n,vc,vcn,vct,t;
  double d,s,m1,m2,m21,R1,R2,normvct,compFn,compvcn,normFt,Ftmax;
  double ERRf=1e-8;
  
  r21=Mass2.r-Mass1.r;
  d=norma(r21);
  s=(Mass2.R-Mass1.R)-d;

  if(s>0){
    m1=Mass1.m; m2=Mass2.m; m21=(m1*m2)/(m1+m2);
    R1=Mass1.R; R2=Mass2.R;
    n=r21/d;

    vc=(Mass2.v-Mass1.v)-(Mass2.w^n)*R2-(Mass1.w^n)*R1;
    compvcn=vc*n;
    vcn=compvcn*n;
    vct=vc-vcn; normvct=norma(vct);

    if(normvct<ERRf)t.cargue(0,0,0);else t=vct/normvct;

    //Normal forces

    compFn=(K1*pow(s,1.5));
    compFn-=m21*sqrt(s)*Gamma*compvcn; if(compFn<0) compFn=0;

    Fn=n*compFn;

    //Tangent forces

    interp+=(vct*dt);
    Ft=interp*(-Kcundall);

    Ftmax=Mu*compFn; normFt=norma(Ft);
    if(normFt>Ftmax) Ft=interp*(-Ftmax/norma(interp));

    F=Fn+Ft;
    Mass2.addforce(F); Mass1.addtorque((n*(-R2)^Ft));
    Mass1.addforce(F*(-1)); Mass1.addtorque((n*(-R1)^((-1)*Ft)));
    interacting=true;
  }
  else if(interacting==true){
    interp.cargue(0,0,0); interacting=false;
    }*/
}

//----------General Functions---------------------------

void movedt(mass *Masses,collider Hooke, double dt)
{
  //move with Omelyan
  int i;
  for(i=0;i<N;i++) Masses[i].move_r(dt,Zeta);
   Hooke.allforces(Masses,dt);
  for(i=0;i<N;i++)Masses[i].move_v(dt,(1-2*Lambda)/2);
  
  for(i=0;i<N;i++)Masses[i].move_r(dt,Xi);
  Hooke.allforces(Masses,dt);
  for(i=0;i<N;i++) Masses[i].move_v(dt,Lambda);

  for(i=0;i<N;i++)Masses[i].move_r(dt,1-2*(Xi+Zeta));
  Hooke.allforces(Masses,dt);
  for(i=0;i<N;i++)Masses[i].move_v(dt,Lambda);

  for(i=0;i<N;i++)Masses[i].move_r(dt,Xi);
  Hooke.allforces(Masses,dt);
  for(i=0;i<N;i++) Masses[i].move_v(dt,(1-2*Lambda)/2);
  
  for(i=0;i<N;i++)Masses[i].move_r(dt,Zeta);
}

//--------------------ANIMACION------------------------------

void initanimation2D(void)
{
  cout << "set term gif animate" << endl;
  cout << "set output 'gify.gif'" << endl; 
  cout << "unset key" << endl;
  cout << "set xrange [" << -L << ":"<<(Nx)*L << "]" << endl;
  cout << "set yrange [" << -L << ":"<<(Ny)*L << "]" << endl;
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
  cout << "set yrange [" << -L << ":"<<(Ny+1)*L << "]" << endl;
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


