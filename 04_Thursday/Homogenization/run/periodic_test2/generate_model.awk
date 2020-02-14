BEGIN{
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#-to double the density of layer, set multi to 2
#-to compute the reference solution, set amp to 0
#
  multi=4.;
# multi=0.4;
  amp=0.4;
#  amp=0.0;
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  npart=200*multi;
  c0=1500.;p0=1000.; 
  E0=p0*c0^2;
  E0m1=1./E0;  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bufpart=00*multi;len=1.; 
  dx=len/npart;
  print npart;
  x1=dx;
  print 0., 2;
  print 0.,p0,c0;
  print x1,p0,c0;
    print 0.,c0>"velocity";
    print x1,c0>"velocity";
    print 0.,c0*c0*p0>"mu";
    print 0.,c0*c0*p0>"mu";
    print 0,p0>"density";
    print x1,p0>"density";
  for (i=2;i<npart;i++) {
    if (i>=bufpart && i <= npart-bufpart){
      x=rand();
      x2= i   *dx;
      if (i%2 == 0 ) { x=1. } else { x =0.} 
           Em1=E0m1*(1.+2*amp*(x-0.5));
           p2=p0*(1.-0.5*amp*(x-0.5));
           c2=sqrt(1./p2/Em1);
#
           Em1=E0m1*(1.+2*amp*(x-0.5));
           p1=p0*(1.-0.5*amp*(x-0.5));
           c1=sqrt(1./p1/Em1);
    }else{
            x2= i   *dx;
            c2=c0;
            p2=p0;
            c1=c0;
            p1=p0;
    }
    print x1, 2;
    print x1,p1,c1;
    print x2,p2,c2;
    print x1,c1>"velocity";
    print x2,c2>"velocity";
    print x1,c1*c1*p1>"mu";
    print x2,c1*c1*p1>"mu";
    print x1,p1>"density";
    print x2,p2>"density";
    x1=x2
  }
  print x1, 2;
  print x1,p0,c0;
  print len,p0,c0;
    print x1,c0>"velocity";
    print len,c0>"velocity";
    print x1,p0>"density";
    print len,p0>"density";
}
#function wt(x)
#{
#  pi=3.14159265358979323844;
#  xx1=00000;
#  xx2=00000;
#  if (x<xx1){
#    w=0.; }
#  else if (x>xx2){
#    w=1.;}
#  else{
#     w=1.-0.5*(1.0+cos(pi*(x-xx1)/(xx2-xx1)));
#  }
#  return w
#}
