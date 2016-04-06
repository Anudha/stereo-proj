#!/usr/bin/python
from __future__ import division
import numpy as np
from Tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from PIL import Image
import PngImagePlugin
import ttk
import sys
from tkFileDialog import *
import os
#import pdb

pi=np.pi

################"
def unique_rows(a):
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    return np.unique(b).view(a.dtype).reshape(-1, a.shape[1])
##################

#Image._initialized=2
###################################################################"
##### Fonction projection  sur l'abaque
####################################################################

def proj(x,y,z): 
  
    if z==1: 
        X=0
        Y=0
    elif z<-0.000001:
        X=250
        Y=250
    else: 
            
        X=x/(1+z)
        Y=y/(1+z)
    
    return np.array([X,Y],float) 
    
###################################################################"
##### Fonction rotation 
####################################################################

def rotation(phi1,phi,phi2):
   phi1=phi1*pi/180;
   phi=phi*pi/180;
   phi2=phi2*pi/180;
   R=np.array([[np.cos(phi1)*np.cos(phi2)-np.cos(phi)*np.sin(phi1)*np.sin(phi2),
            -np.cos(phi)*np.cos(phi2)*np.sin(phi1)-np.cos(phi1)*
            np.sin(phi2),np.sin(phi)*np.sin(phi1)],[np.cos(phi2)*np.sin(phi1)
            +np.cos(phi)*np.cos(phi1)*np.sin(phi2),np.cos(phi)*np.cos(phi1)
            *np.cos(phi2)-np.sin(phi1)*np.sin(phi2), -np.cos(phi1)*np.sin(phi)],
            [np.sin(phi)*np.sin(phi2), np.cos(phi2)*np.sin(phi), np.cos(phi)]],float)
   return R

####################################################################
##### Fonction rotation autour d'un axe 
####################################################################

def Rot(th,a,b,c):
   th=th*pi/180;
   aa=a/np.linalg.norm([a,b,c]);
   bb=b/np.linalg.norm([a,b,c]);
   cc=c/np.linalg.norm([a,b,c]);
   c1=np.array([[1,0,0],[0,1,0],[0,0,1]],float)
   c2=np.array([[aa**2,aa*bb,aa*cc],[bb*aa,bb**2,bb*cc],[cc*aa,
                cc*bb,cc**2]],float)
   c3=np.array([[0,-cc,bb],[cc,0,-aa],[-bb,aa,0]],float)
   R=np.cos(th)*c1+(1-np.cos(th))*c2+np.sin(th)*c3

   return R    


####################################################################
##### Fonction cristal
####################################################################
def crist():
    global axes,axesh,D,Dstar,V
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    e=eval(e_entry.get())
    d2=eval(d_label_var.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    V=a*b*c*np.sqrt(1-(np.cos(alp)**2)-(np.cos(bet))**2-(np.cos(gam))**2+2*np.cos(alp)*np.cos(bet)*np.cos(gam))
    D=np.array([[a,b*np.cos(gam),c*np.cos(bet)],[0,b*np.sin(gam),  c*(np.cos(alp)-np.cos(bet)*np.cos(gam))/np.sin(gam)],[0,0,V/(a*b*np.sin(gam))]])
    Dstar=np.transpose(np.linalg.inv(D))
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
    axes=np.zeros(((2*e+1)**3-1,3))
    axesh=np.zeros(((2*e+1)**3-1,5))
    axesh[:,4]=col_trace.get()
    id=0
    for i in range(-e,e+1):
        for j in range(-e,e+1):
            for k in range(-e,e+1):
                if (i,j,k)!=(0,0,0):
                    d=1/(np.sqrt(np.dot(np.array([i,j,k]),np.dot(np.linalg.inv(G),np.array([i,j,k])))))
                    if d>d2*0.1*np.amax([a,b,c]):
                        if var_uvw.get()==0:                    
                            Ma=np.dot(Dstar,np.array([i,j,k],float))
                            axesh[id,0]=Ma[0]
                            axesh[id,1]=Ma[1]
                            axesh[id,2]=Ma[2]
                            axesh[id,3]=0
                            axes[id,:]=np.array([i,j,k],float)
                        else:
                            Ma=np.dot(D,np.array([i,j,k],float))
                            axesh[id,0]=Ma[0]
                            axesh[id,1]=Ma[1]
                            axesh[id,2]=Ma[2]
                            axesh[id,3]=1
                            axes[id,:]=np.array([i,j,k],float)
                        id=id+1
                    
    return axes,axesh,D,Dstar,V

def dm():
    global dmip
    a=f.add_subplot(111)    
    a.figure.clear()
    a=f.add_subplot(111)
    dmip=dmip-eval(d_entry.get())
    d_label_var.set(dmip)
    crist()
    trace()
    
    return dmip
def dp():
    global dmip
    a=f.add_subplot(111)    
    a.figure.clear()
    a=f.add_subplot(111)
    dmip=dmip+eval(d_entry.get())
    d_label_var.set(dmip)
    crist()    
    trace()
    
    return dmip 
    
####################################################################
##### Fonction angle entre deux directions
####################################################################

def angle():
   global D, Dstar
   c10=eval(c10_entry.get())
   c11=eval(c11_entry.get())
   c12=eval(c12_entry.get())
   c20=eval(c20_entry.get())
   c21=eval(c21_entry.get())
   c22=eval(c22_entry.get())
   c1=np.array([c10,c11,c12])
   c2=np.array([c20,c21,c22])
   if var_uvw.get()==0: 
       c1c=np.dot(Dstar,c1)
       c2c=np.dot(Dstar,c2)
   else:
       c1c=np.dot(Dstar,c1)
       c2c=np.dot(Dstar,c2)
   the=np.arccos(np.dot(c1c,c2c)/(np.linalg.norm(c1c)*np.linalg.norm(c2c)))                   
   thes=str(np.around(the*180/pi,decimals=2))        
   angle_var.set(thes)
               
####################################################################
##### Fonction largeur plan
####################################################################

def largeur_plan():
    global D, Dstar, M   
    plan1=eval(plan1_entry.get())   
    plan2=eval(plan2_entry.get())
    plan3=eval(plan3_entry.get())
    
    n=np.array([plan1,plan2,plan3])
    if var_uvw.get()==0: 
       n2=np.dot(Dstar,n)
       
    else:
        n2=np.dot(D,n)
    nr=np.dot(M,n2)
    la=np.zeros((1,41))
    k=0
    N=np.array([0,0,1])
    for t in range(-40,41,2):
        nri=np.dot(Rot(t,0,1,0),nr)
        angle=np.arccos(np.dot(nri,N)/np.linalg.norm(nri))          
        la[0,k]=np.cos(angle)
        k=k+1
    
    plt.plot(range(-40,41,2),la[0,:])
    plt.xlabel('Tilt angle')
    plt.ylabel('Apparent width')
    plt.grid(True)
    plt.show()

####################################################################
##### Fonction facteur de Schmid
####################################################################
def schmid():
    global D, Dstar,M
    b1=eval(b1_entry.get())   
    b2=eval(b2_entry.get())
    b3=eval(b3_entry.get())
    n1=eval(n1_entry.get())   
    n2=eval(n2_entry.get())
    n3=eval(n3_entry.get())
    n=np.array([n1,n2,n3])
    b=np.array([b1,b2,b3])
    if var_uvw.get()==0: 
       npr=np.dot(Dstar,n)
       bpr=np.dot(Dstar,b)
       
    else:
       npr=np.dot(Dstar,n)
       bpr=np.dot(Dstar,b)
    npr2=np.dot(M,npr)
    bpr2=np.dot(M,bpr)
    T=np.array([0,1,0])
    anglen=np.arccos(np.dot(npr2,T)/np.linalg.norm(npr2))
    angleb=np.arccos(np.dot(bpr2,T)/np.linalg.norm(bpr2))
    s=np.cos(anglen)*np.cos(angleb)
    s2=str(np.around(s,decimals=2))
    schmid_var.set(s2)

####################################################################
##### Fonction facteur de Schmid
####################################################################
def schmid_trace():
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid
    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
        
    tr_schmid=np.vstack((tr_schmid,np.array([pole1,pole2,pole3])))
    trace()
    
def undo_schmid_trace():
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    tr_s=tr_schmid
    for i in range(1,tr_schmid.shape[0]):
		if tr_schmid[i,0]==pole1 and tr_schmid[i,1]==pole2 and tr_schmid[i,2]==pole3:
			tr_s=np.delete(tr_schmid,i,0)
    tr_schmid=tr_s
    trace()

def fact(angle,r,t,n):
	x=r*np.cos(t)/n
	y=r*np.sin(t)/n
	f=np.cos(angle)*2*y/((1+x**2+y**2))
	return f    

def schmid_trace2(C):
    global D, Dstar,M,a
    for h in range(1,tr_schmid.shape[0]):
        b1=C[h,0]   
        b2=C[h,1]   
        b3=C[h,2]   
        b=np.array([b1,b2,b3])
        
        if var_uvw.get()==0: 
            bpr=np.dot(Dstar,b)/np.linalg.norm(np.dot(Dstar,b))
        else:
            bpr=np.dot(Dstar,b)/np.linalg.norm(np.dot(Dstar,b))
		  
        bpr2=np.dot(M,bpr)
        T=np.array([0,1,0])
        angleb=np.arccos(np.dot(bpr2,T)/np.linalg.norm(bpr2))
        n=300
        r=np.linspace(0,n,100)
        t=np.linspace(0,2*pi,100)
        r,t=np.meshgrid(r,t)
        F=fact(angleb,r,t,n)
        lev=[-0.5,-0.4,-0.3,-0.2,0.2,0.3,0.4,0.5]
        CS=a.contour(r*np.cos(t)+300, r*np.sin(t)+300, F,lev,linewidths=2)
        fmt = {}
        strs = [ '('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') -0.5','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') -0.4','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') -0.3','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') -0.2','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') 0.2','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') 0.3','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') 0.4','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') 0.5']
        for l,s in zip( CS.levels, strs ):
			fmt[l] = s
	a.clabel(CS,fmt=fmt,fontsize=10,inline=True)						

    
####################################################################
##### Fonction rotation projection
####################################################################
# on definit maintenant les rotations positives et negatives selon x, y et z
def rxp():
    global x,M,a,trP
    a = f.add_subplot(111)     
    a.figure.clear()
    thx=eval(rx_entry.get())
  
    M=np.dot(Rot(thx,1,0,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_var.set(t)
    x=x+eval(rx_entry.get())
    cx.set(x)
    return x,M
    
    
def rxm():
    global x,M,a
    a = f.add_subplot(111)     
    a.figure.clear()            
    thx=-eval(rx_entry.get())
    
    M=np.dot(Rot(thx,1,0,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_var.set(t)
    x=x-eval(rx_entry.get())
    cx.set(x)
    return x,M

def ryp():
    global y,M,a
    a = f.add_subplot(111)     
    a.figure.clear()
    thy=eval(ry_entry.get())
    
    M=np.dot(Rot(thy,0,1,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_var.set(t)
    y=y+eval(ry_entry.get())
    cy.set(y)
    return y,M
    
def rym():
    global y,M,a
    a = f.add_subplot(111)     
    a.figure.clear()
    thy=-eval(ry_entry.get())
    
    M=np.dot(Rot(thy,0,1,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_var.set(t)
    y=y-eval(ry_entry.get())
    cy.set(y)
    return y,M


def rzp():
    global z,M,a
    a = f.add_subplot(111)     
    a.figure.clear()
    thz=eval(rz_entry.get())
    
    M=np.dot(Rot(thz,0,0,1),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_var.set(t)
    z=z+eval(rz_entry.get())
    cz.set(z)
    return z,M
    
def rzm():
    global z,M,a
    a = f.add_subplot(111)     
    a.figure.clear()
    thz=-eval(rz_entry.get())
    
    M=np.dot(Rot(thz,0,0,1),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    angle_euler_var.set(t)
    z=z-eval(rz_entry.get())
    cz.set(z)
    return z,M

####################################################################
##### Fonction ajouter un pole
####################################################################
def pole(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar
    
#    fp=f.add_subplot(111)
#    
#       
#    Pp=np.zeros((1,2),float)
    if var_hexa.get()==1:
        if var_uvw.get()==1:
            pole1a=2*pole1+pole2
            pole2a=2*pole2+pole1
            pole1=pole1a
            pole2=pole2a
    
    Gs=np.array([pole1,pole2,pole3],float)
    #print(Gs)
    if var_uvw.get()==0:                    
            Gsh=np.dot(Dstar,Gs)/np.linalg.norm(np.dot(Dstar,Gs))
    else:
        Gsh=np.dot(D,Gs)/np.linalg.norm(np.dot(D,Gs))
     
    S=np.dot(M,Gsh)
    if S[2]<0:
        S=-S
        Gsh=-Gsh
        pole1=-pole1
        pole2=-pole2
        pole3=-pole3
#    Pp=proj(S[0],S[1],S[2])*600/2
#    l=str(int(pole1))+str(int(pole2))+str(int(pole3))
#    fp.plot(Pp[0]+600/2,Pp[1]+600/2,'go')    
#    fp.annotate(l,(Pp[0]+600/2,Pp[1]+600/2))
#    fp.axis([0,600,0,600])
#    fp.axis('off')   
#    fp.figure.canvas.draw() 
    
    axes=np.vstack((axes,np.array([pole1,pole2,pole3])))
    axes=np.vstack((axes,np.array([-pole1,-pole2,-pole3])))
    T=np.vstack((T,np.array([S[0],S[1],S[2]])))
    T=np.vstack((T,np.array([-S[0],-S[1],-S[2]])))
    if var_uvw.get()==0 :
        axesh=np.vstack((axesh,np.array([Gsh[0],Gsh[1],Gsh[2],0,col_trace.get()])))
        axesh=np.vstack((axesh,np.array([-Gsh[0],-Gsh[1],-Gsh[2],0,col_trace.get()])))
    else:
        axesh=np.vstack((axesh,np.array([Gsh[0],Gsh[1],Gsh[2],1,col_trace.get()])))
        axesh=np.vstack((axesh,np.array([-Gsh[0],-Gsh[1],-Gsh[2],1,col_trace.get()])))
    return axes,axesh,T

def undo_pole(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar
    
#    fp=f.add_subplot(111)
#    
#       
#    Pp=np.zeros((1,2),float)
    if var_hexa.get()==1:
        if var_uvw.get()==1:
            pole1a=2*pole1+pole2
            pole2a=2*pole2+pole1
            pole1=pole1a
            pole2=pole2a
    
    Gs=np.array([pole1,pole2,pole3],float)
    #print(Gs)
    if var_uvw.get()==0:                    
            Gsh=np.dot(Dstar,Gs)/np.linalg.norm(np.dot(Dstar,Gs))
    else:
        Gsh=np.dot(D,Gs)/np.linalg.norm(np.dot(D,Gs))
     
    S=np.dot(M,Gsh)
    if S[2]<0:
        S=-S
        Gsh=-Gsh
        pole1=-pole1
        pole2=-pole2
        pole3=-pole3
#    Pp=proj(S[0],S[1],S[2])*600/2
#    l=str(int(pole1))+str(int(pole2))+str(int(pole3))
#    fp.plot(Pp[0]+600/2,Pp[1]+600/2,'go')    
#    fp.annotate(l,(Pp[0]+600/2,Pp[1]+600/2))
#    fp.axis([0,600,0,600])
#    fp.axis('off')   
#    fp.figure.canvas.draw() 
    
    ind=np.where((axes[:,0]==pole1) & (axes[:,1]==pole2)& (axes[:,2]==pole3))
    indm=np.where((axes[:,0]==-pole1) & (axes[:,1]==-pole2)& (axes[:,2]==-pole3))
    axes=np.delete(axes,ind,0)
    axes=np.delete(axes,indm,0)
    T=np.delete(T,ind,0)
    T=np.delete(T,indm,0)
    if var_uvw.get()==0 :
        axesh=np.delete(axesh,ind,0)
        axesh=np.delete(axesh,indm,0)
    else:
        axesh=np.delete(axesh,ind,0)
        axesh=np.delete(axesh,indm,0)
    return axes,axesh,T

    
def d(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
    ds=(np.sqrt(np.dot(np.array([pole1,pole2,pole3]),np.dot(np.linalg.inv(G),np.array([pole1,pole2,pole3])))))
    return ds
def addpole_sym():
    global M,axes,axesh,T,V,D,Dstar,G    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])     
    v=d(pole1,pole2,pole3)
    
    pole(pole1,pole2,pole3)
    if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
        pole(pole1,pole2,pole3)
        pole(pole1,pole2,-pole3)
        pole(pole2,pole1,pole3)
        pole(pole2,pole1,-pole3)
        pole(-pole1-pole2,pole2,pole3)
        pole(-pole1-pole2,pole2,-pole3)
        pole(pole1,-pole1-pole2,pole3)
        pole(pole1,-pole1-pole2,-pole3)
        pole(pole2,-pole1-pole2,pole3)
        pole(pole2,-pole1-pole2,-pole3)
        pole(-pole1-pole2,pole1,pole3)
        pole(-pole1-pole2,pole1,-pole3)

    else:
        if np.abs(d(pole1,pole2,-pole3)-v)<0.001:
                pole(pole1,pole2,-pole3)
        if np.abs(d(pole1,-pole2,pole3)-v)<0.001:
                pole(pole1,-pole2,pole3)
        if np.abs(d(-pole1,pole2,pole3)-v)<0.001:
            pole(-pole1,pole2,pole3)
        if np.abs(d(pole2,pole1,pole3)-v)<0.001:
            pole(pole2,pole1,pole3)
        if np.abs(d(pole2,pole1,-pole3)-v)<0.001:
            pole(pole2,pole1,-pole3)
        if np.abs(d(pole2,-pole1,pole3)-v)<0.001:
            pole(pole2,-pole1,pole3)
        if np.abs(d(-pole2,pole1,pole3)-v)<0.001:
            pole(-pole2,pole1,pole3)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            pole(pole2,pole3,pole1)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            pole(pole2,pole3,-pole1)
        if np.abs(d(pole2,-pole3,pole1)-v)<0.001:
            pole(pole2,-pole3,pole1)
        if np.abs(d(-pole2,pole3,pole1)-v)<0.001:
            pole(-pole2,pole3,pole1)
        if np.abs(d(pole1,pole3,pole2)-v)<0.001:
            pole(pole1,pole3,pole2)
        if np.abs(d(pole1,pole3,-pole2)-v)<0.001:
            pole(pole1,pole3,-pole2)
        if np.abs(d(pole1,-pole3,pole2)-v)<0.001:
            pole(pole1,-pole3,pole2)
        if np.abs(d(-pole1,pole3,pole2)-v)<0.001:
            pole(-pole1,pole3,pole2)
        if np.abs(d(pole3,pole1,pole2)-v)<0.001:
            pole(pole3,pole1,pole2)
        if np.abs(d(pole3,pole1,-pole2)-v)<0.001:
            pole(pole3,pole1,-pole2)
        if np.abs(d(pole3,-pole1,pole2)-v)<0.001:
            pole(pole3,-pole1,pole2)
        if np.abs(d(-pole3,pole1,pole2)-v)<0.001:
            pole(-pole3,pole1,pole2)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            pole(pole3,pole2,pole1)
        if np.abs(d(pole3,pole2,-pole1)-v)<0.001:
            pole(pole3,pole2,-pole1)
        if np.abs(d(pole3,-pole2,pole1)-v)<0.001:
            pole(pole3,-pole2,pole1)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            pole(pole3,pole2,pole1)
    trace()

def undo_sym():
    global M,axes,axesh,T,V,D,Dstar,G    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])     
    v=d(pole1,pole2,pole3)
    
    undo_pole(pole1,pole2,pole3)
    if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
        undo_pole(pole1,pole2,pole3)
        undo_pole(pole1,pole2,-pole3)
        undo_pole(pole2,pole1,pole3)
        undo_pole(pole2,pole1,-pole3)
        undo_pole(-pole1-pole2,pole2,pole3)
        undo_pole(-pole1-pole2,pole2,-pole3)
        undo_pole(pole1,-pole1-pole2,pole3)
        undo_pole(pole1,-pole1-pole2,-pole3)
        undo_pole(pole2,-pole1-pole2,pole3)
        undo_pole(pole2,-pole1-pole2,-pole3)
        undo_pole(-pole1-pole2,pole1,pole3)
        undo_pole(-pole1-pole2,pole1,-pole3)

    else:
        if np.abs(d(pole1,pole2,-pole3)-v)<0.001:
                undo_pole(pole1,pole2,-pole3)
        if np.abs(d(pole1,-pole2,pole3)-v)<0.001:
                undo_pole(pole1,-pole2,pole3)
        if np.abs(d(-pole1,pole2,pole3)-v)<0.001:
            undo_pole(-pole1,pole2,pole3)
        if np.abs(d(pole2,pole1,pole3)-v)<0.001:
            undo_pole(pole2,pole1,pole3)
        if np.abs(d(pole2,pole1,-pole3)-v)<0.001:
            undo_pole(pole2,pole1,-pole3)
        if np.abs(d(pole2,-pole1,pole3)-v)<0.001:
            undo_pole(pole2,-pole1,pole3)
        if np.abs(d(-pole2,pole1,pole3)-v)<0.001:
            undo_pole(-pole2,pole1,pole3)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            undo_pole(pole2,pole3,pole1)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            undo_pole(pole2,pole3,-pole1)
        if np.abs(d(pole2,-pole3,pole1)-v)<0.001:
            undo_pole(pole2,-pole3,pole1)
        if np.abs(d(-pole2,pole3,pole1)-v)<0.001:
            undo_pole(-pole2,pole3,pole1)
        if np.abs(d(pole1,pole3,pole2)-v)<0.001:
            undo_pole(pole1,pole3,pole2)
        if np.abs(d(pole1,pole3,-pole2)-v)<0.001:
            undo_pole(pole1,pole3,-pole2)
        if np.abs(d(pole1,-pole3,pole2)-v)<0.001:
            undo_pole(pole1,-pole3,pole2)
        if np.abs(d(-pole1,pole3,pole2)-v)<0.001:
            undo_pole(-pole1,pole3,pole2)
        if np.abs(d(pole3,pole1,pole2)-v)<0.001:
            undo_pole(pole3,pole1,pole2)
        if np.abs(d(pole3,pole1,-pole2)-v)<0.001:
            undo_pole(pole3,pole1,-pole2)
        if np.abs(d(pole3,-pole1,pole2)-v)<0.001:
            undo_pole(pole3,-pole1,pole2)
        if np.abs(d(-pole3,pole1,pole2)-v)<0.001:
            undo_pole(-pole3,pole1,pole2)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            undo_pole(pole3,pole2,pole1)
        if np.abs(d(pole3,pole2,-pole1)-v)<0.001:
            undo_pole(pole3,pole2,-pole1)
        if np.abs(d(pole3,-pole2,pole1)-v)<0.001:
            undo_pole(pole3,-pole2,pole1)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            undo_pole(pole3,pole2,pole1)
    trace()
    
    
def addpole():
    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    pole(pole1,pole2,pole3)
    trace()
    
def undo_addpole():
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid,nn
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    undo_pole(pole1,pole2,pole3)
    
    trace()
####################################################################
##### Fonction tracer plan
####################################################################
def trace_plan(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar,trP
    
    
    pole_i=0
    pole_c=col_trace.get()
    
    if var_hexa.get()==1:
        if var_uvw.get()==1:
            pole1=2*eval(pole1_entry.get())+eval(pole2_entry.get())
            pole2=2*eval(pole2_entry.get())+eval(pole1_entry.get())
            pole3=eval(pole3_entry.get())
            pole_i=1
    
        
    trP=np.vstack((trP,np.array([pole1,pole2,pole3,pole_i,pole_c])))
    b=np.ascontiguousarray(trP).view(np.dtype((np.void, trP.dtype.itemsize * trP.shape[1])))
    
    trP=np.unique(b).view(trP.dtype).reshape(-1, trP.shape[1])
    #print(trP)
    
    
def trace_addplan():
    global M,axes,axesh,T,V,D,Dstar,trP
    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    
    trace_plan(pole1,pole2,pole3)
    trace()
    
def undo_trace_addplan():
    global M,axes,axesh,T,V,D,Dstar,trP
    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    
    undo_trace_plan(pole1,pole2,pole3)
    trace()
    
def undo_trace_plan(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid
    
    ind=np.where((trP[:,0]==pole1) & (trP[:,1]==pole2)& (trP[:,2]==pole3))
    indm=np.where((trP[:,0]==-pole1) & (trP[:,1]==-pole2)& (trP[:,2]==-pole3))
    
    trP=np.delete(trP,ind,0)
    trP=np.delete(trP,indm,0)
    b=np.ascontiguousarray(trP).view(np.dtype((np.void, trP.dtype.itemsize * trP.shape[1])))
    
    trP=np.unique(b).view(trP.dtype).reshape(-1, trP.shape[1])
    

def trace_plan_sym():
    global M,axes,axesh,T,V,D,Dstar,G    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])     
    v=d(pole1,pole2,pole3)
    
    trace_plan(pole1,pole2,pole3)
    if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
        trace_plan(pole1,pole2,pole3)
        trace_plan(pole1,pole2,-pole3)
        trace_plan(pole2,pole1,pole3)
        trace_plan(pole2,pole1,-pole3)
        trace_plan(-pole1-pole2,pole2,pole3)
        trace_plan(-pole1-pole2,pole2,-pole3)
        trace_plan(pole1,-pole1-pole2,pole3)
        trace_plan(pole1,-pole1-pole2,-pole3)
        trace_plan(pole2,-pole1-pole2,pole3)
        trace_plan(pole2,-pole1-pole2,-pole3)
        trace_plan(-pole1-pole2,pole1,pole3)
        trace_plan(-pole1-pole2,pole1,-pole3)

    else:
        if np.abs(d(pole1,pole2,-pole3)-v)<0.001:
                trace_plan(pole1,pole2,-pole3)
        if np.abs(d(pole1,-pole2,pole3)-v)<0.001:
                trace_plan(pole1,-pole2,pole3)
        if np.abs(d(-pole1,pole2,pole3)-v)<0.001:
            trace_plan(-pole1,pole2,pole3)
        if np.abs(d(pole2,pole1,pole3)-v)<0.001:
            trace_plan(pole2,pole1,pole3)
        if np.abs(d(pole2,pole1,-pole3)-v)<0.001:
            trace_plan(pole2,pole1,-pole3)
        if np.abs(d(pole2,-pole1,pole3)-v)<0.001:
            trace_plan(pole2,-pole1,pole3)
        if np.abs(d(-pole2,pole1,pole3)-v)<0.001:
            trace_plan(-pole2,pole1,pole3)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            trace_plan(pole2,pole3,pole1)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            trace_plan(pole2,pole3,-pole1)
        if np.abs(d(pole2,-pole3,pole1)-v)<0.001:
            trace_plan(pole2,-pole3,pole1)
        if np.abs(d(-pole2,pole3,pole1)-v)<0.001:
            trace_plan(-pole2,pole3,pole1)
        if np.abs(d(pole1,pole3,pole2)-v)<0.001:
            trace_plan(pole1,pole3,pole2)
        if np.abs(d(pole1,pole3,-pole2)-v)<0.001:
            trace_plan(pole1,pole3,-pole2)
        if np.abs(d(pole1,-pole3,pole2)-v)<0.001:
            trace_plan(pole1,-pole3,pole2)
        if np.abs(d(-pole1,pole3,pole2)-v)<0.001:
            trace_plan(-pole1,pole3,pole2)
        if np.abs(d(pole3,pole1,pole2)-v)<0.001:
            trace_plan(pole3,pole1,pole2)
        if np.abs(d(pole3,pole1,-pole2)-v)<0.001:
            trace_plan(pole3,pole1,-pole2)
        if np.abs(d(pole3,-pole1,pole2)-v)<0.001:
            trace_plan(pole3,-pole1,pole2)
        if np.abs(d(-pole3,pole1,pole2)-v)<0.001:
            trace_plan(-pole3,pole1,pole2)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            trace_plan(pole3,pole2,pole1)
        if np.abs(d(pole3,pole2,-pole1)-v)<0.001:
            trace_plan(pole3,pole2,-pole1)
        if np.abs(d(pole3,-pole2,pole1)-v)<0.001:
            trace_plan(pole3,-pole2,pole1)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            trace_plan(pole3,pole2,pole1)
    trace()
def undo_trace_plan_sym():
    global M,axes,axesh,T,V,D,Dstar,G    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])     
    v=d(pole1,pole2,pole3)
    
    undo_trace_plan(pole1,pole2,pole3)
    if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
        undo_trace_plan(pole1,pole2,pole3)
        undo_trace_plan(pole1,pole2,-pole3)
        undo_trace_plan(pole2,pole1,pole3)
        undo_trace_plan(pole2,pole1,-pole3)
        undo_trace_plan(-pole1-pole2,pole2,pole3)
        undo_trace_plan(-pole1-pole2,pole2,-pole3)
        undo_trace_plan(pole1,-pole1-pole2,pole3)
        undo_trace_plan(pole1,-pole1-pole2,-pole3)
        undo_trace_plan(pole2,-pole1-pole2,pole3)
        undo_trace_plan(pole2,-pole1-pole2,-pole3)
        undo_trace_plan(-pole1-pole2,pole1,pole3)
        undo_trace_plan(-pole1-pole2,pole1,-pole3)

    else:
        if np.abs(d(pole1,pole2,-pole3)-v)<0.001:
                undo_trace_plan(pole1,pole2,-pole3)
        if np.abs(d(pole1,-pole2,pole3)-v)<0.001:
                undo_trace_plan(pole1,-pole2,pole3)
        if np.abs(d(-pole1,pole2,pole3)-v)<0.001:
            undo_trace_plan(-pole1,pole2,pole3)
        if np.abs(d(pole2,pole1,pole3)-v)<0.001:
            undo_trace_plan(pole2,pole1,pole3)
        if np.abs(d(pole2,pole1,-pole3)-v)<0.001:
            undo_trace_plan(pole2,pole1,-pole3)
        if np.abs(d(pole2,-pole1,pole3)-v)<0.001:
            undo_trace_plan(pole2,-pole1,pole3)
        if np.abs(d(-pole2,pole1,pole3)-v)<0.001:
            undo_trace_plan(-pole2,pole1,pole3)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            undo_trace_plan(pole2,pole3,pole1)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            undo_trace_plan(pole2,pole3,-pole1)
        if np.abs(d(pole2,-pole3,pole1)-v)<0.001:
            undo_trace_plan(pole2,-pole3,pole1)
        if np.abs(d(-pole2,pole3,pole1)-v)<0.001:
            undo_trace_plan(-pole2,pole3,pole1)
        if np.abs(d(pole1,pole3,pole2)-v)<0.001:
            undo_trace_plan(pole1,pole3,pole2)
        if np.abs(d(pole1,pole3,-pole2)-v)<0.001:
            undo_trace_plan(pole1,pole3,-pole2)
        if np.abs(d(pole1,-pole3,pole2)-v)<0.001:
            undo_trace_plan(pole1,-pole3,pole2)
        if np.abs(d(-pole1,pole3,pole2)-v)<0.001:
            undo_trace_plan(-pole1,pole3,pole2)
        if np.abs(d(pole3,pole1,pole2)-v)<0.001:
            undo_trace_plan(pole3,pole1,pole2)
        if np.abs(d(pole3,pole1,-pole2)-v)<0.001:
            undo_trace_plan(pole3,pole1,-pole2)
        if np.abs(d(pole3,-pole1,pole2)-v)<0.001:
            undo_trace_plan(pole3,-pole1,pole2)
        if np.abs(d(-pole3,pole1,pole2)-v)<0.001:
            undo_trace_plan(-pole3,pole1,pole2)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            undo_trace_plan(pole3,pole2,pole1)
        if np.abs(d(pole3,pole2,-pole1)-v)<0.001:
            undo_trace_plan(pole3,pole2,-pole1)
        if np.abs(d(pole3,-pole2,pole1)-v)<0.001:
            undo_trace_plan(pole3,-pole2,pole1)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            undo_trace_plan(pole3,pole2,pole1)
    trace()    
def trace_plan2(B):
    global M,axes,axesh,T,V,D,Dstar,a
    
   
    
    for h in range(0,B.shape[0]):
        pole1=B[h,0]
        pole2=B[h,1]
        pole3=B[h,2]
        Gs=np.array([pole1,pole2,pole3],float)
        if B[h,3]==0:                    
            Gsh=np.dot(Dstar,Gs)/np.linalg.norm(np.dot(Dstar,Gs))
        else:
            Gsh=np.dot(D,Gs)/np.linalg.norm(np.dot(D,Gs))
        S=np.dot(M,Gsh)
        
        if S[2]<0:
            S=-S
            Gsh=-Gsh
            pole1=-pole1
            pole2=-pole2
            pole3=-pole3
        r=np.sqrt(S[0]**2+S[1]**2+S[2]**2)
        A=np.zeros((2,100))
        Q=np.zeros((1,2))
        if S[2]==0:
             t=90
             
        else:
             t=np.arctan2(S[1],S[0])*180/pi
        w=0
        ph=np.arccos(S[2]/r)*180/pi
        for g in np.linspace(-pi,pi-0.00001,100):
            Aa=np.dot(Rot(t,0,0,1),np.dot(Rot(ph,0,1,0),np.array([np.sin(g),np.cos(g),0])))
            A[:,w]=proj(Aa[0],Aa[1],Aa[2])*600/2
            if A[0,w]<>75000:
                Q=np.vstack((Q,A[:,w]))
                w=w+1
        Q=np.delete(Q,0,0)    
        
        if B[h,4]==1:
            a.plot(Q[:,0]+600/2,Q[:,1]+600/2,'g')
        if B[h,4]==2:
            a.plot(Q[:,0]+600/2,Q[:,1]+600/2,'b')
        if B[h,4]==3:
            a.plot(Q[:,0]+600/2,Q[:,1]+600/2,'r')
            
            
        
####################################################################
##### Click a pole
####################################################################    
def click_a_pole(event):
        
    global M,Dstar,D
    x=event.x
    y=event.y
    x=(x-411)*2/620
    y=-(y-400)*2/620
    X=2*x/(1+x**2+y**2)
    Y=2*y/(1+x**2+y**2)
    Z=(-1+x**2+y**2)/(1+x**2+y**2)
    if Z<0:
        X=-X
        Y=-Y
    A=np.dot(np.linalg.inv(M),np.array([X,Y,Z]))
    n=0
    L=np.zeros((3,16**3))                      
                           
                        
    for i in range(-8,9,1):
        for j in range(-8,9,1):
            for k in range(-8,9,1):
                if np.linalg.norm([i,j,k])<>0:
                    if var_uvw.get()==0:
                        Q=np.dot(Dstar,np.array([i,j,k],float))/np.linalg.norm(np.dot(Dstar,np.array([i,j,k],float)))
                        if np.abs(Q[0]-A[0])<0.05 and np.abs(Q[1]-A[1])<0.05 and np.abs(Q[2]-A[2])<0.05:
                            L[:,n]=np.array([i,j,k],float)
                            n=n+1
                           
                    else:
                          
                        Q=np.dot(D,np.array([i,j,k],float))/np.linalg.norm(np.dot(D,np.array([i,j,k],float)))
                        if np.abs(Q[0]-A[0])<0.05 and np.abs(Q[1]-A[1])<0.05 and np.abs(Q[2]-A[2])<0.05:
                            L[:,n]=np.array([i,j,k],float)
                            n=n+1
      

    if np.linalg.norm(L[:,0])<>0:
        if var_hexa.get()==1:
            if var_uvw.get()==1:
                La=(2*L[0,0]-L[1,0])/3
                Lb=(2*L[1,0]-L[0,0])/3
                L[0,0]=La
                L[1,0]=Lb
        #print(L[0,0],L[1,0],L[2,0])
        pole(L[0,0],L[1,0],L[2,0])
        trace()
####################################################################
##### Inclinaison-beta
####################################################################   
def coordinates(event):
#    print(event.x,event.y)
    x=event.x
    y=event.y    
    x=(x-411)*2/620
    y=-(y-400)*2/620    
    long0=90*pi/180
    lat0=0
    r=np.sqrt(x**2+y**2)
    c=2*np.arctan(r)
    longi=(long0+np.arctan2(x*np.sin(c),r*np.cos(lat0)*np.cos(c)-y*np.sin(lat0)*np.sin(c)))*180/pi 
    lat=np.arcsin(np.cos(c)*np.sin(lat0)+y*np.sin(c)*np.cos(lat0)/r)*180/pi
    lat=90-lat
    longi=-longi+180
    if longi>90:
     longi=longi-180
    c=str(np.around(longi,decimals=1))+str(',')+str(np.around(lat,decimals=1))
    coord_label_var.set(c)

####################################################################
##### Fonction rotation autour g
####################################################################



def rotgm():
    global g,M,Dstar
    a = f.add_subplot(111)     
    a.figure.clear()
    thg=-eval(rot_g_entry.get())
    diff1=eval(diff1_entry.get())
    diff2=eval(diff2_entry.get())
    diff3=eval(diff3_entry.get())
    A=np.array([diff1,diff2,diff3])
    Ad=np.dot(Dstar,A)    
    Ap=np.dot(M,Ad)/np.linalg.norm(np.dot(M,Ad))
    #rotproj(thg,Ap[0],Ap[1],Ap[2])
    M=np.dot(Rot(thg,Ap[0],Ap[1],Ap[2]),M)
    trace()    
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_var.set(t)
    g=g-eval(rot_g_entry.get())
    cg.set(g)
    return g,M
    
def rotgp():
    global g,M,D
    a = f.add_subplot(111)     
    a.figure.clear()
    thg=eval(rot_g_entry.get())
    diff1=eval(diff1_entry.get())
    diff2=eval(diff2_entry.get())
    diff3=eval(diff3_entry.get())
    A=np.array([diff1,diff2,diff3])
    Ad=np.dot(Dstar,A)    
    Ap=np.dot(M,Ad)/np.linalg.norm(np.dot(M,Ad))
    #rotproj(thg,Ap[0],Ap[1],Ap[2])
    M=np.dot(Rot(thg,Ap[0],Ap[1],Ap[2]),M)
    trace()    
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_var.set(t)
    g=g+eval(rot_g_entry.get())
    cg.set(g)
    return g,M
####################################################################
##### Fonction principale
####################################################################
def trace():
    global T,x,y,z,axes,axesh,M,trP,a
    a = f.add_subplot(111) 
    a.figure.clear()
    a = f.add_subplot(111)
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
#    print(axesh)
    P=np.zeros((axes.shape[0],2))
    T=np.zeros((axes.shape))
    C=[]

    trace_plan2(trP)
    schmid_trace2(tr_schmid)
    for i in range(0,axes.shape[0]):
        axeshr=np.array([axesh[i,0],axesh[i,1],axesh[i,2]])
#        print(axeshr)
        axeshr=axeshr/np.linalg.norm(axeshr)
        if axesh[i,4]==1:
            C.append('g')
        if axesh[i,4]==2:
            C.append('b')
        if axesh[i,4]==3:
            C.append('r')

               
        T[i,:]=np.dot(M,axeshr)
        P[i,:]=proj(T[i,0],T[i,1],T[i,2])*600/2
        m=np.amax([np.abs(axes[i,0]),np.abs(axes[i,1]),np.abs(axes[i,2])])
        if (np.around(axes[i,0]/m)==axes[i,0]/m) & (np.around(axes[i,1]/m)==axes[i,1]/m) & (np.around(axes[i,2]/m)==axes[i,2]/m):
            if var_hexa.get()==0:
                s=str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(int(axes[i,2]/m))
            if var_hexa.get()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(axes[i,2]/m))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0]/m)-axes[i,1]/m))+str(int(2*(axes[i,1]/m)-axes[i,0]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(3*axes[i,2]/m))+']'
                
        else:
            if var_hexa.get()==0:
                s=str(int(axes[i,0]))+str(int(axes[i,1]))+str(int(axes[i,2]))
            if var_hexa.get()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]))+str(int(axes[i,1]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(axes[i,2]))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0])-axes[i,1]))+str(int(2*(axes[i,1])-axes[i,0]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(3*axes[i,2]))+']'
             
            
        a.annotate(s,(P[i,0]+600/2,P[i,1]+600/2))
       
    if var_carre.get()==0:           
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=eval(size_var.get()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=eval(size_var.get()), facecolors='none', linewidths=1.5)       
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw()
    #a.figure.clear()

    
def princ():
    global T,x,y,z,M,Dstar,D,g,M0,trP,axeshr,nn
    trP=np.zeros((1,5))
    crist() 
    a = f.add_subplot(111)
    a.figure.clear()
    a = f.add_subplot(111)
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    
    diff1=eval(diff1_entry.get())
    diff2=eval(diff2_entry.get())
    diff3=eval(diff3_entry.get())
    tilt=eval(tilt_entry.get())
    inclinaison=eval(inclinaison_entry.get())    
     
    d0=np.array([diff1,diff2,diff3])
    if var_uvw.get()==0: 
       d=np.dot(Dstar,d0)
              
    else:
       d=np.dot(Dstar,d0)
    if diff2==0 and diff1==0:
        normal=np.array([1,0,0])
        ang=pi/2
    else:
        normal=np.array([-d[2],0,d[0]])
        ang=np.arccos(np.dot(d,np.array([0,1,0]))/np.linalg.norm(d))
    
     
    R=np.dot(Rot(-tilt,0,1,0),np.dot(Rot(-inclinaison,0,0,1),Rot(ang*180/pi, normal[0],normal[1],normal[2])))    
    
       
    P=np.zeros((axes.shape[0],2))
    T=np.zeros((axes.shape))
    nn=axes.shape[0]
    C=[]
    for i in range(0,axes.shape[0]):
        
        axeshr=np.array([axesh[i,0],axesh[i,1],axesh[i,2]])
        axeshr=axeshr/np.linalg.norm(axeshr)
        T[i,:]=np.dot(R,axeshr)
        P[i,:]=proj(T[i,0],T[i,1],T[i,2])*600/2
        if col_trace.get()==1:
            C.append('g')
            axesh[i,4]=1
        if col_trace.get()==2:
            C.append('b')
            axesh[i,4]=2
        if col_trace.get()==3:
            C.append('r')
            axesh[i,4]=3
        m=np.amax([np.abs(axes[i,0]),np.abs(axes[i,1]),np.abs(axes[i,2])])
        if (np.around(axes[i,0]/m)==axes[i,0]/m) & (np.around(axes[i,1]/m)==axes[i,1]/m) & (np.around(axes[i,2]/m)==axes[i,2]/m):
            if var_hexa.get()==0:
                s=str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(int(axes[i,2]/m))
            if var_hexa.get()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(axes[i,2]/m))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0]/m)-axes[i,1]/m))+str(int(2*(axes[i,1]/m)-axes[i,0]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(3*axes[i,2]/m))+']'
                
        else:
            if var_hexa.get()==0:
                s=str(int(axes[i,0]))+str(int(axes[i,1]))+str(int(axes[i,2]))
            if var_hexa.get()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]))+str(int(axes[i,1]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(axes[i,2]))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0])-axes[i,1]))+str(int(2*(axes[i,1])-axes[i,0]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(3*axes[i,2]))+']'
                
             
        a.annotate(s,(P[i,0]+600/2,P[i,1]+600/2))
        
    if var_carre.get()==0:           
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=eval(size_var.get()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=eval(size_var.get()), facecolors='none', linewidths=1.5)            
    
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw() 
    
    
    
    x=0             #reinitialise les valeurs de rotations autour de x, y et z quand on clique sur trace
    y=0
    z=0
    g=0
    cx.set(x)
    cy.set(y)
    cz.set(z)
    cg.set(g)
    M=R
    M0=R
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_var.set(t)
   
    return T,x,y,z,g,M,M0 
    

def princ2():
    global T,x,y,z,M,Dstar,D,g,M0,trP,axeshr,nn
    
    trP=np.zeros((1,5))
    a = f.add_subplot(111)
    a.figure.clear()
    a = f.add_subplot(111)
    phi1=eval(phi1_entry.get())
    phi=eval(phi_entry.get())
    phi2=eval(phi2_entry.get())
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    crist()    
    P=np.zeros((axes.shape[0],2))
    T=np.zeros((axes.shape))
    nn=axes.shape[0]
    C=[]    
    
    for i in range(0,axes.shape[0]):
        axeshr=np.array([axesh[i,0],axesh[i,1],axesh[i,2]])
        axeshr=axeshr/np.linalg.norm(axeshr)
        T[i,:]=np.dot(rotation(phi1,phi,phi2),axeshr)
        
        P[i,:]=proj(T[i,0],T[i,1],T[i,2])*600/2
        if col_trace.get()==1:
            C.append('g')
            axesh[i,4]=1
        if col_trace.get()==2:
            C.append('b')
            axesh[i,4]=2
        if col_trace.get()==3:
            C.append('r')
            axesh[i,4]=3

        m=np.amax([np.abs(axes[i,0]),np.abs(axes[i,1]),np.abs(axes[i,2])])
        if (np.around(axes[i,0]/m)==axes[i,0]/m) & (np.around(axes[i,1]/m)==axes[i,1]/m) & (np.around(axes[i,2]/m)==axes[i,2]/m):
            if var_hexa.get()==0:
                s=str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(int(axes[i,2]/m))
            if var_hexa.get()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(axes[i,2]/m))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0]/m)-axes[i,1]/m))+str(int(2*(axes[i,1]/m)-axes[i,0]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(3*axes[i,2]/m))+']'
                
        else:
            if var_hexa.get()==0:
                s=str(int(axes[i,0]))+str(int(axes[i,1]))+str(int(axes[i,2]))
            if var_hexa.get()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]))+str(int(axes[i,1]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(axes[i,2]))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0])-axes[i,1]))+str(int(2*(axes[i,1])-axes[i,0]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(3*axes[i,2]))+']'
                
             
        a.annotate(s,(P[i,0]+600/2,P[i,1]+600/2))
        
               
    if var_carre.get()==0:           
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=eval(size_var.get()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=eval(size_var.get()), facecolors='none', linewidths=1.5)       
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw()  
    
    
    
    x=0             #reinitialise les valeurs de rotations autour de x, y et z quand on clique sur trace
    y=0
    z=0
    cx.set(x)
    cy.set(y)
    cz.set(z)
    M=rotation(phi1,phi,phi2)
    t=str(np.around(phi1,decimals=1))+str(',')+str(np.around(phi,decimals=1))+str(',')+str(np.around(phi2,decimals=1))
    angle_euler_var.set(t)
   
    return T,x,y,z,M 
                    
######################################################################
# GUI
######################################################################
def image_save():
    
    s = asksaveasfile(mode='w', defaultextension=".jpg")    
    if s:    
        f.savefig(s.name)
    #s.close()
####################################################
#fonction d'initialisation
##################################################
def init():
    global x,y,z,cx,cy,cz,cg,g,var_uvw,var_hexa, d_label_var,dmip,angle_euler_var,coord_label_var,angle_var,schmid_var,dhkl_var,trP,tr_schmid,var_carre, col_trace
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    a = f.add_subplot(111)
    a.axis('off')
    a.imshow(img)
    a.figure.canvas.draw()
    #T=np.zeros(axes.shape)
    x=0
    y=0
    z=0
    g=0
    dmip=0
    d_label_var=StringVar()
    d_label_var.set(0)
    cx=StringVar()
    cx.set(x)
    cy=StringVar()
    cy.set(y)    
    cz=StringVar()
    cz.set(z)
    var_uvw=IntVar()
    var_carre=IntVar()
    var_hexa=IntVar()
    cg=StringVar()    
    cg.set(g)
    angle_euler_var=StringVar()
    coord_label_var=StringVar()
    angle_var=StringVar()
    schmid_var=StringVar()
    dhkl_var=DoubleVar()
    trP=np.zeros((1,3))
    tr_schmid=np.zeros((1,3))
    col_trace=IntVar()
    
    return x,y,z,g,cg,cx,cy,cz,var_uvw,d_label_var,dmip,trP,var_carre
   

##############################################################
# fonction pour quitter
#######################################################
def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
#############################################################


root = Tk()
root.wm_title("Stereo-Proj")
root.geometry('1220x798+10+40')
root.configure(bg = '#BDBDBD')
#root.resizable(0,0)
#s=ttk.Style()
#s.theme_use('clam')
style = ttk.Style()
theme = style.theme_use()
default = style.lookup(theme, 'background')



################################################
# Creation d'une zone pour tracer des graphiques
################################################
f = Figure(facecolor='white',figsize=[2,2],dpi=100)

canvas = FigureCanvasTkAgg(f, master=root)

canvas.get_tk_widget().place(x=0,y=0,height=800,width=800)
canvas._tkcanvas.bind('<Button-3>', click_a_pole)
canvas._tkcanvas.bind('<Motion>', coordinates)
canvas.show()
#toolbar = NavigationToolbar2TkAgg( canvas, root )
#toolbar.zoom('off')
#toolbar.update()


###################################################

init()


##############################################
# Boutons
##############################################
phi1_entry = Entry (master=root)
phi1_entry.place(relx=0.82,rely=0.52,relheight=0.03,relwidth=0.04)
phi1_entry.configure(background="white")
phi1_entry.configure(foreground="black")
phi1_entry.configure(highlightcolor="black")
phi1_entry.configure(insertbackground="black")
phi1_entry.configure(selectbackground="#c4c4c4")
phi1_entry.configure(selectforeground="black")


phi_entry = Entry (master=root)
phi_entry.place(relx=0.87,rely=0.52,relheight=0.03,relwidth=0.04)
phi_entry.configure(background="white")
phi_entry.configure(foreground="black")
phi_entry.configure(highlightcolor="black")
phi_entry.configure(insertbackground="black")
phi_entry.configure(selectbackground="#c4c4c4")
phi_entry.configure(selectforeground="black")

label_euler = Label (master=root)
label_euler.place(relx=0.82,rely=0.48,height=19,width=200)
label_euler.configure(activebackground="#cccccc")
label_euler.configure(activeforeground="black")
label_euler.configure(foreground="black")
label_euler.configure(highlightcolor="black")
label_euler.configure(text='''Euler Angles (phi1,phi, phi2)''')

phi2_entry = Entry (master=root)
phi2_entry.place(relx=0.92,rely=0.52,relheight=0.03,relwidth=0.04)
phi2_entry.configure(background="white")
phi2_entry.configure(foreground="black")
phi2_entry.configure(highlightcolor="black")
phi2_entry.configure(insertbackground="black")
phi2_entry.configure(selectbackground="#c4c4c4")
phi2_entry.configure(selectforeground="black")

button_trace = Button (master=root)
button_trace.place(relx=0.87,rely=0.56,height=21,width=49)
button_trace.configure(activebackground="#f9f9f9")
button_trace.configure(activeforeground="black")
button_trace.configure(background="#ff0000")
button_trace.configure(command=princ2)
button_trace.configure(foreground="black")
button_trace.configure(highlightcolor="black")
button_trace.configure(pady="0")
button_trace.configure(text='''PLOT''')



color_button = Checkbutton (master=root)
color_button.place(relx=0.66,rely=0.415,relheight=0.03,relwidth=0.08)
color_button.configure(text='''square/circle''')
color_button.configure(variable=var_carre)

diff1_entry = Entry (master=root)
diff1_entry.place(relx=0.66,rely=0.49,relheight=0.03,relwidth=0.05)

diff1_entry.configure(background="white")
diff1_entry.configure(foreground="black")
diff1_entry.configure(highlightbackground="#e0e0dfdfe3e3")
diff1_entry.configure(highlightcolor="#000000")
diff1_entry.configure(insertbackground="#000000")
diff1_entry.configure(selectbackground="#c4c4c4")
diff1_entry.configure(selectforeground="black")

diff2_entry = Entry (master=root)
diff2_entry.place(relx=0.71,rely=0.49,relheight=0.03,relwidth=0.05)

diff2_entry.configure(background="white")
diff2_entry.configure(foreground="black")
diff2_entry.configure(highlightcolor="black")
diff2_entry.configure(insertbackground="black")
diff2_entry.configure(selectbackground="#c4c4c4")
diff2_entry.configure(selectforeground="black")

diff_label = Label (master=root)
diff_label.place(relx=0.66,rely=0.45,height=19,width=62)
diff_label.configure(activebackground="#cccccc")
diff_label.configure(activeforeground="black")
diff_label.configure(foreground="black")
diff_label.configure(highlightcolor="black")
diff_label.configure(text='''g vector''')


diff3_entry = Entry (master=root)
diff3_entry.place(relx=0.76,rely=0.49,relheight=0.03,relwidth=0.05)

diff3_entry.configure(background="white")
diff3_entry.configure(foreground="black")
diff3_entry.configure(highlightcolor="black")
diff3_entry.configure(insertbackground="black")
diff3_entry.configure(selectbackground="#c4c4c4")
diff3_entry.configure(selectforeground="black")


tilt_entry = Entry (master=root)
tilt_entry.place(relx=0.66,rely=0.55,relheight=0.03,relwidth=0.04)
tilt_entry.configure(background="white")
tilt_entry.configure(foreground="black")
tilt_entry.configure(highlightbackground="#e0e0dfdfe3e3")
tilt_entry.configure(highlightcolor="#000000")
tilt_entry.configure(insertbackground="#000000")
tilt_entry.configure(selectbackground="#c4c4c4")
tilt_entry.configure(selectforeground="black")

inclinaison_entry = Entry (master=root)
inclinaison_entry.place(relx=0.72,rely=0.55,relheight=0.03
    ,relwidth=0.04)
inclinaison_entry.configure(background="white")
inclinaison_entry.configure(foreground="black")
inclinaison_entry.configure(highlightbackground="#e0e0dfdfe3e3")
inclinaison_entry.configure(highlightcolor="#000000")
inclinaison_entry.configure(insertbackground="#000000")
inclinaison_entry.configure(selectbackground="#c4c4c4")
inclinaison_entry.configure(selectforeground="black")

button_trace = Button (master=root)
button_trace.place(relx=0.66,rely=0.61,height=21,width=49)
button_trace.configure(activebackground="#f9f9f9")
button_trace.configure(activeforeground="black")
button_trace.configure(background="#ff0000")
button_trace.configure(command=princ)
button_trace.configure(foreground="black")
button_trace.configure(highlightcolor="black")
button_trace.configure(pady="0")
button_trace.configure(text='''PLOT''')

rx_buttonm = Button (master=root)
rx_buttonm.place(relx=0.83,rely=0.08,height=27,width=27)
rx_buttonm.configure(activebackground="#f9f9f9")
rx_buttonm.configure(activeforeground="black")
rx_buttonm.configure(command=rxm)
rx_buttonm.configure(foreground="black")
rx_buttonm.configure(highlightcolor="black")
rx_buttonm.configure(pady="0")
rx_buttonm.configure(text='''-''')

rz_buttonp = Button (master=root)
rz_buttonp.place(relx=0.92,rely=0.18,height=21,width=17)
rz_buttonp.configure(activebackground="#f9f9f9")
rz_buttonp.configure(activeforeground="black")
rz_buttonp.configure(command=rzp)
rz_buttonp.configure(foreground="black")
rz_buttonp.configure(highlightcolor="black")
rz_buttonp.configure(pady="0")
rz_buttonp.configure(text='''+''')

rz_buttonm = Button (master=root)
rz_buttonm.place(relx=0.84,rely=0.18,height=21,width=13)
rz_buttonm.configure(activebackground="#f9f9f9")
rz_buttonm.configure(activeforeground="black")
rz_buttonm.configure(command=rzm)
rz_buttonm.configure(foreground="black")
rz_buttonm.configure(highlightcolor="black")
rz_buttonm.configure(pady="0")
rz_buttonm.configure(text='''-''')

ry_buttonp = Button (master=root)
ry_buttonp.place(relx=0.92,rely=0.13,height=21,width=17)
ry_buttonp.configure(activebackground="#f9f9f9")
ry_buttonp.configure(activeforeground="black")
ry_buttonp.configure(command=ryp)
ry_buttonp.configure(foreground="black")
ry_buttonp.configure(highlightcolor="black")
ry_buttonp.configure(pady="0")
ry_buttonp.configure(text='''+''')

ry_buttonm = Button (master=root)
ry_buttonm.place(relx=0.84,rely=0.13,height=21,width=13)
ry_buttonm.configure(activebackground="#f9f9f9")
ry_buttonm.configure(activeforeground="black")
ry_buttonm.configure(command=rym)
ry_buttonm.configure(foreground="black")
ry_buttonm.configure(highlightcolor="black")
ry_buttonm.configure(pady="0")
ry_buttonm.configure(text='''-''')

rx_buttonp = Button (master=root)
rx_buttonp.place(relx=0.91,rely=0.08,height=27,width=36)
rx_buttonp.configure(activebackground="#f9f9f9")
rx_buttonp.configure(activeforeground="black")
rx_buttonp.configure(command=rxp)
rx_buttonp.configure(foreground="black")
rx_buttonp.configure(highlightcolor="black")
rx_buttonp.configure(pady="0")
rx_buttonp.configure(text='''+''')

rx_entry = Entry (master=root)
rx_entry.place(relx=0.85,rely=0.08,relheight=0.03,relwidth=0.05)
rx_entry.configure(background="white")
rx_entry.configure(foreground="black")
rx_entry.configure(highlightcolor="black")
rx_entry.configure(insertbackground="black")
rx_entry.configure(selectbackground="#c4c4c4")
rx_entry.configure(selectforeground="black")

ry_entry = Entry (master=root)
ry_entry.place(relx=0.85,rely=0.13,relheight=0.03,relwidth=0.05)
ry_entry.configure(background="white")
ry_entry.configure(foreground="black")
ry_entry.configure(highlightcolor="black")
ry_entry.configure(insertbackground="black")
ry_entry.configure(selectbackground="#c4c4c4")
ry_entry.configure(selectforeground="black")

rz_entry = Entry (master=root)
rz_entry.place(relx=0.85,rely=0.18,relheight=0.03,relwidth=0.05)
rz_entry.configure(background="white")
rz_entry.configure(foreground="black")
rz_entry.configure(highlightcolor="black")
rz_entry.configure(insertbackground="black")
rz_entry.configure(selectbackground="#c4c4c4")
rz_entry.configure(selectforeground="black")

label_euler1 = Label (master=root)
label_euler1.place(relx=0.83,rely=0.03,height=19,width=143)
label_euler1.configure(activebackground="#cccccc")
label_euler1.configure(activeforeground="black")
label_euler1.configure(foreground="black")
label_euler1.configure(highlightcolor="black")
label_euler1.configure(text='''x y and z rotations''')

rx_label = Label (master=root)
rx_label.configure(activebackground="#f9f9f9")
rx_label.configure(activeforeground="black")
rx_label.configure(foreground="black")
rx_label.configure(highlightcolor="black")
rx_label.configure(textvariable=cx)

ry_label = Label (master=root)
ry_label.place(relx=0.95,rely=0.13,height=19,width=34)
ry_label.configure(activebackground="#f9f9f9")
ry_label.configure(activeforeground="black")
ry_label.configure(foreground="black")
ry_label.configure(highlightcolor="black")
ry_label.configure(textvariable=cy)

rz_label = Label (master=root)
rz_label.place(relx=0.95,rely=0.18,height=19,width=34)
rz_label.configure(activebackground="#f9f9f9")
rz_label.configure(activeforeground="black")
rz_label.configure(foreground="black")
rz_label.configure(highlightcolor="black")
rz_label.configure(textvariable=cz)

Cristal_label = Label (master=root)
Cristal_label.place(relx=0.66,rely=0.03,height=19,width=142)
Cristal_label.configure(text='''Crystal parameters''')

a_cristal_label = Label (master=root)
a_cristal_label.place(relx=0.68,rely=0.06,height=19,width=12)
a_cristal_label.configure(text='''a''')

b_cristal_label = Label (master=root)
b_cristal_label.place(relx=0.68,rely=0.1,height=19,width=12)
b_cristal_label.configure(activebackground="#f9f9f9")
b_cristal_label.configure(activeforeground="black")
b_cristal_label.configure(foreground="black")
b_cristal_label.configure(highlightcolor="black")
b_cristal_label.configure(text='''b''')

c_cristal_label = Label (master=root)
c_cristal_label.place(relx=0.68,rely=0.14,height=19,width=11)
c_cristal_label.configure(activebackground="#f9f9f9")
c_cristal_label.configure(activeforeground="black")
c_cristal_label.configure(foreground="black")
c_cristal_label.configure(highlightcolor="black")
c_cristal_label.configure(text='''c''')

alp_cristal_label = Label (master=root)
alp_cristal_label.place(relx=0.67,rely=0.19,height=19,width=32)
alp_cristal_label.configure(activebackground="#f9f9f9")
alp_cristal_label.configure(activeforeground="black")
alp_cristal_label.configure(foreground="black")
alp_cristal_label.configure(highlightcolor="black")
alp_cristal_label.configure(text='''alpha''')

bet_cristal_label = Label (master=root)
bet_cristal_label.place(relx=0.67,rely=0.23,height=19,width=28)
bet_cristal_label.configure(activebackground="#f9f9f9")
bet_cristal_label.configure(activeforeground="black")
bet_cristal_label.configure(foreground="black")
bet_cristal_label.configure(highlightcolor="black")
bet_cristal_label.configure(text='''beta''')

gam_cristal_label = Label (master=root)
gam_cristal_label.place(relx=0.66,rely=0.26,height=19,width=45)
gam_cristal_label.configure(activebackground="#f9f9f9")
gam_cristal_label.configure(activeforeground="black")
gam_cristal_label.configure(foreground="black")
gam_cristal_label.configure(highlightcolor="black")
gam_cristal_label.configure(text='''gamma''')

a_entry = Entry (master=root)
a_entry.place(relx=0.7,rely=0.06,relheight=0.03,relwidth=0.06)
a_entry.configure(background="white")
a_entry.configure(insertbackground="black")

b_entry = Entry (master=root)
b_entry.place(relx=0.7,rely=0.1,relheight=0.03,relwidth=0.06)
b_entry.configure(background="white")
b_entry.configure(foreground="black")
b_entry.configure(highlightcolor="black")
b_entry.configure(insertbackground="black")
b_entry.configure(selectbackground="#c4c4c4")
b_entry.configure(selectforeground="black")

c_entry = Entry (master=root)
c_entry.place(relx=0.7,rely=0.14,relheight=0.03,relwidth=0.06)
c_entry.configure(background="white")
c_entry.configure(foreground="black")
c_entry.configure(highlightcolor="black")
c_entry.configure(insertbackground="black")
c_entry.configure(selectbackground="#c4c4c4")
c_entry.configure(selectforeground="black")

alp_entry = Entry (master=root)
alp_entry.place(relx=0.7,rely=0.18,relheight=0.03,relwidth=0.06)
alp_entry.configure(background="white")
alp_entry.configure(foreground="black")
alp_entry.configure(highlightcolor="black")
alp_entry.configure(insertbackground="black")
alp_entry.configure(selectbackground="#c4c4c4")
alp_entry.configure(selectforeground="black")

bet_entry = Entry (master=root)
bet_entry.place(relx=0.7,rely=0.23,relheight=0.03,relwidth=0.06)
bet_entry.configure(background="white")
bet_entry.configure(foreground="black")
bet_entry.configure(highlightcolor="black")
bet_entry.configure(insertbackground="black")
bet_entry.configure(selectbackground="#c4c4c4")
bet_entry.configure(selectforeground="black")

gam_entry = Entry (master=root)
gam_entry.place(relx=0.7,rely=0.26,relheight=0.03,relwidth=0.06)
gam_entry.configure(background="white")
gam_entry.configure(foreground="black")
gam_entry.configure(highlightcolor="black")
gam_entry.configure(insertbackground="black")
gam_entry.configure(selectbackground="#c4c4c4")
gam_entry.configure(selectforeground="black")

uvw_button = Checkbutton (master=root)
uvw_button.place(relx=0.71,rely=0.61,relheight=0.03,relwidth=0.04)
uvw_button.configure(text='''uvw''')
uvw_button.configure(variable=var_uvw)

hexa_button = Checkbutton (master=root)
hexa_button.place(relx=0.76,rely=0.61,relheight=0.03,relwidth=0.06)
hexa_button.configure(text='''Hexa''')
hexa_button.configure(variable=var_hexa)

e_label = Label (master=root)
e_label.place(relx=0.66,rely=0.31,height=19,width=66)
e_label.configure(text='''max indices''')

e_entry = Entry (master=root)
e_entry.place(relx=0.72,rely=0.31,relheight=0.03,relwidth=0.05)
e_entry.configure(background="white")
e_entry.configure(insertbackground="black")

e2_label = Label (master=root)
e2_label.place(relx=0.68,rely=0.36,height=19,width=12)
e2_label.configure(text='''d''')

dm_button = Button (master=root)
dm_button.place(relx=0.7,rely=0.36,height=21,width=13)
dm_button.configure(activebackground="#f9f9f9")
dm_button.configure(activeforeground="black")
dm_button.configure(command=dm)
dm_button.configure(foreground="black")
dm_button.configure(highlightcolor="black")
dm_button.configure(pady="0")
dm_button.configure(text='''-''')

d_entry = Entry (master=root)
d_entry.place(relx=0.72,rely=0.36,relheight=0.02,relwidth=0.04)
d_entry.configure(background="white")
d_entry.configure(foreground="black")
d_entry.configure(highlightcolor="black")
d_entry.configure(insertbackground="black")
d_entry.configure(selectbackground="#c4c4c4")
d_entry.configure(selectforeground="black")

dp_button = Button (master=root)
dp_button.place(relx=0.76,rely=0.36,height=21,width=17)
dp_button.configure(activebackground="#f9f9f9")
dp_button.configure(activeforeground="black")
dp_button.configure(command=dp)
dp_button.configure(foreground="black")
dp_button.configure(highlightcolor="black")
dp_button.configure(pady="0")
dp_button.configure(text='''+''')

d_label = Label (master=root)
d_label.place(relx=0.73,rely=0.39,height=19,width=16)
d_label.configure(textvariable=d_label_var)

label_addpole = Label (master=root)
label_addpole.place(relx=0.81,rely=0.25,height=19,width=80)
label_addpole.configure(activebackground="#cccccc")
label_addpole.configure(activeforeground="black")
label_addpole.configure(foreground="black")
label_addpole.configure(highlightcolor="black")
label_addpole.configure(text='''Add a pole''')

pole1_entry = Entry (master=root)
pole1_entry.place(relx=0.81,rely=0.29,relheight=0.02,relwidth=0.04)

pole1_entry.configure(background="white")
pole1_entry.configure(foreground="black")
pole1_entry.configure(highlightcolor="black")
pole1_entry.configure(insertbackground="black")
pole1_entry.configure(selectbackground="#c4c4c4")
pole1_entry.configure(selectforeground="black")

pole2_entry = Entry (master=root)
pole2_entry.place(relx=0.86,rely=0.29,relheight=0.02,relwidth=0.04)

pole2_entry.configure(background="white")
pole2_entry.configure(foreground="black")
pole2_entry.configure(highlightcolor="black")
pole2_entry.configure(insertbackground="black")
pole2_entry.configure(selectbackground="#c4c4c4")
pole2_entry.configure(selectforeground="black")

pole3_entry = Entry (master=root)
pole3_entry.place(relx=0.91,rely=0.29,relheight=0.02,relwidth=0.04)

pole3_entry.configure(background="white")
pole3_entry.configure(foreground="black")
pole3_entry.configure(highlightcolor="black")
pole3_entry.configure(insertbackground="black")
pole3_entry.configure(selectbackground="#c4c4c4")
pole3_entry.configure(selectforeground="black")

addpole_button = Button (master=root)
addpole_button.place(relx=0.81,rely=0.34,height=31,width=57)
addpole_button.configure(activebackground="#f9f9f9")
addpole_button.configure(activeforeground="black")
addpole_button.configure(command=addpole)
addpole_button.configure(foreground="black")
addpole_button.configure(highlightcolor="black")
addpole_button.configure(pady="0")
addpole_button.configure(text='''Add''')

undo_addpole_button = Button (master=root)
undo_addpole_button.place(relx=0.855,rely=0.34,height=31,width=11)
undo_addpole_button.configure(command=undo_addpole)
undo_addpole_button.configure(pady="0")
undo_addpole_button.configure(text='''-''')

angle_euler_label = Label (master=root)
angle_euler_label.place(relx=0.84,rely=0.64,height=19,width=128)
angle_euler_label.configure(textvariable=angle_euler_var)

sym_button = Button (master=root)
sym_button.place(relx=0.87,rely=0.34,height=31,width=71)
sym_button.configure(command=addpole_sym)
sym_button.configure(pady="0")
sym_button.configure(text='''Symmetry''')

undo_sym_button = Button (master=root)
undo_sym_button.place(relx=0.925,rely=0.34,height=31,width=11)
undo_sym_button.configure(command=undo_sym)
undo_sym_button.configure(pady="0")
undo_sym_button.configure(text='''-''')

trace_plan_button = Button (master=root)
trace_plan_button.place(relx=0.81,rely=0.39,height=31,width=71)
trace_plan_button.configure(command=trace_addplan)
trace_plan_button.configure(pady="0")
trace_plan_button.configure(text='''Plane''')

trace_plan_sym_button = Button (master=root)
trace_plan_sym_button.place(relx=0.81,rely=0.43,height=31,width=97)
trace_plan_sym_button.configure(command=trace_plan_sym)
trace_plan_sym_button.configure(pady="0")
trace_plan_sym_button.configure(text='''Sym Planes''')

undo_trace_plan_button = Button (master=root)
undo_trace_plan_button.place(relx=0.865,rely=0.39,height=31,width=11)
undo_trace_plan_button.configure(command=undo_trace_addplan)
undo_trace_plan_button.configure(pady="0")
undo_trace_plan_button.configure(text='''-''')

undo_trace_plan_sym_button = Button (master=root)
undo_trace_plan_sym_button.place(relx=0.89,rely=0.43,height=31,width=11)
undo_trace_plan_sym_button.configure(command=undo_trace_plan_sym)
undo_trace_plan_sym_button.configure(pady="0")
undo_trace_plan_sym_button.configure(text='''-''')


trace_schmid = Button (master=root)
trace_schmid.place(relx=0.88,rely=0.39,height=31,width=87)
trace_schmid.configure(command=schmid_trace)
trace_schmid.configure(pady="0")
trace_schmid.configure(text='''Iso-Schmid''')

undo_trace_schmid_button = Button (master=root)
undo_trace_schmid_button.place(relx=0.95,rely=0.39,height=31,width=11)
undo_trace_schmid_button.configure(command=undo_schmid_trace)
undo_trace_schmid_button.configure(pady="0")
undo_trace_schmid_button.configure(text='''-''')

coord_label = Label (master=root)
coord_label.place(relx=0.84,rely=0.95,height=19,width=126)
coord_label.configure(textvariable=coord_label_var)

label_euler2 = Label (master=root)
label_euler2.place(relx=0.84,rely=0.6,height=19,width=80)
label_euler2.configure(activebackground="#cccccc")
label_euler2.configure(activeforeground="black")
label_euler2.configure(foreground="black")
label_euler2.configure(highlightcolor="black")
label_euler2.configure(text='''Phi1,Phi,Phi2''')

label_coord = Label (master=root)
label_coord.place(relx=0.84,rely=0.92,height=19,width=88)
label_coord.configure(activebackground="#cccccc")
label_coord.configure(activeforeground="black")
label_coord.configure(foreground="black")
label_coord.configure(highlightcolor="black")
label_coord.configure(text='''Tilt, inclination''')

c10_entry = Entry (master=root)
c10_entry.place(relx=0.66,rely=0.78,relheight=0.02,relwidth=0.04)
c10_entry.configure(background="white")

c11_entry = Entry (master=root)
c11_entry.place(relx=0.66,rely=0.8,relheight=0.02,relwidth=0.04)
c11_entry.configure(background="white")

c12_entry = Entry (master=root)
c12_entry.place(relx=0.66,rely=0.83,relheight=0.02,relwidth=0.04)
c12_entry.configure(background="white")

angle_button = Button (master=root)
angle_button.place(relx=0.68,rely=0.86,height=21,width=71)
angle_button.configure(command=angle)
angle_button.configure(pady="0")
angle_button.configure(text='''Angle''')


color_trace_vert=Radiobutton(master=root, text="green", variable=col_trace, value=1).place(relx=0.67,rely=0.93,height=21,width=73)
color_trace_bleu=Radiobutton(master=root, text="blue", variable=col_trace, value=2).place(relx=0.73,rely=0.93,height=21,width=71)
color_trace_rouge=Radiobutton(master=root, text="red", variable=col_trace, value=3).place(relx=0.78,rely=0.93,height=21,width=71)

size_var_label= Label (master=root)
size_var_label.place(relx=0.67,rely=0.96,height=19,width=100)
size_var_label.configure(text='''Marker size''')

size_var=Entry(master=root)
size_var.place(relx=0.76,rely=0.96,relheight=0.02,relwidth=0.03)


c20_entry = Entry (master=root)
c20_entry.place(relx=0.7,rely=0.78,relheight=0.02,relwidth=0.04)
c20_entry.configure(background="white")

c21_entry = Entry (master=root)
c21_entry.place(relx=0.7,rely=0.8,relheight=0.02,relwidth=0.04)
c21_entry.configure(background="white")

c22_entry = Entry (master=root)
c22_entry.place(relx=0.7,rely=0.83,relheight=0.02,relwidth=0.04)
c22_entry.configure(background="white")

angle_label = Label (master=root)
angle_label.place(relx=0.69,rely=0.9,height=19,width=58)
angle_label.configure(textvariable=angle_var)

plan1_entry = Entry (master=root)
plan1_entry.place(relx=0.78,rely=0.75,relheight=0.02,relwidth=0.04)

plan1_entry.configure(background="white")

plan2_entry = Entry (master=root)
plan2_entry.place(relx=0.78,rely=0.78,relheight=0.02,relwidth=0.04)

plan2_entry.configure(background="white")

plan3_entry = Entry (master=root)
plan3_entry.place(relx=0.78,rely=0.81,relheight=0.02,relwidth=0.04)

plan3_entry.configure(background="white")

largeur_button = Button (master=root)
largeur_button.place(relx=0.77,rely=0.85,height=21,width=66)
largeur_button.configure(command=largeur_plan)
largeur_button.configure(pady="0")
largeur_button.configure(text='''Width''')

b1_entry = Entry (master=root)
b1_entry.place(relx=0.85,rely=0.75,relheight=0.02,relwidth=0.04)
b1_entry.configure(background="white")

b2_entry = Entry (master=root)
b2_entry.place(relx=0.85,rely=0.78,relheight=0.02,relwidth=0.04)
b2_entry.configure(background="white")

b3_entry = Entry (master=root)
b3_entry.place(relx=0.85,rely=0.8,relheight=0.02,relwidth=0.04)
b3_entry.configure(background="white")

n1_entry = Entry (master=root)
n1_entry.place(relx=0.91,rely=0.75,relheight=0.02,relwidth=0.04)
n1_entry.configure(background="white")

n2_entry = Entry (master=root)
n2_entry.place(relx=0.91,rely=0.78,relheight=0.02,relwidth=0.04)
n2_entry.configure(background="white")

n3_entry = Entry (master=root)
n3_entry.place(relx=0.91,rely=0.8,relheight=0.02,relwidth=0.04)
n3_entry.configure(background="white")

schmid_button = Button (master=root)
schmid_button.place(relx=0.85,rely=0.84,height=21,width=117)
schmid_button.configure(command=schmid)
schmid_button.configure(pady="0")
schmid_button.configure(text='''Schmid factor''')

schmid_label = Label (master=root)
schmid_label.place(relx=0.88,rely=0.87,height=19,width=56)
schmid_label.configure(textvariable=schmid_var)

burgers_label = Label (master=root)
burgers_label.place(relx=0.85,rely=0.72,height=19,width=60)
burgers_label.configure(text='''b vector''')

normal_label = Label (master=root)
normal_label.place(relx=0.91,rely=0.72,height=19,width=60)
normal_label.configure(text='''n vector''')



tilt_label = Label (master=root)
tilt_label.place(relx=0.66,rely=0.53,height=19,width=20)
tilt_label.configure(activebackground="#cccccc")
tilt_label.configure(activeforeground="black")
tilt_label.configure(foreground="black")
tilt_label.configure(highlightcolor="black")
tilt_label.configure(text='''Tilt''')

diff_label2 = Label (master=root)
diff_label2.place(relx=0.72,rely=0.53,height=19,width=65)
diff_label2.configure(activebackground="#cccccc")
diff_label2.configure(activeforeground="black")
diff_label2.configure(foreground="black")
diff_label2.configure(highlightcolor="black")
diff_label2.configure(text='''Inclination''')

rot_gm_button = Button (master=root)
rot_gm_button.place(relx=0.66,rely=0.71,height=21,width=18)
rot_gm_button.configure(activebackground="#f9f9f9")
rot_gm_button.configure(activeforeground="black")
rot_gm_button.configure(command=rotgm)
rot_gm_button.configure(foreground="black")
rot_gm_button.configure(highlightcolor="black")
rot_gm_button.configure(pady="0")
rot_gm_button.configure(text='''-''')

rot_gp_button = Button (master=root)
rot_gp_button.place(relx=0.74,rely=0.71,height=21,width=17)
rot_gp_button.configure(activebackground="#f9f9f9")
rot_gp_button.configure(activeforeground="black")
rot_gp_button.configure(command=rotgp)
rot_gp_button.configure(foreground="black")
rot_gp_button.configure(highlightcolor="black")
rot_gp_button.configure(pady="0")
rot_gp_button.configure(text='''+''')

rot_g_entry = Entry (master=root)
rot_g_entry.place(relx=0.69,rely=0.71,relheight=0.02,relwidth=0.04)

rot_g_entry.configure(background="white")
rot_g_entry.configure(foreground="black")
rot_g_entry.configure(highlightcolor="black")
rot_g_entry.configure(insertbackground="black")
rot_g_entry.configure(selectbackground="#c4c4c4")
rot_g_entry.configure(selectforeground="black")

rot_diff_label = Label (master=root)
rot_diff_label.place(relx=0.66,rely=0.68,height=19,width=190)
rot_diff_label.configure(activebackground="#cccccc")
rot_diff_label.configure(activeforeground="black")
rot_diff_label.configure(foreground="black")
rot_diff_label.configure(highlightcolor="black")
rot_diff_label.configure(text='''Rotation around g vector''')

rg_label = Label (master=root)
rg_label.place(relx=0.76,rely=0.71,height=19,width=36)
rg_label.configure(activebackground="#f9f9f9")
rg_label.configure(activeforeground="black")
rg_label.configure(foreground="black")
rg_label.configure(highlightcolor="black")
rg_label.configure(textvariable=cg)


def dhkl():
    i=eval(pole1_entry.get())
    j=eval(pole2_entry.get())
    k=eval(pole3_entry.get())
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
    d=np.around(1/(np.sqrt(np.dot(np.array([i,j,k]),np.dot(np.linalg.inv(G),np.array([i,j,k]))))), decimals=3)
    dhkl_var.set(d)
    return          

norm_button = Button (master=root)
norm_button.place(relx=0.94,rely=0.32,height=31,width=71)
norm_button.configure(command=dhkl)
norm_button.configure(pady="0")
norm_button.configure(text='''dhkl''')

dhkl_label = Label (master=root)
dhkl_label.place(relx=0.94,rely=0.36,height=19,width=46)
dhkl_label.configure(textvariable=dhkl_var)

menu = Menu(master=root)
filemenu = Menu(menu, tearoff=0)
menu.add_cascade(label="Save", menu=filemenu)

root.config(menu=menu)
filemenu.add_command(label="Save figure", command=image_save) 

def structure(i0):
    global x0, var_hexa, d_label_var, e_entry
    
    d_entry.delete(0,END)
    e_entry.delete(0,END)
    a_entry.delete(0,END)
    a_entry.insert(1,eval(x0[i0][1]))
    b_entry.delete(0,END)    
    b_entry.insert(1,eval(x0[i0][2]))
    c_entry.delete(0,END)    
    c_entry.insert(1,eval(x0[i0][3]))
    alp_entry.delete(0,END)    
    alp_entry.insert(1,eval(x0[i0][4]))
    bet_entry.delete(0,END)    
    bet_entry.insert(1,eval(x0[i0][5]))
    gam_entry.delete(0,END)    
    gam_entry.insert(1,eval(x0[i0][6]))
    if eval(x0[i0][4])==90 and eval(x0[i0][5])==90 and eval(x0[i0][6])==120 :
        var_hexa.set(1)
        e_entry.insert(2,2)
        d_label_var.set(3)
    else:
        var_hexa.set(0)
        d_entry.insert(1,1)
        e_entry.insert(1,1)
    

def createstructure(i):
    return lambda:structure(i)    
    
cristalmenu=Menu(menu,tearoff=0)
menu.add_cascade(label="Structures", menu=cristalmenu)
file_struct=open(os.path.join(os.path.dirname(__file__), 'structure.txt') ,"r")

x0=[]
i=0
for line in file_struct:
    x0.append(map(str, line.split()))
    cristalmenu.add_command(label=x0[i][0], command=createstructure(i))
    i=i+1
  
file_struct.close()        
a_entry.insert(1,1)
b_entry.insert(1,1)
c_entry.insert(1,1)
alp_entry.insert(90,90)
bet_entry.insert(90,90)
gam_entry.insert(90,90)
e_entry.insert(1,1)
d_entry.insert(1,1)
rx_entry.insert(10,10)
ry_entry.insert(10,10)
rz_entry.insert(10,10)

phi1_entry.insert(0,0)
phi_entry.insert(0,0)
phi2_entry.insert(0,0)

size_var.insert(40,40)
col_trace.set(2)
mainloop()     
