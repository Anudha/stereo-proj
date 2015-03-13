#!/usr/bin/python
from __future__ import division
import numpy as np
from PyQt4 import QtGui, QtCore
import sys
import random
import sys
import os
import Image
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import pyplot as plt


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



def color_trace():
        color_trace=1
        if color_trace_bleu.isChecked():
                color_trace=1
        if color_trace_bleu.isChecked():
                color_trace=2
        if color_trace_rouge.isChecked():
                color_trace=3
        return color_trace

def var_uvw():
        var_uvw=0
        if uvw_button.isChecked():
                var_uvw=1
        
        return var_uvw

def var_hexa():
        var_hexa=0
        if hexa_button.isChecked():
                var_hexa=1
        
        return var_hexa

def var_carre():
        var_carre=0
        if style_box.isChecked():
                var_carre=1
        
        return var_carre

####################################################################
##### Fonction cristal
####################################################################
def crist():
    global axes,axesh,D,Dstar,V
    a=np.float(a_entry.text())
    b=np.float(b_entry.text())
    c=np.float(c_entry.text())
    alpha=np.float(alpha_entry.text())
    beta=np.float(beta_entry.text())
    gamma=np.float(gamma_entry.text())
    e=np.int(e_entry.text())
    d2=np.float(d_label_var.text())
    alpha=alpha*pi/180;
    beta=beta*pi/180;
    gamma=gamma*pi/180;
    V=a*b*c*np.sqrt(1-(np.cos(alpha)**2)-(np.cos(beta))**2-(np.cos(gamma))**2+2*b*c*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    D=np.array([[a,b*np.cos(gamma),c*np.cos(beta)],[0,b*np.sin(gamma),  c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)],[0,0,V/(a*b*np.sin(gamma))]])
    Dstar=np.transpose(np.linalg.inv(D))
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])    
    axes=np.zeros(((2*e+1)**3-1,3))
    axesh=np.zeros(((2*e+1)**3-1,5))
    axesh[:,4]=color_trace()
    id=0
    for i in range(-e,e+1):
        for j in range(-e,e+1):
            for k in range(-e,e+1):
                if (i,j,k)!=(0,0,0):
                    d=1/(np.sqrt(np.dot(np.array([i,j,k]),np.dot(np.linalg.inv(G),np.array([i,j,k])))))
                    if d>d2*0.1*np.amax([a,b,c]):
                        if var_uvw()==0:                    
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
    a=figure.add_subplot(111)    
    a.figure.clear()
    a=figure.add_subplot(111)
    dmip=dmip-np.float(d_entry.text())
    d_label_var.setText(str(dmip))
    crist()
    trace()
    
    return dmip
def dp():
    global dmip
    a=figure.add_subplot(111)    
    a.figure.clear()
    a=figure.add_subplot(111)
    dmip=dmip+np.float(d_entry.text())
    d_label_var.setText(str(dmip))
    crist()    
    trace()
    
    return dmip 
    

####################################################################
##### Fonction largeur plan
####################################################################
class Window(QtGui.QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
 
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
 
         
        self.toolbar = NavigationToolbar(self.canvas, self)
        
                
        gridLayout = QtGui.QGridLayout()
        self.lineEdit = QtGui.QLineEdit()
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        gridLayout.addWidget(self.lineEdit, 1, 0, 1, 1)
        self.lineEdit_2 = QtGui.QLineEdit()
        gridLayout.addWidget(self.lineEdit_2, 1, 1, 1, 1)
        self.lineEdit_3 = QtGui.QLineEdit()
        gridLayout.addWidget(self.lineEdit_3, 1, 2, 1, 1)
        
        self.buttonBox = QtGui.QDialogButtonBox()
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        gridLayout.addWidget(self.buttonBox, 2, 0, 1, 3)
        gridLayout.addWidget(self.canvas, 0, 0, 1, 3)
        gridLayout.addWidget(self.toolbar, 3, 0, 1, 3)
        
        self.buttonBox.rejected.connect(self.close)
        self.buttonBox.accepted.connect(self.plot)
              
        
        self.setLayout(gridLayout)
 
    def home(self):
        self.toolbar.home()
    def zoom(self):
        self.toolbar.zoom()
    def pan(self):
        self.toolbar.pan()
         
    def plot(self):
        ''' plot some random stuff '''
        global D, Dstar, M 
    
    
        plan1=np.float(self.lineEdit.text())   
        plan2=np.float(self.lineEdit_2.text())
        plan3=np.float(self.lineEdit_3.text())
        
        n=np.array([plan1,plan2,plan3])
        if var_uvw()==0: 
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
        
#        plt.plot(range(-40,41,2),la[0,:])
#        plt.xlabel('Tilt angle')
#        plt.ylabel('Apparent width')
#        plt.grid(True)
#        plt.show()
        
        ax = self.figure.add_subplot(111)
        ax.hold(False)
        ax.plot(range(-40,41,2),la[0,:])
        self.canvas.draw()

    
####################################################################
##### Fonction facteur de Schmid
####################################################################
def schmid():
    global D, Dstar,M
    b1=np.float(b1_entry.text())   
    b2=np.float(b2_entry.text())
    b3=np.float(b3_entry.text())
    n1=np.float(n1_entry.text())   
    n2=np.float(n2_entry.text())
    n3=np.float(n3_entry.text())
    n=np.array([n1,n2,n3])
    b=np.array([b1,b2,b3])
    if var_uvw()==0: 
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
    
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
        
    tr_schmid=np.vstack((tr_schmid,np.array([pole1,pole2,pole3])))
    trace()
    
def undo_schmid_trace():
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
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
        
        if var_uvw()==0: 
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
    a = figure.add_subplot(111)     
    a.figure.clear()
    thx=np.float(rx_entry.text())
  
    M=np.dot(Rot(thx,1,0,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_label.setText(t)
    x=x+np.float(rx_entry.text())
    rx_label.setText(str(x))
    return x,M
    
    
def rxm():
    global x,M,a
    a = figure.add_subplot(111)     
    a.figure.clear()            
    thx=-np.float(rx_entry.text())
    
    M=np.dot(Rot(thx,1,0,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_label.setText(t)
    x=x-np.float(rx_entry.text())
    rx_label.setText(str(x))
    return x,M

def ryp():
    global y,M,a
    a = figure.add_subplot(111)     
    a.figure.clear()
    thy=np.float(ry_entry.text())
    
    M=np.dot(Rot(thy,0,1,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_label.setText(t)
    y=y+np.float(ry_entry.text())
    ry_label.setText(str(y))
    return y,M
    
def rym():
    global y,M,a
    a = figure.add_subplot(111)     
    a.figure.clear()
    thy=-np.float(ry_entry.text())
    
    M=np.dot(Rot(thy,0,1,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_label.setText(t)
    y=y-np.float(ry_entry.text())
    ry_label.setText(str(y))
    return y,M


def rzp():
    global z,M,a
    a = figure.add_subplot(111)     
    a.figure.clear()
    thz=np.float(rz_entry.text())
    
    M=np.dot(Rot(thz,0,0,1),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_label.setText(t)
    z=z+np.float(rz_entry.text())
    rz_label.setText(str(z))
    return z,M
    
def rzm():
    global z,M,a
    a = figure.add_subplot(111)     
    a.figure.clear()
    thz=-np.float(rz_entry.text())
    
    M=np.dot(Rot(thz,0,0,1),M)
    trace()
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    angle_euler_label.setText(t)
    z=z-np.float(rz_entry.text())
    rz_label.setText(str(z))
    return z,M

####################################################################
##### Fonction ajouter un pole
####################################################################
def pole(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar
    
#    fp=figure.add_subplot(111)
#    
#       
#    Pp=np.zeros((1,2),float)
    if var_hexa()==1:
        if var_uvw()==1:
            pole1a=2*pole1+pole2
            pole2a=2*pole2+pole1
            pole1=pole1a
            pole2=pole2a
    
    Gs=np.array([pole1,pole2,pole3],float)
    #print(Gs)
    if var_uvw()==0:                    
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
    if var_uvw()==0 :
        axesh=np.vstack((axesh,np.array([Gsh[0],Gsh[1],Gsh[2],0,color_trace()])))
        axesh=np.vstack((axesh,np.array([-Gsh[0],-Gsh[1],-Gsh[2],0,color_trace()])))
    else:
        axesh=np.vstack((axesh,np.array([Gsh[0],Gsh[1],Gsh[2],1,color_trace()])))
        axesh=np.vstack((axesh,np.array([-Gsh[0],-Gsh[1],-Gsh[2],1,color_trace()])))
    return axes,axesh,T

def undo_pole(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar
    
#    fp=figure.add_subplot(111)
#    
#       
#    Pp=np.zeros((1,2),float)
    if var_hexa()==1:
        if var_uvw()==1:
            pole1a=2*pole1+pole2
            pole2a=2*pole2+pole1
            pole1=pole1a
            pole2=pole2a
    
    Gs=np.array([pole1,pole2,pole3],float)
    #print(Gs)
    if var_uvw()==0:                    
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
    if var_uvw()==0 :
        axesh=np.delete(axesh,ind,0)
        axesh=np.delete(axesh,indm,0)
    else:
        axesh=np.delete(axesh,ind,0)
        axesh=np.delete(axesh,indm,0)
    return axes,axesh,T

    
def d(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar
    a=np.float(a_entry.text())
    b=np.float(b_entry.text())
    c=np.float(c_entry.text())
    alpha=np.float(alpha_entry.text())
    beta=np.float(beta_entry.text())
    gamma=np.float(gamma_entry.text())
    alpha=alpha*pi/180;
    beta=beta*pi/180;
    gamma=gamma*pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])    
    ds=(np.sqrt(np.dot(np.array([pole1,pole2,pole3]),np.dot(np.linalg.inv(G),np.array([pole1,pole2,pole3])))))
    return ds
def addpole_sym():
    global M,axes,axesh,T,V,D,Dstar,G    
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
    a=np.float(a_entry.text())
    b=np.float(b_entry.text())
    c=np.float(c_entry.text())
    alpha=np.float(alpha_entry.text())
    beta=np.float(beta_entry.text())
    gamma=np.float(gamma_entry.text())
    alpha=alpha*pi/180;
    beta=beta*pi/180;
    gamma=gamma*pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])     
    v=d(pole1,pole2,pole3)
    
    pole(pole1,pole2,pole3)
    if np.abs(alpha-pi/2)<0.001 and np.abs(beta-pi/2)<0.001 and np.abs(gamma-2*pi/3)<0.001:
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
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
    a=np.float(a_entry.text())
    b=np.float(b_entry.text())
    c=np.float(c_entry.text())
    alpha=np.float(alpha_entry.text())
    beta=np.float(beta_entry.text())
    gamma=np.float(gamma_entry.text())
    alpha=alpha*pi/180;
    beta=beta*pi/180;
    gamma=gamma*pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])     
    v=d(pole1,pole2,pole3)
    
    undo_pole(pole1,pole2,pole3)
    if np.abs(alpha-pi/2)<0.001 and np.abs(beta-pi/2)<0.001 and np.abs(gamma-2*pi/3)<0.001:
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
    
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
    pole(pole1,pole2,pole3)
    trace()
    
def undo_addpole():
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid,nn
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
    undo_pole(pole1,pole2,pole3)
    
    trace()
####################################################################
##### Fonction tracer plan
####################################################################
def trace_plan(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar,trP
    
    
    pole_i=0
    pole_c=color_trace()
    
    if var_hexa()==1:
        if var_uvw()==1:
            pole1=2*np.float(pole1_entry.text())+np.float(pole2_entry.text())
            pole2=2*np.float(pole2_entry.text())+np.float(pole1_entry.text())
            pole3=np.float(pole3_entry.text())
            pole_i=1
    
        
    trP=np.vstack((trP,np.array([pole1,pole2,pole3,pole_i,pole_c])))
    b=np.ascontiguousarray(trP).view(np.dtype((np.void, trP.dtype.itemsize * trP.shape[1])))
    
    trP=np.unique(b).view(trP.dtype).reshape(-1, trP.shape[1])
    #print(trP)
    
    
def trace_addplan():
    global M,axes,axesh,T,V,D,Dstar,trP
    
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
    
    trace_plan(pole1,pole2,pole3)
    trace()
    
def undo_trace_addplan():
    global M,axes,axesh,T,V,D,Dstar,trP
    
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
    
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
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
    a=np.float(a_entry.text())
    b=np.float(b_entry.text())
    c=np.float(c_entry.text())
    alpha=np.float(alpha_entry.text())
    beta=np.float(beta_entry.text())
    gamma=np.float(gamma_entry.text())
    alpha=alpha*pi/180;
    beta=beta*pi/180;
    gamma=gamma*pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])     
    v=d(pole1,pole2,pole3)
    
    trace_plan(pole1,pole2,pole3)
    if np.abs(alpha-pi/2)<0.001 and np.abs(beta-pi/2)<0.001 and np.abs(gamma-2*pi/3)<0.001:
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
    pole1=np.float(pole1_entry.text())
    pole2=np.float(pole2_entry.text())
    pole3=np.float(pole3_entry.text())
    a=np.float(a_entry.text())
    b=np.float(b_entry.text())
    c=np.float(c_entry.text())
    alpha=np.float(alpha_entry.text())
    beta=np.float(beta_entry.text())
    gamma=np.float(gamma_entry.text())
    alpha=alpha*pi/180;
    beta=beta*pi/180;
    gamma=gamma*pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])     
    v=d(pole1,pole2,pole3)
    
    undo_trace_plan(pole1,pole2,pole3)
    if np.abs(alpha-pi/2)<0.001 and np.abs(beta-pi/2)<0.001 and np.abs(gamma-2*pi/3)<0.001:
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
    if event.button==3:
            x=event.x
            y=event.y
            x=(x-411)*2/620
            y=(y-400)*2/620
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
                            if var_uvw()==0:
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
                if var_hexa()==1:
                    if var_uvw()==1:
                        La=(2*L[0,0]-L[1,0])/3
                        Lb=(2*L[1,0]-L[0,0])/3
                        L[0,0]=La
                        L[1,0]=Lb
                #print(L[0,0],L[1,0],L[2,0])
                pole(L[0,0],L[1,0],L[2,0])
                trace()
        

####################################################################
##### Fonction rotation autour g
####################################################################


def rotgm():
    global g,M,Dstar
    a = figure.add_subplot(111)     
    a.figure.clear()
    thg=-np.float(rot_g_entry.text())
    diff1=np.float(diff1_entry.text())
    diff2=np.float(diff2_entry.text())
    diff3=np.float(diff3_entry.text())
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
    
    angle_euler_label.setText(t)
    g=g-np.float(rot_g_entry.text())
    rg_label.setText(str(g))
    return g,M
    
def rotgp():
    global g,M,D
    a = figure.add_subplot(111)     
    a.figure.clear()
    thg=np.float(rot_g_entry.text())
    diff1=np.float(diff1_entry.text())
    diff2=np.float(diff2_entry.text())
    diff3=np.float(diff3_entry.text())
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
    
    angle_euler_label.setText(t)
    g=g+np.float(rot_g_entry.text())
    rg_label.setText(str(g))
    return g,M
####################################################################
##### Fonction principale
####################################################################
def trace():
    global T,x,y,z,axes,axesh,M,trP,a
    a = figure.add_subplot(111) 
    a.figure.clear()
    a = figure.add_subplot(111)
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
            if var_hexa()==0:
                s=str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(int(axes[i,2]/m))
            if var_hexa()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(axes[i,2]/m))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0]/m)-axes[i,1]/m))+str(int(2*(axes[i,1]/m)-axes[i,0]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(3*axes[i,2]/m))+']'
                
        else:
            if var_hexa()==0:
                s=str(int(axes[i,0]))+str(int(axes[i,1]))+str(int(axes[i,2]))
            if var_hexa()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]))+str(int(axes[i,1]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(axes[i,2]))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0])-axes[i,1]))+str(int(2*(axes[i,1])-axes[i,0]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(3*axes[i,2]))+']'
             
            
        a.annotate(s,(P[i,0]+600/2,P[i,1]+600/2))
       
    if var_carre()==0:           
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=np.float(size_var.text()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=np.float(size_var.text()), facecolors='none', linewidths=1.5)       
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw()
    #a.figure.clear()

    
def princ():
    global T,x,y,z,M,Dstar,D,g,M0,trP,axeshr,nn
    trP=np.zeros((1,5))
    crist() 
    a = figure.add_subplot(111)
    a.figure.clear()
    a = figure.add_subplot(111)
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    
    diff1=np.float(diff1_entry.text())
    diff2=np.float(diff2_entry.text())
    diff3=np.float(diff3_entry.text())
    tilt=np.float(tilt_entry.text())
    inclinaison=np.float(inclinaison_entry.text())    
     
    d0=np.array([diff1,diff2,diff3])
    if var_uvw()==0: 
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
        if color_trace()==1:
            C.append('g')
            axesh[i,4]=1
        if color_trace()==2:
            C.append('b')
            axesh[i,4]=2
        if color_trace()==3:
            C.append('r')
            axesh[i,4]=3
        m=np.amax([np.abs(axes[i,0]),np.abs(axes[i,1]),np.abs(axes[i,2])])
        if (np.around(axes[i,0]/m)==axes[i,0]/m) & (np.around(axes[i,1]/m)==axes[i,1]/m) & (np.around(axes[i,2]/m)==axes[i,2]/m):
            if var_hexa()==0:
                s=str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(int(axes[i,2]/m))
            if var_hexa()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(axes[i,2]/m))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0]/m)-axes[i,1]/m))+str(int(2*(axes[i,1]/m)-axes[i,0]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(3*axes[i,2]/m))+']'
                
        else:
            if var_hexa()==0:
                s=str(int(axes[i,0]))+str(int(axes[i,1]))+str(int(axes[i,2]))
            if var_hexa()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]))+str(int(axes[i,1]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(axes[i,2]))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0])-axes[i,1]))+str(int(2*(axes[i,1])-axes[i,0]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(3*axes[i,2]))+']'
                
             
        a.annotate(s,(P[i,0]+600/2,P[i,1]+600/2))
        
    if var_carre()==0:           
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=np.float(size_var.text()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=np.float(size_var.text()), facecolors='none', linewidths=1.5)            
    
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw() 
    
    
    
    x=0             #reinitialise les valeurs de rotations autour de x, y et z quand on clique sur trace
    y=0
    z=0
    g=0
    rx_label.setText(str(x))
    ry_label.setText(str(y))
    rz_label.setText(str(z))
    rg_label.setText(str(g))
    M=R
    M0=R
    phir=np.arccos(M[2,2])*180/pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    angle_euler_label.setText(t)
   
    return T,x,y,z,g,M,M0 
    

def princ2():
    global T,x,y,z,M,Dstar,D,g,M0,trP,axeshr,nn
    
    trP=np.zeros((1,5))
    a = figure.add_subplot(111)
    a.figure.clear()
    a = figure.add_subplot(111)
    phi1=np.float(phi1_entry.text())
    phi=np.float(phi_entry.text())
    phi2=np.float(phi2_entry.text())
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
        if color_trace()==1:
            C.append('g')
            axesh[i,4]=1
        if color_trace()==2:
            C.append('b')
            axesh[i,4]=2
        if color_trace()==3:
            C.append('r')
            axesh[i,4]=3

        m=np.amax([np.abs(axes[i,0]),np.abs(axes[i,1]),np.abs(axes[i,2])])
        if (np.around(axes[i,0]/m)==axes[i,0]/m) & (np.around(axes[i,1]/m)==axes[i,1]/m) & (np.around(axes[i,2]/m)==axes[i,2]/m):
            if var_hexa()==0:
                s=str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(int(axes[i,2]/m))
            if var_hexa()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]/m))+str(int(axes[i,1]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(axes[i,2]/m))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0]/m)-axes[i,1]/m))+str(int(2*(axes[i,1]/m)-axes[i,0]/m))+str(-int(axes[i,1]/m)-int(axes[i,0]/m))+str(int(3*axes[i,2]/m))+']'
                
        else:
            if var_hexa()==0:
                s=str(int(axes[i,0]))+str(int(axes[i,1]))+str(int(axes[i,2]))
            if var_hexa()==1:
                if axesh[i,3]==0:
                    s='('+str(int(axes[i,0]))+str(int(axes[i,1]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(axes[i,2]))+')'
                if axesh[i,3]==1:
                    s='['+str(int(2*(axes[i,0])-axes[i,1]))+str(int(2*(axes[i,1])-axes[i,0]))+str(-int(axes[i,1])-int(axes[i,0]))+str(int(3*axes[i,2]))+']'
                
             
        a.annotate(s,(P[i,0]+600/2,P[i,1]+600/2))
        
               
    if var_carre()==0:           
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=np.float(size_var.text()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=np.float(size_var.text()), facecolors='none', linewidths=1.5)       
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw()  
    
    
    
    x=0             #reinitialise les valeurs de rotations autour de x, y et z quand on clique sur trace
    y=0
    z=0
    rx_label.setText(str(x))
    ry_label.setText(str(y))
    rz_label.setText(str(z))
    M=rotation(phi1,phi,phi2)
    t=str(np.around(phi1,decimals=1))+str(',')+str(np.around(phi,decimals=1))+str(',')+str(np.around(phi2,decimals=1))
    angle_euler_label.setText(t)
   
    return T,x,y,z,M 

def dhkl():
    i=np.float(pole1_entry.text())
    j=np.float(pole2_entry.text())
    k=np.float(pole3_entry.text())
    a=np.float(a_entry.text())
    b=np.float(b_entry.text())
    c=np.float(c_entry.text())
    alp=np.float(alpha_entry.text())
    bet=np.float(beta_entry.text())
    gam=np.float(gamma_entry.text())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
    d=np.around(1/(np.sqrt(np.dot(np.array([i,j,k]),np.dot(np.linalg.inv(G),np.array([i,j,k]))))), decimals=3)
    dhkl_label.setText(str(d))
    return          

################################################################""
###############                      GUI
##############################################################"
try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

app = QtGui.QApplication(sys.argv)
MainWindow=QtGui.QMainWindow()
MainWindow.setWindowTitle('Application simple')

centralwidget = QtGui.QWidget(MainWindow)
centralwidget.setObjectName(_fromUtf8("centralwidget"))
gridLayout_7 = QtGui.QGridLayout(centralwidget)
gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
gridLayout = QtGui.QGridLayout()
gridLayout.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
gridLayout.setObjectName(_fromUtf8("gridLayout"))
groupBox_5 = QtGui.QGroupBox(centralwidget)
groupBox_5.setMaximumSize(QtCore.QSize(150, 16777215))
groupBox_5.setObjectName(_fromUtf8("groupBox_5"))
gridLayout_6 = QtGui.QGridLayout(groupBox_5)
gridLayout_6.setObjectName(_fromUtf8("gridLayout_6"))
lab_coord = QtGui.QLabel(groupBox_5)
lab_coord.setObjectName(_fromUtf8("lab_coord"))
gridLayout_6.addWidget(lab_coord, 4, 0, 1, 3)
button_trace2 = QtGui.QPushButton(groupBox_5)
button_trace2.setObjectName(_fromUtf8("button_trace2"))
gridLayout_6.addWidget(button_trace2, 1, 0, 1, 2)
phi1_entry = QtGui.QLineEdit(groupBox_5)
phi1_entry.setObjectName(_fromUtf8("phi1_entry"))
gridLayout_6.addWidget(phi1_entry, 0, 0, 1, 1)
angle_euler_label = QtGui.QLabel(groupBox_5)
angle_euler_label.setText(_fromUtf8(""))
angle_euler_label.setObjectName(_fromUtf8("angle_euler_label"))
gridLayout_6.addWidget(angle_euler_label, 3, 0, 1, 3)
phi2_entry = QtGui.QLineEdit(groupBox_5)
phi2_entry.setObjectName(_fromUtf8("phi2_entry"))
gridLayout_6.addWidget(phi2_entry, 0, 2, 1, 1)
phi_entry = QtGui.QLineEdit(groupBox_5)
phi_entry.setObjectName(_fromUtf8("phi_entry"))
gridLayout_6.addWidget(phi_entry, 0, 1, 1, 1)
lab_euler2 = QtGui.QLabel(groupBox_5)
lab_euler2.setObjectName(_fromUtf8("lab_euler2"))
gridLayout_6.addWidget(lab_euler2, 2, 0, 1, 3)
coord_label = QtGui.QLabel(groupBox_5)
coord_label.setText(_fromUtf8(""))
coord_label.setObjectName(_fromUtf8("coord_label"))
gridLayout_6.addWidget(coord_label, 5, 0, 1, 2)
gridLayout.addWidget(groupBox_5, 3, 3, 1, 1)
groupBox_3 = QtGui.QGroupBox(centralwidget)
groupBox_3.setMaximumSize(QtCore.QSize(150, 16777215))
groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
gridLayout_4 = QtGui.QGridLayout(groupBox_3)
gridLayout_4.setMargin(5)
gridLayout_4.setSpacing(5)
gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
rx_buttonp = QtGui.QPushButton(groupBox_3)
rx_buttonp.setObjectName(_fromUtf8("rx_buttonp"))
gridLayout_4.addWidget(rx_buttonp, 3, 3, 1, 1)
rx_label = QtGui.QLabel(groupBox_3)
rx_label.setObjectName(_fromUtf8("rx_label"))
gridLayout_4.addWidget(rx_label, 3, 0, 1, 1)
rz_buttonm = QtGui.QPushButton(groupBox_3)
rz_buttonm.setObjectName(_fromUtf8("rz_buttonm"))
gridLayout_4.addWidget(rz_buttonm, 11, 1, 1, 1)
rz_buttonp = QtGui.QPushButton(groupBox_3)
rz_buttonp.setObjectName(_fromUtf8("rz_buttonp"))
gridLayout_4.addWidget(rz_buttonp, 11, 3, 1, 1)
rz_entry = QtGui.QLineEdit(groupBox_3)
rz_entry.setObjectName(_fromUtf8("rz_entry"))
gridLayout_4.addWidget(rz_entry, 11, 2, 1, 1)
ry_buttonp = QtGui.QPushButton(groupBox_3)
ry_buttonp.setObjectName(_fromUtf8("ry_buttonp"))
gridLayout_4.addWidget(ry_buttonp, 7, 3, 1, 1)
ry_entry = QtGui.QLineEdit(groupBox_3)
ry_entry.setObjectName(_fromUtf8("ry_entry"))
gridLayout_4.addWidget(ry_entry, 7, 2, 1, 1)
rz_label = QtGui.QLabel(groupBox_3)
rz_label.setObjectName(_fromUtf8("rz_label"))
gridLayout_4.addWidget(rz_label, 11, 0, 1, 1)
ry_label = QtGui.QLabel(groupBox_3)
ry_label.setObjectName(_fromUtf8("ry_label"))
gridLayout_4.addWidget(ry_label, 7, 0, 1, 1)
rx_buttonm = QtGui.QPushButton(groupBox_3)
rx_buttonm.setObjectName(_fromUtf8("rx_buttonm"))
gridLayout_4.addWidget(rx_buttonm, 3, 1, 1, 1)
rx_entry = QtGui.QLineEdit(groupBox_3)
rx_entry.setObjectName(_fromUtf8("rx_entry"))
gridLayout_4.addWidget(rx_entry, 3, 2, 1, 1)
ry_buttonm = QtGui.QPushButton(groupBox_3)
ry_buttonm.setObjectName(_fromUtf8("ry_buttonm"))
gridLayout_4.addWidget(ry_buttonm, 7, 1, 1, 1)
rx_label_2 = QtGui.QLabel(groupBox_3)
rx_label_2.setText(_fromUtf8(""))
rx_label_2.setObjectName(_fromUtf8("rx_label_2"))
gridLayout_4.addWidget(rx_label_2, 5, 2, 1, 1)
ry_label_2 = QtGui.QLabel(groupBox_3)
ry_label_2.setText(_fromUtf8(""))
ry_label_2.setObjectName(_fromUtf8("ry_label_2"))
gridLayout_4.addWidget(ry_label_2, 8, 2, 1, 1)
rz_label_2 = QtGui.QLabel(groupBox_3)
rz_label_2.setText(_fromUtf8(""))
rz_label_2.setObjectName(_fromUtf8("rz_label_2"))
gridLayout_4.addWidget(rz_label_2, 12, 2, 1, 1)
gridLayout.addWidget(groupBox_3, 1, 3, 1, 1)
groupBox_4 = QtGui.QGroupBox(centralwidget)
groupBox_4.setMaximumSize(QtCore.QSize(150, 16777215))
groupBox_4.setObjectName(_fromUtf8("groupBox_4"))
gridLayout_5 = QtGui.QGridLayout(groupBox_4)
gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
undo_trace_schmid = QtGui.QPushButton(groupBox_4)
undo_trace_schmid.setObjectName(_fromUtf8("undo_trace_schmid"))
gridLayout_5.addWidget(undo_trace_schmid, 5, 2, 1, 1)
norm_button = QtGui.QPushButton(groupBox_4)
norm_button.setObjectName(_fromUtf8("norm_button"))
gridLayout_5.addWidget(norm_button, 7, 0, 1, 1)
undo_sym_button = QtGui.QPushButton(groupBox_4)
undo_sym_button.setObjectName(_fromUtf8("undo_sym_button"))
gridLayout_5.addWidget(undo_sym_button, 2, 2, 1, 1)
trace_schmid_button = QtGui.QPushButton(groupBox_4)
trace_schmid_button.setObjectName(_fromUtf8("trace_schmid_button"))
gridLayout_5.addWidget(trace_schmid_button, 5, 0, 1, 2)
addpole_button = QtGui.QPushButton(groupBox_4)
addpole_button.setObjectName(_fromUtf8("addpole_button"))
gridLayout_5.addWidget(addpole_button, 1, 0, 1, 2)
dhkl_label = QtGui.QLabel(groupBox_4)
dhkl_label.setText(_fromUtf8(""))
dhkl_label.setObjectName(_fromUtf8("dhkl_label"))
gridLayout_5.addWidget(dhkl_label, 7, 1, 1, 2)
trace_plan_button = QtGui.QPushButton(groupBox_4)
trace_plan_button.setObjectName(_fromUtf8("trace_plan_button"))
gridLayout_5.addWidget(trace_plan_button, 3, 0, 1, 2)
pole3_entry = QtGui.QLineEdit(groupBox_4)
pole3_entry.setObjectName(_fromUtf8("pole3_entry"))
gridLayout_5.addWidget(pole3_entry, 0, 2, 1, 1)
undo_trace_plan_sym_button = QtGui.QPushButton(groupBox_4)
undo_trace_plan_sym_button.setObjectName(_fromUtf8("undo_trace_plan_sym_button"))
gridLayout_5.addWidget(undo_trace_plan_sym_button, 4, 2, 1, 1)
trace_plan_sym_button = QtGui.QPushButton(groupBox_4)
trace_plan_sym_button.setObjectName(_fromUtf8("trace_plan_sym_button"))
gridLayout_5.addWidget(trace_plan_sym_button, 4, 0, 1, 2)
sym_button = QtGui.QPushButton(groupBox_4)
sym_button.setObjectName(_fromUtf8("sym_button"))
gridLayout_5.addWidget(sym_button, 2, 0, 1, 2)
undo_trace_plan_button = QtGui.QPushButton(groupBox_4)
undo_trace_plan_button.setObjectName(_fromUtf8("undo_trace_plan_button"))
gridLayout_5.addWidget(undo_trace_plan_button, 3, 2, 1, 1)
undo_addpole_button = QtGui.QPushButton(groupBox_4)
undo_addpole_button.setObjectName(_fromUtf8("undo_addpole_button"))
gridLayout_5.addWidget(undo_addpole_button, 1, 2, 1, 1)
pole2_entry = QtGui.QLineEdit(groupBox_4)
pole2_entry.setObjectName(_fromUtf8("pole2_entry"))
gridLayout_5.addWidget(pole2_entry, 0, 1, 1, 1)
pole1_entry = QtGui.QLineEdit(groupBox_4)
pole1_entry.setObjectName(_fromUtf8("pole1_entry"))
gridLayout_5.addWidget(pole1_entry, 0, 0, 1, 1)
largeur_button = QtGui.QPushButton(groupBox_4)
largeur_button.setObjectName(_fromUtf8("largeur_button"))
gridLayout_5.addWidget(largeur_button, 6, 0, 1, 3)
gridLayout.addWidget(groupBox_4, 2, 3, 1, 1)
groupBox_2 = QtGui.QGroupBox(centralwidget)
groupBox_2.setEnabled(True)
groupBox_2.setMaximumSize(QtCore.QSize(150, 16777215))
groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
gridLayout_3 = QtGui.QGridLayout(groupBox_2)
gridLayout_3.setMargin(5)
gridLayout_3.setSpacing(5)
gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
inclinaison_entry = QtGui.QLineEdit(groupBox_2)
inclinaison_entry.setObjectName(_fromUtf8("inclinaison_entry"))
gridLayout_3.addWidget(inclinaison_entry, 10, 3, 1, 1)
button_trace = QtGui.QPushButton(groupBox_2)
button_trace.setObjectName(_fromUtf8("button_trace"))
gridLayout_3.addWidget(button_trace, 11, 0, 1, 3)
diff1_entry = QtGui.QLineEdit(groupBox_2)
diff1_entry.setObjectName(_fromUtf8("diff1_entry"))
gridLayout_3.addWidget(diff1_entry, 3, 0, 1, 1)
inclination_label = QtGui.QLabel(groupBox_2)
inclination_label.setObjectName(_fromUtf8("inclination_label"))
gridLayout_3.addWidget(inclination_label, 10, 0, 1, 3)
rot_gp_button = QtGui.QPushButton(groupBox_2)
rot_gp_button.setObjectName(_fromUtf8("rot_gp_button"))
gridLayout_3.addWidget(rot_gp_button, 15, 3, 1, 1)
rot_gm_button = QtGui.QPushButton(groupBox_2)
rot_gm_button.setObjectName(_fromUtf8("rot_gm_button"))
gridLayout_3.addWidget(rot_gm_button, 15, 0, 1, 1)
diff3_entry = QtGui.QLineEdit(groupBox_2)
diff3_entry.setObjectName(_fromUtf8("diff3_entry"))
gridLayout_3.addWidget(diff3_entry, 3, 3, 1, 1)
tilt_label = QtGui.QLabel(groupBox_2)
tilt_label.setObjectName(_fromUtf8("tilt_label"))
gridLayout_3.addWidget(tilt_label, 4, 0, 1, 1)
rot_g_entry = QtGui.QLineEdit(groupBox_2)
rot_g_entry.setObjectName(_fromUtf8("rot_g_entry"))
gridLayout_3.addWidget(rot_g_entry, 15, 2, 1, 1)
uvw_button = QtGui.QCheckBox(groupBox_2)
uvw_button.setObjectName(_fromUtf8("uvw_button"))
gridLayout_3.addWidget(uvw_button, 12, 0, 1, 3)
rg_label = QtGui.QLabel(groupBox_2)
rg_label.setText(_fromUtf8(""))
rg_label.setObjectName(_fromUtf8("rg_label"))
gridLayout_3.addWidget(rg_label, 16, 2, 1, 1)
tilt_entry = QtGui.QLineEdit(groupBox_2)
tilt_entry.setObjectName(_fromUtf8("tilt_entry"))
gridLayout_3.addWidget(tilt_entry, 4, 3, 1, 1)
rot_diff_label = QtGui.QLabel(groupBox_2)
rot_diff_label.setObjectName(_fromUtf8("rot_diff_label"))
gridLayout_3.addWidget(rot_diff_label, 14, 0, 1, 4)
diff2_entry = QtGui.QLineEdit(groupBox_2)
diff2_entry.setObjectName(_fromUtf8("diff2_entry"))
gridLayout_3.addWidget(diff2_entry, 3, 2, 1, 1)
diff_label = QtGui.QLabel(groupBox_2)
diff_label.setObjectName(_fromUtf8("diff_label"))
gridLayout_3.addWidget(diff_label, 0, 0, 1, 3)
hexa_button = QtGui.QCheckBox(groupBox_2)
hexa_button.setObjectName(_fromUtf8("hexa_button"))
gridLayout_3.addWidget(hexa_button, 13, 0, 1, 3)
gridLayout.addWidget(groupBox_2, 2, 1, 1, 1)
groupBox = QtGui.QGroupBox(centralwidget)
groupBox.setMaximumSize(QtCore.QSize(150, 16777215))
groupBox.setObjectName(_fromUtf8("groupBox"))
gridLayout_2 = QtGui.QGridLayout(groupBox)
gridLayout_2.setMargin(5)
gridLayout_2.setSpacing(5)
gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
a_entry = QtGui.QLineEdit(groupBox)
a_entry.setObjectName(_fromUtf8("a_entry"))
gridLayout_2.addWidget(a_entry, 0, 1, 1, 1)
b_label = QtGui.QLabel(groupBox)
b_label.setObjectName(_fromUtf8("b_label"))
gridLayout_2.addWidget(b_label, 1, 0, 1, 1)
d_label = QtGui.QLabel(groupBox)
d_label.setObjectName(_fromUtf8("d_label"))
gridLayout_2.addWidget(d_label, 7, 0, 1, 1)
a_label = QtGui.QLabel(groupBox)
a_label.setObjectName(_fromUtf8("a_label"))
gridLayout_2.addWidget(a_label, 0, 0, 1, 1)
d_entry = QtGui.QLineEdit(groupBox)
d_entry.setObjectName(_fromUtf8("d_entry"))
gridLayout_2.addWidget(d_entry, 7, 2, 1, 1)
dp_button = QtGui.QPushButton(groupBox)
dp_button.setObjectName(_fromUtf8("dp_button"))
gridLayout_2.addWidget(dp_button, 7, 4, 1, 1)
c_label = QtGui.QLabel(groupBox)
c_label.setObjectName(_fromUtf8("c_label"))
gridLayout_2.addWidget(c_label, 2, 0, 1, 1)
dm_button = QtGui.QPushButton(groupBox)
dm_button.setObjectName(_fromUtf8("dm_button"))
gridLayout_2.addWidget(dm_button, 7, 1, 1, 1)
b_entry = QtGui.QLineEdit(groupBox)
b_entry.setObjectName(_fromUtf8("b_entry"))
gridLayout_2.addWidget(b_entry, 1, 1, 1, 1)
c_entry = QtGui.QLineEdit(groupBox)
c_entry.setObjectName(_fromUtf8("c_entry"))
gridLayout_2.addWidget(c_entry, 2, 1, 1, 1)
d_label_var = QtGui.QLabel(groupBox)
d_label_var.setText(_fromUtf8(""))
d_label_var.setObjectName(_fromUtf8("d_label_var"))
gridLayout_2.addWidget(d_label_var, 7, 5, 1, 1)
gamma_entry = QtGui.QLineEdit(groupBox)
gamma_entry.setObjectName(_fromUtf8("gamma_entry"))
gridLayout_2.addWidget(gamma_entry, 2, 5, 1, 1)
alpha_entry = QtGui.QLineEdit(groupBox)
alpha_entry.setObjectName(_fromUtf8("alpha_entry"))
gridLayout_2.addWidget(alpha_entry, 0, 5, 1, 1)
alpha_label = QtGui.QLabel(groupBox)
alpha_label.setObjectName(_fromUtf8("alpha_label"))
gridLayout_2.addWidget(alpha_label, 0, 2, 1, 3)
beta_label = QtGui.QLabel(groupBox)
beta_label.setObjectName(_fromUtf8("beta_label"))
gridLayout_2.addWidget(beta_label, 1, 2, 1, 2)
gamma_label = QtGui.QLabel(groupBox)
gamma_label.setObjectName(_fromUtf8("gamma_label"))
gridLayout_2.addWidget(gamma_label, 2, 2, 1, 3)
beta_entry = QtGui.QLineEdit(groupBox)
beta_entry.setObjectName(_fromUtf8("beta_entry"))
gridLayout_2.addWidget(beta_entry, 1, 5, 1, 1)
e_entry = QtGui.QLineEdit(groupBox)
e_entry.setObjectName(_fromUtf8("e_entry"))
gridLayout_2.addWidget(e_entry, 6, 5, 1, 1)
e_label = QtGui.QLabel(groupBox)
e_label.setObjectName(_fromUtf8("e_label"))
gridLayout_2.addWidget(e_label, 6, 0, 1, 5)
gridLayout.addWidget(groupBox, 1, 1, 1, 1)
groupBox_9 = QtGui.QGroupBox(centralwidget)
groupBox_9.setMaximumSize(QtCore.QSize(170, 16777215))
groupBox_9.setObjectName(_fromUtf8("groupBox_9"))
gridLayout_10 = QtGui.QGridLayout(groupBox_9)
gridLayout_10.setObjectName(_fromUtf8("gridLayout_10"))
size_var_label = QtGui.QLabel(groupBox_9)
size_var_label.setObjectName(_fromUtf8("size_var_label"))
gridLayout_10.addWidget(size_var_label, 4, 0, 1, 2)
color_trace_rouge = QtGui.QRadioButton(groupBox_9)
color_trace_rouge.setObjectName(_fromUtf8("color_trace_rouge"))
gridLayout_10.addWidget(color_trace_rouge, 2, 0, 1, 1)
size_var = QtGui.QLineEdit(groupBox_9)
size_var.setObjectName(_fromUtf8("size_var"))
gridLayout_10.addWidget(size_var, 4, 2, 1, 1)
color_trace_bleu = QtGui.QRadioButton(groupBox_9)
color_trace_bleu.setObjectName(_fromUtf8("color_trace_bleu"))
gridLayout_10.addWidget(color_trace_bleu, 0, 0, 1, 1)
color_trace_vert = QtGui.QRadioButton(groupBox_9)
color_trace_vert.setObjectName(_fromUtf8("color_trace_vert"))
gridLayout_10.addWidget(color_trace_vert, 0, 1, 1, 2)
style_box = QtGui.QCheckBox(groupBox_9)
style_box.setObjectName(_fromUtf8("style_box"))
gridLayout_10.addWidget(style_box, 3, 0, 1, 3)
gridLayout.addWidget(groupBox_9, 3, 1, 1, 1)
gridLayout_7.addLayout(gridLayout, 0, 1, 1, 1)
figure=plt.figure(facecolor='white',figsize=[2,2],dpi=100)
canvas = FigureCanvas(figure)
canvas.setMinimumSize(QtCore.QSize(800, 800))
gridLayout_7.addWidget(canvas, 0, 0, 1, 1)
MainWindow.setCentralWidget(centralwidget)
menubar = QtGui.QMenuBar(MainWindow)
#menubar.setGeometry(QtCore.QRect(0, 0, 1152, 25))
menubar.setObjectName(_fromUtf8("menubar"))
menuSave = QtGui.QMenu(menubar)
menuSave.setObjectName(_fromUtf8("menuSave"))
menuStructure = QtGui.QMenu(menubar)
menuStructure.setObjectName(_fromUtf8("menuStructure"))
menuAngle = QtGui.QMenu(menubar)
menuAngle.setObjectName(_fromUtf8("menuAngle"))
menuSchmid_factor = QtGui.QMenu(menubar)
menuSchmid_factor.setObjectName(_fromUtf8("menuSchmid_factor"))
menuResol = QtGui.QMenu(menubar)
menuResol.setObjectName(_fromUtf8("menuResol"))
MainWindow.setMenuBar(menubar)
statusbar = QtGui.QStatusBar(MainWindow)
statusbar.setObjectName(_fromUtf8("statusbar"))
MainWindow.setStatusBar(statusbar)
actionSave_figure = QtGui.QAction(MainWindow)
actionSave_figure.setObjectName(_fromUtf8("actionSave_figure"))
actionCalculate_Schmid_factor = QtGui.QAction(MainWindow)
actionCalculate_Schmid_factor.setObjectName(_fromUtf8("actionCalculate_Schmid_factor"))
actionCalculate_angle = QtGui.QAction(MainWindow)
actionCalculate_angle.setObjectName(_fromUtf8("actionCalculate_angle"))
menuSave.addAction(actionSave_figure)
menuAngle.addAction(actionCalculate_angle)
menuSchmid_factor.addAction(actionCalculate_Schmid_factor)
menubar.addAction(menuSave.menuAction())
menubar.addAction(menuStructure.menuAction())
menubar.addAction(menuAngle.menuAction())
menubar.addAction(menuSchmid_factor.menuAction())
menubar.addAction(menuResol.menuAction())

MainWindow.setWindowTitle(QtGui.QApplication.translate("Stereo-Proj", "Stereo-Proj", None, QtGui.QApplication.UnicodeUTF8))
groupBox_5.setTitle(QtGui.QApplication.translate("MainWindow", "Euler Angles", None, QtGui.QApplication.UnicodeUTF8))
lab_coord.setText(QtGui.QApplication.translate("MainWindow", "Tilt, Inclination", None, QtGui.QApplication.UnicodeUTF8))
button_trace2.setText(QtGui.QApplication.translate("MainWindow", "PLOT", None, QtGui.QApplication.UnicodeUTF8))
lab_euler2.setText(QtGui.QApplication.translate("MainWindow", "Phi 1, Phi, Phi2", None, QtGui.QApplication.UnicodeUTF8))
groupBox_3.setTitle(QtGui.QApplication.translate("MainWindow", "Rotation", None, QtGui.QApplication.UnicodeUTF8))
rx_buttonp.setText(QtGui.QApplication.translate("MainWindow", "+", None, QtGui.QApplication.UnicodeUTF8))
rx_label.setText(QtGui.QApplication.translate("MainWindow", "x", None, QtGui.QApplication.UnicodeUTF8))
rz_buttonm.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
rz_buttonp.setText(QtGui.QApplication.translate("MainWindow", "+", None, QtGui.QApplication.UnicodeUTF8))
ry_buttonp.setText(QtGui.QApplication.translate("MainWindow", "+", None, QtGui.QApplication.UnicodeUTF8))
rz_label.setText(QtGui.QApplication.translate("MainWindow", "z", None, QtGui.QApplication.UnicodeUTF8))
ry_label.setText(QtGui.QApplication.translate("MainWindow", "y", None, QtGui.QApplication.UnicodeUTF8))
rx_buttonm.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
ry_buttonm.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
groupBox_4.setTitle(QtGui.QApplication.translate("MainWindow", "Pole/Plane", None, QtGui.QApplication.UnicodeUTF8))
undo_trace_schmid.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
norm_button.setText(QtGui.QApplication.translate("MainWindow", "dhkl", None, QtGui.QApplication.UnicodeUTF8))
undo_sym_button.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
trace_schmid_button.setText(QtGui.QApplication.translate("MainWindow", "Schmid", None, QtGui.QApplication.UnicodeUTF8))
addpole_button.setText(QtGui.QApplication.translate("MainWindow", "Add", None, QtGui.QApplication.UnicodeUTF8))
trace_plan_button.setText(QtGui.QApplication.translate("MainWindow", " Plane", None, QtGui.QApplication.UnicodeUTF8))
undo_trace_plan_sym_button.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
trace_plan_sym_button.setText(QtGui.QApplication.translate("MainWindow", "Sym Plane", None, QtGui.QApplication.UnicodeUTF8))
sym_button.setText(QtGui.QApplication.translate("MainWindow", "Symmetry", None, QtGui.QApplication.UnicodeUTF8))
undo_trace_plan_button.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
undo_addpole_button.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
largeur_button.setText(QtGui.QApplication.translate("MainWindow", "Width", None, QtGui.QApplication.UnicodeUTF8))
groupBox_2.setTitle(QtGui.QApplication.translate("MainWindow", "Axis/Rotation", None, QtGui.QApplication.UnicodeUTF8))
button_trace.setText(QtGui.QApplication.translate("MainWindow", "PLOT", None, QtGui.QApplication.UnicodeUTF8))
inclination_label.setText(QtGui.QApplication.translate("MainWindow", "Inclination", None, QtGui.QApplication.UnicodeUTF8))
rot_gp_button.setText(QtGui.QApplication.translate("MainWindow", "+", None, QtGui.QApplication.UnicodeUTF8))
rot_gm_button.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
tilt_label.setText(QtGui.QApplication.translate("MainWindow", "Tilt", None, QtGui.QApplication.UnicodeUTF8))
uvw_button.setText(QtGui.QApplication.translate("MainWindow", "uvw", None, QtGui.QApplication.UnicodeUTF8))
rot_diff_label.setText(QtGui.QApplication.translate("MainWindow", "Rotation along g", None, QtGui.QApplication.UnicodeUTF8))
diff_label.setText(QtGui.QApplication.translate("MainWindow", "g-vector", None, QtGui.QApplication.UnicodeUTF8))
hexa_button.setText(QtGui.QApplication.translate("MainWindow", "hexa", None, QtGui.QApplication.UnicodeUTF8))
groupBox.setTitle(QtGui.QApplication.translate("MainWindow", "Crystal Parameters", None, QtGui.QApplication.UnicodeUTF8))
b_label.setText(QtGui.QApplication.translate("MainWindow", "b", None, QtGui.QApplication.UnicodeUTF8))
d_label.setText(QtGui.QApplication.translate("MainWindow", "d", None, QtGui.QApplication.UnicodeUTF8))
a_label.setText(QtGui.QApplication.translate("MainWindow", "a", None, QtGui.QApplication.UnicodeUTF8))
dp_button.setText(QtGui.QApplication.translate("MainWindow", "+", None, QtGui.QApplication.UnicodeUTF8))
c_label.setText(QtGui.QApplication.translate("MainWindow", "c", None, QtGui.QApplication.UnicodeUTF8))
dm_button.setText(QtGui.QApplication.translate("MainWindow", "-", None, QtGui.QApplication.UnicodeUTF8))
alpha_label.setText(QtGui.QApplication.translate("MainWindow", "alpha", None, QtGui.QApplication.UnicodeUTF8))
beta_label.setText(QtGui.QApplication.translate("MainWindow", "beta", None, QtGui.QApplication.UnicodeUTF8))
gamma_label.setText(QtGui.QApplication.translate("MainWindow", "gamma", None, QtGui.QApplication.UnicodeUTF8))
e_label.setText(QtGui.QApplication.translate("MainWindow", "max indices", None, QtGui.QApplication.UnicodeUTF8))
groupBox_9.setTitle(QtGui.QApplication.translate("MainWindow", "Layout", None, QtGui.QApplication.UnicodeUTF8))
size_var_label.setText(QtGui.QApplication.translate("MainWindow", "Marker size", None, QtGui.QApplication.UnicodeUTF8))
color_trace_rouge.setText(QtGui.QApplication.translate("MainWindow", "red", None, QtGui.QApplication.UnicodeUTF8))
color_trace_bleu.setText(QtGui.QApplication.translate("MainWindow", "blue", None, QtGui.QApplication.UnicodeUTF8))
color_trace_vert.setText(QtGui.QApplication.translate("MainWindow", "green", None, QtGui.QApplication.UnicodeUTF8))
style_box.setText(QtGui.QApplication.translate("MainWindow", "open/filled", None, QtGui.QApplication.UnicodeUTF8))
menuSave.setTitle(QtGui.QApplication.translate("MainWindow", "Save", None, QtGui.QApplication.UnicodeUTF8))
menuStructure.setTitle(QtGui.QApplication.translate("MainWindow", "Structure", None, QtGui.QApplication.UnicodeUTF8))
menuAngle.setTitle(QtGui.QApplication.translate("MainWindow", "Angle", None, QtGui.QApplication.UnicodeUTF8))
menuSchmid_factor.setTitle(QtGui.QApplication.translate("MainWindow", "Schmid factor", None, QtGui.QApplication.UnicodeUTF8))
menuResol.setTitle(QtGui.QApplication.translate("Index", "Resolution", None, QtGui.QApplication.UnicodeUTF8))
actionSave_figure.setText(QtGui.QApplication.translate("MainWindow", "Save figure", None, QtGui.QApplication.UnicodeUTF8))
actionCalculate_Schmid_factor.setText(QtGui.QApplication.translate("MainWindow", "calculate Schmid factor", None, QtGui.QApplication.UnicodeUTF8))
actionCalculate_angle.setText(QtGui.QApplication.translate("MainWindow", "Calculate angle", None, QtGui.QApplication.UnicodeUTF8))


MainWindow.setTabOrder( a_entry,  b_entry)
MainWindow.setTabOrder( b_entry,  c_entry)
MainWindow.setTabOrder( c_entry,  alpha_entry)
MainWindow.setTabOrder( alpha_entry,  beta_entry)
MainWindow.setTabOrder( beta_entry,  gamma_entry)
MainWindow.setTabOrder( gamma_entry,  e_entry)
MainWindow.setTabOrder( e_entry,  rx_buttonm)
MainWindow.setTabOrder( rx_buttonm,  rx_entry)
MainWindow.setTabOrder( rx_entry,  rx_buttonp)
MainWindow.setTabOrder( rx_buttonp,  ry_buttonm)
MainWindow.setTabOrder( ry_buttonm,  ry_entry)
MainWindow.setTabOrder( ry_entry,  ry_buttonp)
MainWindow.setTabOrder( ry_buttonp,  rz_buttonm)
MainWindow.setTabOrder( rz_buttonm,  rz_entry)
MainWindow.setTabOrder( rz_entry,  rz_buttonp)
MainWindow.setTabOrder( rz_buttonp,  d_entry)
MainWindow.setTabOrder( d_entry,  dm_button)
MainWindow.setTabOrder( dm_button,  dp_button)
MainWindow.setTabOrder( dp_button,  diff1_entry)
MainWindow.setTabOrder( diff1_entry,  diff2_entry)
MainWindow.setTabOrder( diff2_entry,  diff3_entry)
MainWindow.setTabOrder( diff3_entry,  tilt_entry)
MainWindow.setTabOrder( tilt_entry,  inclinaison_entry)
MainWindow.setTabOrder( inclinaison_entry,  button_trace)
MainWindow.setTabOrder( button_trace,  uvw_button)
MainWindow.setTabOrder( uvw_button,  hexa_button)
MainWindow.setTabOrder( hexa_button,  rot_gm_button)
MainWindow.setTabOrder( rot_gm_button,  rot_g_entry)
MainWindow.setTabOrder( rot_g_entry,  rot_gp_button)
MainWindow.setTabOrder( rot_gp_button,  pole1_entry)
MainWindow.setTabOrder( pole1_entry,  pole2_entry)
MainWindow.setTabOrder( pole2_entry,  pole3_entry)
MainWindow.setTabOrder( pole3_entry,  addpole_button)
MainWindow.setTabOrder( addpole_button,  undo_addpole_button)
MainWindow.setTabOrder( undo_addpole_button,  sym_button)
MainWindow.setTabOrder( sym_button,  undo_sym_button)
MainWindow.setTabOrder( undo_sym_button,  trace_plan_button)
MainWindow.setTabOrder( trace_plan_button,  undo_trace_plan_button)
MainWindow.setTabOrder( undo_trace_plan_button,  trace_plan_sym_button)
MainWindow.setTabOrder( trace_plan_sym_button,  undo_trace_plan_sym_button)
MainWindow.setTabOrder( undo_trace_plan_sym_button,  trace_schmid_button)
MainWindow.setTabOrder( trace_schmid_button,  undo_trace_schmid)
MainWindow.setTabOrder( undo_trace_schmid,  largeur_button)
MainWindow.setTabOrder( largeur_button,  norm_button)
MainWindow.setTabOrder( norm_button,  phi1_entry)
MainWindow.setTabOrder( phi1_entry,  phi_entry)
MainWindow.setTabOrder( phi_entry,  phi2_entry)
MainWindow.setTabOrder( phi2_entry,  button_trace2)
MainWindow.setTabOrder( button_trace2,  color_trace_bleu)
MainWindow.setTabOrder( color_trace_bleu,  color_trace_vert)
MainWindow.setTabOrder( color_trace_vert,  color_trace_rouge)
MainWindow.setTabOrder( color_trace_rouge,  style_box)
MainWindow.setTabOrder( style_box,  size_var)



def structure(item):
    global x0, var_hexa, d_label_var, e_entry
    
    a_entry.setText(str(item[1]))
    b_entry.setText(str(item[2]))
    c_entry.setText(str(item[3]))
    alpha_entry.setText(str(item[4]))
    beta_entry.setText(str(item[5]))
    gamma_entry.setText(str(item[6]))
    if eval(item[4])==90 and eval(item[5])==90 and eval(item[6])==120 :
        hexa_button.setChecked(True)
        e_entry.setText('2')
        d_label_var.setText('3')
    else:
        d_entry.setText('1')
        e_entry.setText('1')
    
    
    
file_struct=open(os.path.join(os.path.dirname(__file__), 'structure.txt') ,"r")

x0=[]

for line in file_struct:
    x0.append(map(str, line.split()))

i=0
file_struct.close()            

for item in x0:
    entry = menuStructure.addAction(item[0])
    MainWindow.connect(entry,QtCore.SIGNAL('triggered()'), lambda item=item: structure(item))
    i=i+1


####################################################################
##### Fonction angle entre deux directions
####################################################################


class AngleCalc(QtGui.QDialog):
    global Dstar
    def __init__(self, parent=None):
        
        super(AngleCalc, self).__init__(parent)
        
        self.gridLayout = QtGui.QGridLayout(self)
        self.gridLayout.setMargin(10)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.c10 = QtGui.QLineEdit(self)
        self.c10.setObjectName(_fromUtf8("c10"))
        self.gridLayout.addWidget(self.c10, 1, 0, 1, 1)
        self.c11 = QtGui.QLineEdit(self)
        self.c11.setObjectName(_fromUtf8("c11"))
        self.gridLayout.addWidget(self.c11, 2, 0, 1, 1)
        self.c12 = QtGui.QLineEdit(self)
        self.c12.setObjectName(_fromUtf8("c12"))
        self.gridLayout.addWidget(self.c12, 3, 0, 1, 1)
        self.c20 = QtGui.QLineEdit(self)
        self.c20.setObjectName(_fromUtf8("c20"))
        self.gridLayout.addWidget(self.c20, 1, 1, 1, 1)
        self.c21 = QtGui.QLineEdit(self)
        self.c21.setObjectName(_fromUtf8("c21"))
        self.gridLayout.addWidget(self.c21, 2, 1, 1, 1)
        self.c22 = QtGui.QLineEdit(self)
        self.c22.setObjectName(_fromUtf8("c22"))
        self.gridLayout.addWidget(self.c22, 3, 1, 1, 1)
        
        
        self.label = QtGui.QLabel(self)
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 4, 0, 1, 2)
        
        
        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout.addWidget(self.buttonBox, 5, 0, 1, 2)
        self.label_2 = QtGui.QLabel(self)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_2.setText("d1")
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.label_3 = QtGui.QLabel(self)
        self.label_3.setText("d2")
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 0, 1, 1, 1)
        
        


        self.buttonBox.rejected.connect(self.close)
        self.buttonBox.accepted.connect(self.angle)
        
    def angle(self):
       global Dstar
       c100=np.float(self.c10.text())
       c110=np.float(self.c11.text())
       c120=np.float(self.c12.text())
       c200=np.float(self.c20.text())
       c210=np.float(self.c21.text())
       c220=np.float(self.c22.text())
       c1=np.array([c100,c110,c120])
       c2=np.array([c200,c210,c220])
       if uvw_button.isChecked==True: 
           c1c=np.dot(Dstar,c1)
           c2c=np.dot(Dstar,c2)
       else:
           c1c=np.dot(Dstar,c1)
           c2c=np.dot(Dstar,c2)
       the=np.arccos(np.dot(c1c,c2c)/(np.linalg.norm(c1c)*np.linalg.norm(c2c)))                   
       thes=str(np.around(the*180/pi,decimals=2))        
       self.label.setText(thes)
                        
        
dialogAngleCalc = AngleCalc()
dialogAngleCalc.setWindowTitle("Calculate Angle")
MainWindow.connect(actionCalculate_angle, QtCore.SIGNAL('triggered()'), dialogAngleCalc.exec_)  
 
##################################################

##### Schmid factor calculation

################################################### 

class SchmidCalc(QtGui.QDialog):
    global Dstar
    def __init__(self, parent=None):
        
        super(SchmidCalc, self).__init__(parent)
        
        
        self.gridLayout = QtGui.QGridLayout(self)
        self.gridLayout.setMargin(10)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.c10 = QtGui.QLineEdit(self)
        self.c10.setObjectName(_fromUtf8("c10"))
        self.gridLayout.addWidget(self.c10, 1, 0, 1, 1)
        self.c11 = QtGui.QLineEdit(self)
        self.c11.setObjectName(_fromUtf8("c11"))
        self.gridLayout.addWidget(self.c11, 2, 0, 1, 1)
        self.c12 = QtGui.QLineEdit(self)
        self.c12.setObjectName(_fromUtf8("c12"))
        self.gridLayout.addWidget(self.c12, 3, 0, 1, 1)
        self.c20 = QtGui.QLineEdit(self)
        self.c20.setObjectName(_fromUtf8("c20"))
        self.gridLayout.addWidget(self.c20, 1, 1, 1, 1)
        self.c21 = QtGui.QLineEdit(self)
        self.c21.setObjectName(_fromUtf8("c21"))
        self.gridLayout.addWidget(self.c21, 2, 1, 1, 1)
        self.c22 = QtGui.QLineEdit(self)
        self.c22.setObjectName(_fromUtf8("c22"))
        self.gridLayout.addWidget(self.c22, 3, 1, 1, 1)
        
        
        self.label = QtGui.QLabel(self)
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 4, 0, 1, 2)
        
        
        self.buttonBox = QtGui.QDialogButtonBox(self)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout.addWidget(self.buttonBox, 5, 0, 1, 2)
        self.label_2 = QtGui.QLabel(self)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_2.setText("b")
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.label_3 = QtGui.QLabel(self)
        self.label_3.setText("n")
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 0, 1, 1, 1)
        self.buttonBox.rejected.connect(self.close)
        self.buttonBox.accepted.connect(self.schmid)
        
    def schmid(self):
       global D, Dstar,M
       c100=np.float(self.c10.text())
       c110=np.float(self.c11.text())
       c120=np.float(self.c12.text())
       c200=np.float(self.c20.text())
       c210=np.float(self.c21.text())
       c220=np.float(self.c22.text())
       n=np.array([c100,c110,c120])
       b=np.array([c200,c210,c220])
       npr=np.dot(Dstar,n)
       bpr=np.dot(Dstar,b)
       npr2=np.dot(M,npr)
       bpr2=np.dot(M,bpr)
       T=np.array([0,1,0])
       anglen=np.arccos(np.dot(npr2,T)/np.linalg.norm(npr2))
       angleb=np.arccos(np.dot(bpr2,T)/np.linalg.norm(bpr2))
       s=np.cos(anglen)*np.cos(angleb)
       s2=str(np.around(s,decimals=2))
       self.label.setText(s2)

                        
 
dialogSchmidCalc = SchmidCalc()
dialogSchmidCalc.setWindowTitle("Calculate Schmid Factor")
MainWindow.connect(actionCalculate_Schmid_factor, QtCore.SIGNAL('triggered()'), dialogSchmidCalc.exec_) 



def image_save():
    filename=QtGui.QFileDialog.getSaveFileName( MainWindow,"Save file", "", ".png")
    pixmap = QtGui.QPixmap.grabWidget(canvas)
    pixmap.save(str(filename)+".png")
     


MainWindow.connect(actionSave_figure, QtCore.SIGNAL('triggered()'), image_save) 

              
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
    lat=90+lat
    longi=-longi+180
    if longi>90:
     longi=longi-180
    c=str(np.around(longi,decimals=1))+str(',')+str(np.around(lat,decimals=1))
    coord_label.setText(str(c))


#########################################################################


def resol(item):
        
    #print(MainWindow.size(), np.int(item[0]))
    canvas.setMinimumSize(QtCore.QSize(np.int(item[0]), np.int(item[2])))
    #canvas.resize(np.int(item[0]), np.int(item[2]))    
    MainWindow.resize(0,0)
    print(MainWindow.size())
    

file_resol=open(os.path.join(os.path.dirname(__file__), 'resolution.txt') ,"r")

xr=[]

for line in file_resol:
    xr.append(map(str, line.split()))

file_resol.close()            

for itemr in xr:
    entry = menuResol.addAction(itemr[0]+itemr[1]+itemr[2])
    MainWindow.connect(entry,QtCore.SIGNAL('triggered()'), lambda itemr=itemr: resol(itemr))    
 #######################################################################


fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
img=Image.open(fn)
a = figure.add_subplot(111)
a.axis('off')
a.imshow(img)
a.figure.canvas.draw()

MainWindow.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding) 


button_trace2.clicked.connect(princ2)
button_trace.clicked.connect(princ)

rx_buttonp.clicked.connect(rxp)
rx_buttonm.clicked.connect(rxm)
ry_buttonp.clicked.connect(ryp)
ry_buttonm.clicked.connect(rym)
rz_buttonp.clicked.connect(rzp)
rz_buttonm.clicked.connect(rzm)

rot_gm_button.clicked.connect(rotgm)
rot_gp_button.clicked.connect(rotgp)

addpole_button.clicked.connect(addpole)
undo_addpole_button.clicked.connect(undo_addpole)
sym_button.clicked.connect(addpole_sym)
undo_sym_button.clicked.connect(undo_sym)
trace_plan_button.clicked.connect(trace_addplan)
undo_trace_plan_button.clicked.connect(undo_trace_addplan)
trace_plan_sym_button.clicked.connect(trace_plan_sym)
undo_trace_plan_sym_button.clicked.connect(undo_trace_plan_sym)
trace_schmid_button.clicked.connect(schmid_trace)
undo_trace_schmid.clicked.connect(undo_schmid_trace)
largeur_plan=Window()
largeur_plan.setWindowTitle("Apparent Width")
MainWindow.connect(largeur_button, QtCore.SIGNAL('clicked()'), largeur_plan.exec_) 
norm_button.clicked.connect(dhkl)
dm_button.clicked.connect(dm)
dp_button.clicked.connect(dp)

figure.canvas.mpl_connect('motion_notify_event', coordinates)
figure.canvas.mpl_connect('button_press_event', click_a_pole)


dmip=0
tr_schmid=np.zeros((1,3))
x=0
y=0
z=0
g=0

color_trace_bleu.setChecked(True)

d_label_var.setText('0')
alpha_entry.setText('90')
beta_entry.setText('90')
gamma_entry.setText('90')

a_entry.setText('1')
b_entry.setText('1')
c_entry.setText('1')

phi1_entry.setText('0')
phi_entry.setText('0')
phi2_entry.setText('0')

e_entry.setText('1')

rx_label.setText('0.0')
ry_label.setText('0.0')
rz_label.setText('0.0')
rg_label.setText('0.0')
angle_euler_label.setText(' ')
size_var.setText('40')

e_entry.setText('1')
rx_entry.setText('5')
ry_entry.setText('5')
rz_entry.setText('5')

d_entry.setText('1')
rot_g_entry.setText('5')



MainWindow.show()
app.exec_()



