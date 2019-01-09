#!/usr/bin/python

######################################################################
#
#
# Stereo-Proj is a python utility to plot stereographic projetion of a given crystal. It is designed
# to be used in electron microscopy experiments.
# Authors: F. Mompiou, CEMES-CNRS
#
#######################################################################


from __future__ import division
import numpy as np
from PyQt4 import QtGui, QtCore
import sys
import random
import os
from PIL  import Image
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import pyplot as plt
import matplotlib as mpl
import stereoprojUI

                 
#font size on plot 
mpl.rcParams['font.size'] = 12

################
#       Misc
################

def unique_rows(a):
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    return np.unique(b).view(a.dtype).reshape(-1, a.shape[1])



###################################################################"
##### Projection
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
##### Rotation Euler
####################################################################

def rotation(phi1,phi,phi2):
   phi1=phi1*np.pi/180;
   phi=phi*np.pi/180;
   phi2=phi2*np.pi/180;
   R=np.array([[np.cos(phi1)*np.cos(phi2)-np.cos(phi)*np.sin(phi1)*np.sin(phi2),
            -np.cos(phi)*np.cos(phi2)*np.sin(phi1)-np.cos(phi1)*
            np.sin(phi2),np.sin(phi)*np.sin(phi1)],[np.cos(phi2)*np.sin(phi1)
            +np.cos(phi)*np.cos(phi1)*np.sin(phi2),np.cos(phi)*np.cos(phi1)
            *np.cos(phi2)-np.sin(phi1)*np.sin(phi2), -np.cos(phi1)*np.sin(phi)],
            [np.sin(phi)*np.sin(phi2), np.cos(phi2)*np.sin(phi), np.cos(phi)]],float)
   return R

####################################################################
##### Rotation around a given axis
####################################################################

def Rot(th,a,b,c):
   th=th*np.pi/180;
   aa=a/np.linalg.norm([a,b,c]);
   bb=b/np.linalg.norm([a,b,c]);
   cc=c/np.linalg.norm([a,b,c]);
   c1=np.array([[1,0,0],[0,1,0],[0,0,1]],float)
   c2=np.array([[aa**2,aa*bb,aa*cc],[bb*aa,bb**2,bb*cc],[cc*aa,
                cc*bb,cc**2]],float)
   c3=np.array([[0,-cc,bb],[cc,0,-aa],[-bb,aa,0]],float)
   R=np.cos(th)*c1+(1-np.cos(th))*c2+np.sin(th)*c3

   return R    

#######################
#
# Layout functions
#
#######################

def color_trace():
        color_trace=1
        if ui.color_trace_bleu.isChecked():
                color_trace=1
        if ui.color_trace_bleu.isChecked():
                color_trace=2
        if ui.color_trace_rouge.isChecked():
                color_trace=3
        return color_trace

def var_uvw():
        var_uvw=0
        if ui.uvw_button.isChecked():
                var_uvw=1
        
        return var_uvw

def var_hexa():
        var_hexa=0
        if ui.hexa_button.isChecked():
                var_hexa=1
        
        return var_hexa

def var_carre():
        var_carre=0
        if ui.style_box.isChecked():
                var_carre=1
        
        return var_carre

####################################################################
#
#  Crystal definition
#
####################################################################

def crist():
    global axes,axesh,D,Dstar,V
    a=np.float(ui.a_entry.text())
    b=np.float(ui.b_entry.text())
    c=np.float(ui.c_entry.text())
    alpha=np.float(ui.alpha_entry.text())
    beta=np.float(ui.beta_entry.text())
    gamma=np.float(ui.gamma_entry.text())
    e=np.int(ui.e_entry.text())
    d2=np.float(ui.d_label_var.text())
    alpha=alpha*np.pi/180;
    beta=beta*np.pi/180;
    gamma=gamma*np.pi/180;
    V=a*b*c*np.sqrt(1-(np.cos(alpha)**2)-(np.cos(beta))**2-(np.cos(gamma))**2+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
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

######################################################
#
# Reduce number of poles/directions as a function of d-spacing (plus or minus)
#
#######################################################

def dm():
    global dmip,a,minx,maxx,miny,maxy
    
    dmip=dmip-np.float(ui.d_entry.text())
    ui.d_label_var.setText(str(dmip))
    crist()
    trace()
    
    return dmip
    
def dp():
    global dmip, a

    dmip=dmip+np.float(ui.d_entry.text())
    ui.d_label_var.setText(str(dmip))
    crist()    
    trace()
    
    return dmip 
    

#############################################################################
#
#  Apparent width class for dialog box: plot the width of a plane of given normal hkl with the tilt alpha angle
#
##############################################################################

class LargeurPlanCalc(QtGui.QDialog):
    def __init__(self, parent=None):
        super(LargeurPlanCalc, self).__init__(parent)
 
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
 
    def plot(self):
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

        ax = self.figure.add_subplot(111)
        ax.hold(False)
        ax.plot(range(-40,41,2),la[0,:])
        self.canvas.draw()

    
####################################################################
#
# Schmid factor
#
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
#
#  Plot iso-schmid factor, ie for a given plan the locus of b with a given schmid factor (Oy direction
# assumed to be the straining axis
#
####################################################################

def schmid_trace():
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid,a,minx,maxx,miny,maxy
    
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
        
    tr_schmid=np.vstack((tr_schmid,np.array([pole1,pole2,pole3])))
    trace()
    
def undo_schmid_trace():
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid,a,minx,maxx,miny,maxy
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
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
        t=np.linspace(0,2*np.pi,100)
        r,t=np.meshgrid(r,t)
        F=fact(angleb,r,t,n)
        lev=[-0.5,-0.4,-0.3,-0.2,0.2,0.3,0.4,0.5]
        CS=a.contour(r*np.cos(t)+300, r*np.sin(t)+300, F,lev,linewidths=2)
        fmt = {}
        strs = [ '('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') -0.5','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') -0.4','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') -0.3','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') -0.2','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') 0.2','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') 0.3','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') 0.4','('+str(np.int(b1))+str(np.int(b2))+str(np.int(b3))+') 0.5']
        for l,s in zip( CS.levels, strs ):
			fmt[l] = s
	a.clabel(CS,fmt=fmt,fontsize=10,inline=True)						

    
###########################################################################
#
# Rotation of the sample. If Lock Axes is off rotation are along y,x,z directions. If not, the y and z axes 
# of the sample are locked to the crystal axes when the check box is ticked. It mimics double-tilt holder (rotation of alpha along fixed x and rotation of beta along the beta tilt moving axis) or  tilt-rotation holder  (rotation of alpha along fixed # x and rotation of z along the z-rotation moving axis).
#
###########################################################################"

def lock():
	global M, var_lock,M_lock
        
        if ui.lock_checkButton.isChecked():
                var_lock=1
                M_lock=M
        else:
        	var_lock,M_lock=0,0

        return var_lock,M_lock


def rot_alpha_p():
    global angle_alpha,M,a,trP

    tha=np.float(ui.angle_alpha_entry.text())
    M=np.dot(Rot(tha,0,1,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    ui.angle_euler_label.setText(t)
    angle_alpha=angle_alpha+np.float(ui.angle_alpha_entry.text())
    ui.angle_alpha_label_2.setText(str(angle_alpha))
    return angle_alpha,M
    
    
def rot_alpha_m():
    global angle_alpha,M,a,trP

    tha=-np.float(ui.angle_alpha_entry.text())
    M=np.dot(Rot(tha,0,1,0),M)
    trace()
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    ui.angle_euler_label.setText(t)
    angle_alpha=angle_alpha-np.float(ui.angle_alpha_entry.text())
    ui.angle_alpha_label_2.setText(str(angle_alpha))
    return angle_alpha,M

    
def rot_beta_m():
    global angle_beta,M,angle_alpha, angle_z, var_lock, M_lock
   
    if var_lock==0:
    	AxeY=np.array([1,0,0])
    else:
   	A=np.dot(np.linalg.inv(M_lock),np.array([1,0,0]))
	C=np.dot(np.linalg.inv(Dstar),A)
	AxeY=C/np.linalg.norm(C)
    	AxeY=np.dot(M,AxeY)
    
    thb=-np.float(ui.angle_beta_entry.text())
    M=np.dot(Rot(thb,AxeY[0],AxeY[1],AxeY[2]),M)
    trace()
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    ui.angle_euler_label.setText(t)
    angle_beta=angle_beta-np.float(ui.angle_beta_entry.text())
    ui.angle_beta_label_2.setText(str(angle_beta))
    return angle_beta,M   
   
def rot_beta_p():
    global angle_beta,M,angle_alpha, angle_z, var_lock, M_lock
   
    if var_lock==0:
    	AxeY=np.array([1,0,0])
    else:
   	A=np.dot(np.linalg.inv(M_lock),np.array([1,0,0]))
	C=np.dot(np.linalg.inv(Dstar),A)
	AxeY=C/np.linalg.norm(C)
    	AxeY=np.dot(M,AxeY)
    
    thb=np.float(ui.angle_beta_entry.text())
    M=np.dot(Rot(thb,AxeY[0],AxeY[1],AxeY[2]),M)
    trace()
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    ui.angle_euler_label.setText(t)
    angle_beta=angle_beta+np.float(ui.angle_beta_entry.text())
    ui.angle_beta_label_2.setText(str(angle_beta))
    return angle_beta,M   

def rot_z_m():
    global angle_beta,M,angle_alpha, angle_z, var_lock, M_lock
   
    if var_lock==0:
    	AxeZ=np.array([0,0,1])
    else:
   	A=np.dot(np.linalg.inv(M_lock),np.array([0,0,1]))
	C=np.dot(np.linalg.inv(Dstar),A)
	AxeZ=C/np.linalg.norm(C)
    	AxeZ=np.dot(M,AxeZ)
    
    thz=-np.float(ui.angle_z_entry.text())
    M=np.dot(Rot(thz,AxeZ[0],AxeZ[1],AxeZ[2]),M)
    trace()
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    ui.angle_euler_label.setText(t)
    angle_z=angle_z-np.float(ui.angle_z_entry.text())
    ui.angle_z_label_2.setText(str(angle_z))
    return angle_z,M      
   
def rot_z_p():
    global angle_beta,M,angle_alpha, angle_z, var_lock, M_lock
   
    if var_lock==0:
    	AxeZ=np.array([0,0,1])
    else:
   	A=np.dot(np.linalg.inv(M_lock),np.array([0,0,1]))
	C=np.dot(np.linalg.inv(Dstar),A)
	AxeZ=C/np.linalg.norm(C)
    	AxeZ=np.dot(M,AxeZ)
    
    thz=np.float(ui.angle_z_entry.text())
    M=np.dot(Rot(thz,AxeZ[0],AxeZ[1],AxeZ[2]),M)
    trace()
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    ui.angle_euler_label.setText(t)
    angle_z=angle_z+np.float(ui.angle_z_entry.text())
    ui.angle_z_label_2.setText(str(angle_z))
    return angle_z,M      

####################################################################
#
# Rotate around a given pole
#
####################################################################


def rotgm():
    global g,M,Dstar,a

    thg=-np.float(ui.rot_g_entry.text())
    diff1=np.float(ui.diff1_entry.text())
    diff2=np.float(ui.diff2_entry.text())
    diff3=np.float(ui.diff3_entry.text())
    A=np.array([diff1,diff2,diff3])
    Ad=np.dot(Dstar,A)    
    Ap=np.dot(M,Ad)/np.linalg.norm(np.dot(M,Ad))
    M=np.dot(Rot(thg,Ap[0],Ap[1],Ap[2]),M)
    trace()    
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    ui.angle_euler_label.setText(t)
    g=g-np.float(ui.rot_g_entry.text())
    ui.rg_label.setText(str(g))
    return g,M
    
def rotgp():
    global g,M,D

    thg=np.float(ui.rot_g_entry.text())
    diff1=np.float(ui.diff1_entry.text())
    diff2=np.float(ui.diff2_entry.text())
    diff3=np.float(ui.diff3_entry.text())
    A=np.array([diff1,diff2,diff3])
    Ad=np.dot(Dstar,A)    
    Ap=np.dot(M,Ad)/np.linalg.norm(np.dot(M,Ad))
    M=np.dot(Rot(thg,Ap[0],Ap[1],Ap[2]),M)
    trace()    
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    ui.angle_euler_label.setText(t)
    g=g+np.float(ui.rot_g_entry.text())
    ui.rg_label.setText(str(g))
    return g,M

###################################################
#
# Mirror the sample
#
#############################################

#def mirror():
#    global M,a,trP
##    a = figure.add_subplot(111)     
##    a.figure.clear()
#    
#    M_r=np.array([[1,0,0],[0,1,0],[0,0,-1]])
#    M=np.dot(M_r,M)
#    trace()
#    phir=np.arccos(M[2,2])*180/np.pi
#    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
#    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
#    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
#    
#    ui.angle_euler_label.setText(t)
#    
#    return M


####################################################################
#
# Add a given pole and equivalent ones
#
####################################################################

def pole(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar
    
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
    a=np.float(ui.a_entry.text())
    b=np.float(ui.b_entry.text())
    c=np.float(ui.c_entry.text())
    alpha=np.float(ui.alpha_entry.text())
    beta=np.float(ui.beta_entry.text())
    gamma=np.float(ui.gamma_entry.text())
    alpha=alpha*np.pi/180;
    beta=beta*np.pi/180;
    gamma=gamma*np.pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])    
    ds=(np.sqrt(np.dot(np.array([pole1,pole2,pole3]),np.dot(np.linalg.inv(G),np.array([pole1,pole2,pole3])))))
    return ds
    
def addpole_sym():
    global M,axes,axesh,T,V,D,Dstar,G    
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
    a=np.float(ui.a_entry.text())
    b=np.float(ui.b_entry.text())
    c=np.float(ui.c_entry.text())
    alpha=np.float(ui.alpha_entry.text())
    beta=np.float(ui.beta_entry.text())
    gamma=np.float(ui.gamma_entry.text())
    alpha=alpha*np.pi/180;
    beta=beta*np.pi/180;
    gamma=gamma*np.pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])     
    v=d(pole1,pole2,pole3)
    
    pole(pole1,pole2,pole3)
    if np.abs(alpha-np.pi/2)<0.001 and np.abs(beta-np.pi/2)<0.001 and np.abs(gamma-2*np.pi/3)<0.001:
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
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
    a=np.float(ui.a_entry.text())
    b=np.float(ui.b_entry.text())
    c=np.float(ui.c_entry.text())
    alpha=np.float(ui.alpha_entry.text())
    beta=np.float(ui.beta_entry.text())
    gamma=np.float(ui.gamma_entry.text())
    alpha=alpha*np.pi/180;
    beta=beta*np.pi/180;
    gamma=gamma*np.pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])     
    v=d(pole1,pole2,pole3)
    
    undo_pole(pole1,pole2,pole3)
    if np.abs(alpha-np.pi/2)<0.001 and np.abs(beta-np.pi/2)<0.001 and np.abs(gamma-2*np.pi/3)<0.001:
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
    
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
    pole(pole1,pole2,pole3)
    trace()
    
def undo_addpole():
    global M,axes,axesh,T,V,D,Dstar,trP,tr_schmid,nn
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
    undo_pole(pole1,pole2,pole3)
    
    trace()
    
####################################################################
#
# Plot a given plane and equivalent ones
#
####################################################################

def trace_plan(pole1,pole2,pole3):
    global M,axes,axesh,T,V,D,Dstar,trP
    
    pole_i=0
    pole_c=color_trace()
    
    if var_hexa()==1:
        if var_uvw()==1:
            pole1=2*np.float(ui.pole1_entry.text())+np.float(ui.pole2_entry.text())
            pole2=2*np.float(ui.pole2_entry.text())+np.float(ui.pole1_entry.text())
            pole3=np.float(ui.pole3_entry.text())
            pole_i=1
    
        
    trP=np.vstack((trP,np.array([pole1,pole2,pole3,pole_i,pole_c])))
    b=np.ascontiguousarray(trP).view(np.dtype((np.void, trP.dtype.itemsize * trP.shape[1])))
    
    trP=np.unique(b).view(trP.dtype).reshape(-1, trP.shape[1])
    #print(trP)
    
    
def trace_addplan():
    global M,axes,axesh,T,V,D,Dstar,trP
    
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
    
    trace_plan(pole1,pole2,pole3)
    trace()
    
def undo_trace_addplan():
    global M,axes,axesh,T,V,D,Dstar,trP
    
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
    
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
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
    a=np.float(ui.a_entry.text())
    b=np.float(ui.b_entry.text())
    c=np.float(ui.c_entry.text())
    alpha=np.float(ui.alpha_entry.text())
    beta=np.float(ui.beta_entry.text())
    gamma=np.float(ui.gamma_entry.text())
    alpha=alpha*np.pi/180;
    beta=beta*np.pi/180;
    gamma=gamma*np.pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])     
    v=d(pole1,pole2,pole3)
    
    trace_plan(pole1,pole2,pole3)
    if np.abs(alpha-np.pi/2)<0.001 and np.abs(beta-np.pi/2)<0.001 and np.abs(gamma-2*np.pi/3)<0.001:
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
    pole1=np.float(ui.pole1_entry.text())
    pole2=np.float(ui.pole2_entry.text())
    pole3=np.float(ui.pole3_entry.text())
    a=np.float(ui.a_entry.text())
    b=np.float(ui.b_entry.text())
    c=np.float(ui.c_entry.text())
    alpha=np.float(ui.alpha_entry.text())
    beta=np.float(ui.beta_entry.text())
    gamma=np.float(ui.gamma_entry.text())
    alpha=alpha*np.pi/180;
    beta=beta*np.pi/180;
    gamma=gamma*np.pi/180;
    G=np.array([[a**2,a*b*np.cos(gamma),a*c*np.cos(beta)],[a*b*np.cos(gamma),b**2,b*c*np.cos(alpha)],[a*c*np.cos(beta),b*c*np.cos(alpha),c**2]])     
    v=d(pole1,pole2,pole3)
    
    undo_trace_plan(pole1,pole2,pole3)
    if np.abs(alpha-np.pi/2)<0.001 and np.abs(beta-np.pi/2)<0.001 and np.abs(gamma-2*np.pi/3)<0.001:
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
             t=np.arctan2(S[1],S[0])*180/np.pi
        w=0
        ph=np.arccos(S[2]/r)*180/np.pi
        for g in np.linspace(-np.pi,np.pi-0.00001,100):
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
#
# Click a pole 
#
####################################################################    

def click_a_pole(event):
        
    global M,Dstar,D,minx,maxx,miny,maxy,a
      
    if event.button==3:
            x=event.xdata
            y=event.ydata
             
            x=(x-300)/300
            y=(y-300)/300
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
#
# Inclinaison-beta indicator when the mouse is on the stereo
#
#################################################################### 
  
def coordinates(event):

    if event.xdata and event.ydata:
	    x=event.xdata
	    y=event.ydata    
	    x=(x-300)/300
	    y=(y-300)/300
	    long0=90*np.pi/180
	    lat0=0
	    r=np.sqrt(x**2+y**2)
	    c=2*np.arctan(r)
	    longi=(long0+np.arctan2(x*np.sin(c),r*np.cos(lat0)*np.cos(c)-y*np.sin(lat0)*np.sin(c)))*180/np.pi 
	    lat=np.arcsin(np.cos(c)*np.sin(lat0)+y*np.sin(c)*np.cos(lat0)/r)*180/np.pi
	    lat=90+lat
	    longi=-longi+180
	    if longi>90:
	     longi=longi-180
	    c=str(np.around(longi,decimals=1))+str(',')+str(np.around(lat,decimals=1))
	    ui.coord_label.setText(str(c))

########################################################
#
# Calculate interplanar distance 
#
#######################################################

def dhkl():
    i=np.float(ui.pole1_entry.text())
    j=np.float(ui.pole2_entry.text())
    k=np.float(ui.pole3_entry.text())
    a=np.float(ui.a_entry.text())
    b=np.float(ui.b_entry.text())
    c=np.float(ui.c_entry.text())
    alp=np.float(ui.alpha_entry.text())
    bet=np.float(ui.beta_entry.text())
    gam=np.float(ui.gamma_entry.text())
    alp=alp*np.pi/180;
    bet=bet*np.pi/180;
    gam=gam*np.pi/180;
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
    d=np.around(1/(np.sqrt(np.dot(np.array([i,j,k]),np.dot(np.linalg.inv(G),np.array([i,j,k]))))), decimals=3)
    ui.dhkl_label.setText(str(d))
    return        

####################################################################
#
# Reset view after zoom
#
#################################################################### 

def reset_view():
	global a
	
	a.axis([minx,maxx,miny,maxy])
	trace()
	
####################################################################
#
# Enable or disable Wulff net
#
#################################################################### 

def wulff():
	global a
	if ui.wulff_button.isChecked():
		fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
		img=Image.open(fn)
		img= np.array(img)
	else:
		img = 255*np.ones([600,600,3],dtype=np.uint8)
		circle = plt.Circle((300, 300), 300, color='black',fill=False)
		a.add_artist(circle)
		a.plot(300,300,'+',markersize=10,mew=3,color='black')
	
	a.imshow(img,interpolation="bicubic")
	a.axis('off')
    	a.figure.canvas.draw()  
	


 #######################################################################               
#######################################################################
#
# Main
#
#####################################################################    


####################################################################
#
# Refresh action on stereo
#
####################################################################

def trace():
    global T,x,y,z,axes,axesh,M,trP,a
    minx,maxx=a.get_xlim()
    miny,maxy=a.get_ylim()
    a = figure.add_subplot(111) 
    a.figure.clear()
    a = figure.add_subplot(111)
    P=np.zeros((axes.shape[0],2))
    T=np.zeros((axes.shape))
    C=[]
    
    trace_plan2(trP)			
    schmid_trace2(tr_schmid)
    
    for i in range(0,axes.shape[0]):
        axeshr=np.array([axesh[i,0],axesh[i,1],axesh[i,2]])

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
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=np.float(ui.size_var.text()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=np.float(ui.size_var.text()), facecolors='none', linewidths=1.5)       
    
    a.axis([minx,maxx,miny,maxy])
    wulff()
    
 ####################################
 #
 # Initial plot from a given diffraction
 #
 ####################################
 
def princ():
    global T,angle_alpha, angle_beta, angle_z,M,Dstar,D,g,M0,trP,axeshr,nn,a,minx,maxx,miny,maxy
    trP=np.zeros((1,5))
    crist() 
    a = figure.add_subplot(111)
    a.figure.clear()
    a = figure.add_subplot(111)
    
    
    diff1=np.float(ui.diff1_entry.text())
    diff2=np.float(ui.diff2_entry.text())
    diff3=np.float(ui.diff3_entry.text())
    tilt=np.float(ui.tilt_entry.text())
    inclinaison=np.float(ui.inclinaison_entry.text())    
     
    d0=np.array([diff1,diff2,diff3])
    if var_uvw()==0: 
       d=np.dot(Dstar,d0)
              
    else:
       d=np.dot(Dstar,d0)
    if diff2==0 and diff1==0:
        normal=np.array([1,0,0])
        ang=np.pi/2
    else:
        normal=np.array([-d[2],0,d[0]])
        ang=np.arccos(np.dot(d,np.array([0,1,0]))/np.linalg.norm(d))
    
    R=np.dot(Rot(-tilt,0,1,0),np.dot(Rot(-inclinaison,0,0,1),Rot(ang*180/np.pi, normal[0],normal[1],normal[2])))    
       
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
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=np.float(ui.size_var.text()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=np.float(ui.size_var.text()), facecolors='none', linewidths=1.5)            
    
    minx,maxx=-2,602
    miny,maxy=-2,602
    a.axis([minx,maxx,miny,maxy])
    wulff()
    
    angle_alpha=0             
    angle_beta=0
    angle_z=0
    g=0
    ui.angle_alpha_label_2.setText('0.0')
    ui.angle_beta_label_2.setText('0.0')
    ui.angle_z_label_2.setText('0.0')
    ui.angle_beta_label_2.setText('0.0')
    ui.angle_z_label_2.setText('0.0')    
    ui.rg_label.setText('0.0')
    M=R
    M0=R
    phir=np.arccos(M[2,2])*180/np.pi
    phi2r=np.arctan2(M[2,0],M[2,1])*180/np.pi
    phi1r=np.arctan2(M[0,2],-M[1,2])*180/np.pi
    t=str(np.around(phi1r,decimals=1))+str(',')+str(np.around(phir,decimals=1))+str(',')+str(np.around(phi2r,decimals=1))
    
    ui.angle_euler_label.setText(t)
    return T,angle_alpha,angle_beta,angle_z,g,M,M0
    
##############################################"
#
# Plot from Euler angle
#
##################################################"

def princ2():
    global T,angle_alpha,angle_beta,angle_z,M,Dstar,D,g,M0,trP,a,axeshr,nn,minx,maxx,miny,maxy
    
    trP=np.zeros((1,5))
    a = figure.add_subplot(111)
    a.figure.clear()
    a = figure.add_subplot(111)
    phi1=np.float(ui.phi1_entry.text())
    phi=np.float(ui.phi_entry.text())
    phi2=np.float(ui.phi2_entry.text())
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=np.array(Image.open(fn))
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
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,c=C,s=np.float(ui.size_var.text()))
    else:
        a.scatter(P[:,0]+600/2,P[:,1]+600/2,edgecolor=C, s=np.float(ui.size_var.text()), facecolors='none', linewidths=1.5)       
    minx,maxx=-2,602
    miny,maxy=-2,602
    a.axis([minx,maxx,miny,maxy])
    wulff()
    
    angle_alpha=0             
    angle_beta=0
    angle_z=0
    g=0
    ui.angle_alpha_label_2.setText('0.0')
    ui.angle_beta_label_2.setText('0.0')
    ui.angle_z_label_2.setText('0.0')
    ui.angle_beta_label_2.setText('0.0')
    ui.angle_z_label_2.setText('0.0')    
    ui.rg_label.setText('0.0')
    M=rotation(phi1,phi,phi2)
    t=str(np.around(phi1,decimals=1))+str(',')+str(np.around(phi,decimals=1))+str(',')+str(np.around(phi2,decimals=1))
    ui.angle_euler_label.setText(t)
    
    return T,angle_alpha,angle_beta,angle_z,g,M

  
  
#######################################################################
#######################################################################
#
# GUI
#
#######################################################################




######################################################
#
# Menu
#
##########################################################

###########################################################
#
# Structure
#
##############################################################

def structure(item):
    global x0, var_hexa, d_label_var, e_entry
    
    ui.a_entry.setText(str(item[1]))
    ui.b_entry.setText(str(item[2]))
    ui.c_entry.setText(str(item[3]))
    ui.alpha_entry.setText(str(item[4]))
    ui.beta_entry.setText(str(item[5]))
    ui.gamma_entry.setText(str(item[6]))
    if eval(item[4])==90 and eval(item[5])==90 and eval(item[6])==120 :
        ui.hexa_button.setChecked(True)
        ui.e_entry.setText('2')
        ui.d_label_var.setText('3')
    else:
        ui.d_entry.setText('1')
        ui.e_entry.setText('1')
    
 

####################################################################
#
# Class for dialog box for measuring angle between two poles
#
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
       if ui.uvw_button.isChecked==True: 
           c1c=np.dot(Dstar,c1)
           c2c=np.dot(Dstar,c2)
       else:
           c1c=np.dot(Dstar,c1)
           c2c=np.dot(Dstar,c2)
       the=np.arccos(np.dot(c1c,c2c)/(np.linalg.norm(c1c)*np.linalg.norm(c2c)))                   
       thes=str(np.around(the*180/np.pi,decimals=2))        
       self.label.setText(thes)
                        
        

 
##################################################
#
# Class for dialog box for Schmid factor calculation
#
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

                        
 #######################################
 #
 # Save stereo as png
 #
 ############################################"

def image_save():
    filename=QtGui.QFileDialog.getSaveFileName( Index,"Save file", "", ".png")
    pixmap = QtGui.QPixmap.grabWidget(canvas,55,49,710,710)
    pixmap.save(str(filename)+".png")
     
     
##################################################
#
# Class for dialog box for calculating x,y,z directions
#
################################################### 

class xyzCalc(QtGui.QDialog):

       
    def __init__(self, parent=None):
        
        super(xyzCalc, self).__init__(parent)
        
	self.gridLayout = QtGui.QGridLayout(self)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.X_label = QtGui.QLabel(self)
        self.X_label.setObjectName(_fromUtf8("X_label"))
        self.X_label.setText(_fromUtf8("X"))
        self.gridLayout.addWidget(self.X_label, 0, 0, 1, 1)
        self.Z_label = QtGui.QLabel(self)
        self.Z_label.setObjectName(_fromUtf8("Z_label"))
        self.Z_label.setText(_fromUtf8("Z"))
        self.gridLayout.addWidget(self.Z_label, 4, 0, 1, 1)
        self.X_label2 = QtGui.QLabel(self)
        self.X_label2.setText(_fromUtf8(""))
        self.X_label2.setObjectName(_fromUtf8("X_label2"))
        self.gridLayout.addWidget(self.X_label2, 1, 0, 1, 1)
        self.Y_label = QtGui.QLabel(self)
        self.Y_label.setObjectName(_fromUtf8("Y_label"))
        self.Y_label.setText(_fromUtf8("Y"))
        self.gridLayout.addWidget(self.Y_label, 2, 0, 1, 1)
        self.Y_label2 = QtGui.QLabel(self)
        self.Y_label2.setText(_fromUtf8(""))
        self.Y_label2.setObjectName(_fromUtf8("Y_label2"))
        self.gridLayout.addWidget(self.Y_label2, 3, 0, 1, 1)
        self.Z_label2 = QtGui.QLabel(self)
        self.Z_label2.setText(_fromUtf8(""))
        self.Z_label2.setObjectName(_fromUtf8("Z_label2"))
        self.gridLayout.addWidget(self.Z_label2, 5, 0, 1, 1)
        self.buttonBox = QtGui.QPushButton(self)
        self.buttonBox.setText("Update")
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout.addWidget(self.buttonBox, 6, 0, 1, 1)
        self.buttonBox.clicked.connect(self.center)
   
    def center(self):
	global D, Dstar,M
	A=np.dot(np.linalg.inv(M),np.array([0,0,1]))
	A2=np.dot(np.linalg.inv(M),np.array([1,0,0]))
	A3=np.dot(np.linalg.inv(M),np.array([0,1,0]))
	C=np.dot(np.linalg.inv(Dstar),A)
	Zp=C/np.linalg.norm(C)
	C2=np.dot(np.linalg.inv(Dstar),A2)
	Xp=C2/np.linalg.norm(C2)
	C3=np.dot(np.linalg.inv(Dstar),A3)
	Yp=C3/np.linalg.norm(C3)    
	
	self.X_label2.setText(str(Xp[0])+', '+str(Xp[1])+', '+str(Xp[2]))
	self.Y_label2.setText(str(Yp[0])+', '+str(Yp[1])+', '+str(Yp[2]))
	self.Z_label2.setText(str(Zp[0])+', '+str(Zp[1])+', '+str(Zp[2]))
            

##################################################
#
# Add matplotlib toolbar to zoom and pan
#
################################################### 


class NavigationToolbar(NavigationToolbar):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Pan', 'Zoom')]
    def set_message(self, msg):
        pass


#############################################################
#
# Launch
#
#############################################################"
    
try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)
    
if __name__ == "__main__":
    
	app = QtGui.QApplication(sys.argv)
	Index = QtGui.QMainWindow()	
	ui = stereoprojUI.Ui_StereoProj()
	ui.setupUi(Index)
	figure=plt.figure()
	canvas=FigureCanvas(figure)
	ui.mplvl.addWidget(canvas)
	toolbar = NavigationToolbar(canvas, canvas)
	toolbar.setMinimumWidth(601)

# Read structure file
	
	file_struct=open(os.path.join(os.path.dirname(__file__), 'structure.txt') ,"r")
	x0=[]
	for line in file_struct:
	    x0.append(map(str, line.split()))
	i=0
	file_struct.close()            
	for item in x0:
	    entry = ui.menuStructure.addAction(item[0])
	    Index.connect(entry,QtCore.SIGNAL('triggered()'), lambda item=item: structure(item))
	    i=i+1

# Connect dialog boxes and buttons

	dialogSchmidCalc = SchmidCalc()
	dialogSchmidCalc.setWindowTitle("Calculate Schmid Factor")
	Index.connect(ui.actionCalculate_Schmid_factor, QtCore.SIGNAL('triggered()'), dialogSchmidCalc.show) 

	Index.connect(ui.actionSave_figure, QtCore.SIGNAL('triggered()'), image_save) 

	dialogAngleCalc = AngleCalc()
	dialogAngleCalc.setWindowTitle("Calculate Angle")
	Index.connect(ui.actionCalculate_angle, QtCore.SIGNAL('triggered()'), dialogAngleCalc.show) 	

	dialog_largeur_plan=LargeurPlanCalc()
	dialog_largeur_plan.setWindowTitle("Apparent Width")
	Index.connect(ui.actionCalculate_apparent_width, QtCore.SIGNAL('triggered()'), dialog_largeur_plan.show) 	
	dialog_xyz=xyzCalc()
	dialog_xyz.setWindowTitle("X,Y,Z directions")
	Index.connect(ui.actionCalculate_xyz, QtCore.SIGNAL('triggered()'), dialog_xyz.show) 
	
	ui.button_trace2.clicked.connect(princ2)
	ui.button_trace.clicked.connect(princ)
	ui.angle_alpha_buttonp.clicked.connect(rot_alpha_p)
	ui.angle_alpha_buttonm.clicked.connect(rot_alpha_m)
	ui.angle_beta_buttonp.clicked.connect(rot_beta_p)
	ui.angle_beta_buttonm.clicked.connect(rot_beta_m)
	ui.angle_z_buttonp.clicked.connect(rot_z_p)
	ui.angle_z_buttonm.clicked.connect(rot_z_m)
	ui.rot_gm_button.clicked.connect(rotgm)
	ui.rot_gp_button.clicked.connect(rotgp)
	ui.lock_checkButton.stateChanged.connect(lock)
	ui.addpole_button.clicked.connect(addpole)
	ui.undo_addpole_button.clicked.connect(undo_addpole)
	ui.sym_button.clicked.connect(addpole_sym)
	ui.undo_sym_button.clicked.connect(undo_sym)
	ui.trace_plan_button.clicked.connect(trace_addplan)
	ui.undo_trace_plan_button.clicked.connect(undo_trace_addplan)
	ui.trace_plan_sym_button.clicked.connect(trace_plan_sym)
	ui.undo_trace_plan_sym_button.clicked.connect(undo_trace_plan_sym)
	ui.trace_schmid_button.clicked.connect(schmid_trace)
	ui.undo_trace_schmid.clicked.connect(undo_schmid_trace)
	ui.norm_button.clicked.connect(dhkl)
	ui.dm_button.clicked.connect(dm)
	ui.dp_button.clicked.connect(dp)
	ui.reset_view_button.clicked.connect(reset_view)
	figure.canvas.mpl_connect('motion_notify_event', coordinates)
	figure.canvas.mpl_connect('button_press_event', click_a_pole)

# Initialize variables
	
	dmip=0
	tr_schmid=np.zeros((1,3))
	var_lock=0
	ui.lock_checkButton.setChecked(False)
	ui.color_trace_bleu.setChecked(True)
	ui.wulff_button.setChecked(True)
	ui.wulff_button.setChecked(True)
	ui.d_label_var.setText('0')
	ui.alpha_entry.setText('90')
	ui.beta_entry.setText('90')
	ui.gamma_entry.setText('90')
	ui.a_entry.setText('1')
	ui.b_entry.setText('1')
	ui.c_entry.setText('1')
	ui.phi1_entry.setText('0')
	ui.phi_entry.setText('0')
	ui.phi2_entry.setText('0')
	ui.e_entry.setText('1')
	ui.rg_label.setText('0.0')
	ui.angle_euler_label.setText(' ')
	ui.size_var.setText('40')
	ui.e_entry.setText('1')
	ui.angle_alpha_entry.setText('5')
	ui.angle_beta_entry.setText('5')
	ui.angle_z_entry.setText('5')
	ui.angle_beta_entry.setText('5')
	ui.angle_z_entry.setText('5')	
	ui.d_entry.setText('1')
	ui.rot_g_entry.setText('5')
	
	a = figure.add_subplot(111)
	wulff()	
	Index.show()
	sys.exit(app.exec_())




