# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'angle.ui'
#
# Created by: PyQt4 UI code generator 4.12.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
import numpy as np

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

class Ui_Angle(object):
    def setupUi(self, Angle):
        Angle.setObjectName(_fromUtf8("Angle"))
        Angle.resize(288, 220)
        self.layoutWidget = QtGui.QWidget(Angle)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 10, 266, 205))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.angle_label = QtGui.QLabel(self.layoutWidget)
        self.angle_label.setText(_fromUtf8(""))
        self.angle_label.setObjectName(_fromUtf8("angle_label"))
        self.gridLayout.addWidget(self.angle_label, 4, 0, 1, 2)
        self.n22 = QtGui.QLineEdit(self.layoutWidget)
        self.n22.setObjectName(_fromUtf8("n22"))
        self.gridLayout.addWidget(self.n22, 3, 1, 1, 1)
        self.n12 = QtGui.QLineEdit(self.layoutWidget)
        self.n12.setObjectName(_fromUtf8("n12"))
        self.gridLayout.addWidget(self.n12, 3, 0, 1, 1)
        self.n10 = QtGui.QLineEdit(self.layoutWidget)
        self.n10.setObjectName(_fromUtf8("n10"))
        self.gridLayout.addWidget(self.n10, 1, 0, 1, 1)
        self.n2_label = QtGui.QLabel(self.layoutWidget)
        self.n2_label.setObjectName(_fromUtf8("n2_label"))
        self.gridLayout.addWidget(self.n2_label, 0, 1, 1, 1)
        self.buttonBox = QtGui.QDialogButtonBox(self.layoutWidget)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout.addWidget(self.buttonBox, 5, 0, 1, 2)
        self.n21 = QtGui.QLineEdit(self.layoutWidget)
        self.n21.setObjectName(_fromUtf8("n21"))
        self.gridLayout.addWidget(self.n21, 2, 1, 1, 1)
        self.n11 = QtGui.QLineEdit(self.layoutWidget)
        self.n11.setObjectName(_fromUtf8("n11"))
        self.gridLayout.addWidget(self.n11, 2, 0, 1, 1)
        self.n20 = QtGui.QLineEdit(self.layoutWidget)
        self.n20.setObjectName(_fromUtf8("n20"))
        self.gridLayout.addWidget(self.n20, 1, 1, 1, 1)
        self.n1_label = QtGui.QLabel(self.layoutWidget)
        self.n1_label.setObjectName(_fromUtf8("n1_label"))
        self.gridLayout.addWidget(self.n1_label, 0, 0, 1, 1)

	
        self.retranslateUi(Angle)
        QtCore.QMetaObject.connectSlotsByName(Angle)
        Angle.setTabOrder(self.n10, self.n11)
        Angle.setTabOrder(self.n11, self.n12)
        Angle.setTabOrder(self.n12, self.n20)
        Angle.setTabOrder(self.n20, self.n21)
        Angle.setTabOrder(self.n21, self.n22)
        Angle.setTabOrder(self.n22, self.buttonBox)

    def retranslateUi(self, Angle):
        Angle.setWindowTitle(_translate("Angle", "Angle", None))
        self.n2_label.setText(_translate("Angle", "n2", None))
        self.n1_label.setText(_translate("Angle", "n1", None))
        
       
       


