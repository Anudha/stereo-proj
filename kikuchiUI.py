# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'kikuchi.ui'
#
# Created by: PyQt4 UI code generator 4.12.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

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

class Ui_Kikuchi(object):
    def setupUi(self, Kikuchi):
        Kikuchi.setObjectName(_fromUtf8("Kikuchi"))
        Kikuchi.setWindowModality(QtCore.Qt.NonModal)
        Kikuchi.resize(809, 793)
        self.layoutWidget = QtGui.QWidget(Kikuchi)
        self.layoutWidget.setGeometry(QtCore.QRect(9, 9, 791, 771))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.angle_label = QtGui.QLabel(self.layoutWidget)
        self.angle_label.setObjectName(_fromUtf8("angle_label"))
        self.gridLayout.addWidget(self.angle_label, 1, 2, 1, 1)
        self.E_label = QtGui.QLabel(self.layoutWidget)
        self.E_label.setObjectName(_fromUtf8("E_label"))
        self.gridLayout.addWidget(self.E_label, 1, 0, 1, 1)
        self.buttonBox = QtGui.QDialogButtonBox(self.layoutWidget)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout.addWidget(self.buttonBox, 2, 2, 1, 2)
        self.mplwindow = QtGui.QWidget(self.layoutWidget)
        self.mplwindow.setObjectName(_fromUtf8("mplwindow"))
        self.gridLayoutWidget = QtGui.QWidget(self.mplwindow)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(0, 0, 781, 691))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.mplvl = QtGui.QGridLayout(self.gridLayoutWidget)
        self.mplvl.setMargin(0)
        self.mplvl.setObjectName(_fromUtf8("mplvl"))
        self.gridLayout.addWidget(self.mplwindow, 0, 0, 1, 4)
        self.E_entry = QtGui.QLineEdit(self.layoutWidget)
        self.E_entry.setObjectName(_fromUtf8("E_entry"))
        self.gridLayout.addWidget(self.E_entry, 1, 1, 1, 1)
        self.angle_entry = QtGui.QLineEdit(self.layoutWidget)
        self.angle_entry.setObjectName(_fromUtf8("angle_entry"))
        self.gridLayout.addWidget(self.angle_entry, 1, 3, 1, 1)
        self.label_checkBox = QtGui.QCheckBox(self.layoutWidget)
        self.label_checkBox.setObjectName(_fromUtf8("label_checkBox"))
        self.gridLayout.addWidget(self.label_checkBox, 2, 0, 1, 2)

        self.retranslateUi(Kikuchi)
        QtCore.QMetaObject.connectSlotsByName(Kikuchi)
        Kikuchi.setTabOrder(self.E_entry, self.buttonBox)

    def retranslateUi(self, Kikuchi):
        Kikuchi.setWindowTitle(_translate("Kikuchi", "Kikuchi lines", None))
        self.angle_label.setText(_translate("Kikuchi", "Aperture angle", None))
        self.E_label.setText(_translate("Kikuchi", "E (kV)", None))
        self.label_checkBox.setText(_translate("Kikuchi", "labels", None))

