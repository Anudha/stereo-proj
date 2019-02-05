# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'xyz.ui'
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

class Ui_xyz_dialog(object):
    def setupUi(self, xyz_dialog):
        xyz_dialog.setObjectName(_fromUtf8("xyz_dialog"))
        xyz_dialog.resize(399, 379)
        self.gridLayoutWidget = QtGui.QWidget(xyz_dialog)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 10, 381, 366))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.Y_label = QtGui.QLabel(self.gridLayoutWidget)
        self.Y_label.setObjectName(_fromUtf8("Y_label"))
        self.gridLayout.addWidget(self.Y_label, 2, 0, 1, 1)
        self.X_label = QtGui.QLabel(self.gridLayoutWidget)
        self.X_label.setObjectName(_fromUtf8("X_label"))
        self.gridLayout.addWidget(self.X_label, 0, 0, 1, 1)
        self.Z_label = QtGui.QLabel(self.gridLayoutWidget)
        self.Z_label.setObjectName(_fromUtf8("Z_label"))
        self.gridLayout.addWidget(self.Z_label, 4, 0, 1, 1)
        self.xyz_button = QtGui.QPushButton(self.gridLayoutWidget)
        self.xyz_button.setObjectName(_fromUtf8("xyz_button"))
        self.gridLayout.addWidget(self.xyz_button, 6, 0, 1, 1)
        self.X_text = QtGui.QTextEdit(self.gridLayoutWidget)
        self.X_text.setObjectName(_fromUtf8("X_text"))
        self.gridLayout.addWidget(self.X_text, 1, 0, 1, 1)
        self.Y_text = QtGui.QTextEdit(self.gridLayoutWidget)
        self.Y_text.setObjectName(_fromUtf8("Y_text"))
        self.gridLayout.addWidget(self.Y_text, 3, 0, 1, 1)
        self.Z_text = QtGui.QTextEdit(self.gridLayoutWidget)
        self.Z_text.setObjectName(_fromUtf8("Z_text"))
        self.gridLayout.addWidget(self.Z_text, 5, 0, 1, 1)

        self.retranslateUi(xyz_dialog)
        QtCore.QMetaObject.connectSlotsByName(xyz_dialog)

    def retranslateUi(self, xyz_dialog):
        xyz_dialog.setWindowTitle(_translate("xyz_dialog", "xyz directions", None))
        self.Y_label.setText(_translate("xyz_dialog", "Y", None))
        self.X_label.setText(_translate("xyz_dialog", "X", None))
        self.Z_label.setText(_translate("xyz_dialog", "Z", None))
        self.xyz_button.setText(_translate("xyz_dialog", "Update", None))

