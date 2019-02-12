# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'hkl.ui'
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

class Ui_hkl_uvw(object):
    def setupUi(self, hkl_uvw):
        hkl_uvw.setObjectName(_fromUtf8("hkl_uvw"))
        hkl_uvw.resize(290, 483)
        self.layoutWidget = QtGui.QWidget(hkl_uvw)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 10, 266, 462))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.v = QtGui.QLineEdit(self.layoutWidget)
        self.v.setObjectName(_fromUtf8("v"))
        self.gridLayout.addWidget(self.v, 9, 0, 1, 1)
        self.hkl_entry = QtGui.QLineEdit(self.layoutWidget)
        self.hkl_entry.setObjectName(_fromUtf8("hkl_entry"))
        self.gridLayout.addWidget(self.hkl_entry, 13, 0, 1, 1)
        self.uvw_entry = QtGui.QLineEdit(self.layoutWidget)
        self.uvw_entry.setObjectName(_fromUtf8("uvw_entry"))
        self.gridLayout.addWidget(self.uvw_entry, 5, 0, 1, 1)
        self.pushButton_to_hkl = QtGui.QPushButton(self.layoutWidget)
        self.pushButton_to_hkl.setObjectName(_fromUtf8("pushButton_to_hkl"))
        self.gridLayout.addWidget(self.pushButton_to_hkl, 11, 0, 1, 1)
        self.pushButton_copy1 = QtGui.QPushButton(self.layoutWidget)
        self.pushButton_copy1.setObjectName(_fromUtf8("pushButton_copy1"))
        self.gridLayout.addWidget(self.pushButton_copy1, 6, 0, 1, 1)
        self.l = QtGui.QLineEdit(self.layoutWidget)
        self.l.setObjectName(_fromUtf8("l"))
        self.gridLayout.addWidget(self.l, 3, 0, 1, 1)
        self.w = QtGui.QLineEdit(self.layoutWidget)
        self.w.setObjectName(_fromUtf8("w"))
        self.gridLayout.addWidget(self.w, 10, 0, 1, 1)
        self.h = QtGui.QLineEdit(self.layoutWidget)
        self.h.setObjectName(_fromUtf8("h"))
        self.gridLayout.addWidget(self.h, 1, 0, 1, 1)
        self.uvw_label = QtGui.QLabel(self.layoutWidget)
        self.uvw_label.setObjectName(_fromUtf8("uvw_label"))
        self.gridLayout.addWidget(self.uvw_label, 7, 0, 1, 1)
        self.u = QtGui.QLineEdit(self.layoutWidget)
        self.u.setObjectName(_fromUtf8("u"))
        self.gridLayout.addWidget(self.u, 8, 0, 1, 1)
        self.k = QtGui.QLineEdit(self.layoutWidget)
        self.k.setObjectName(_fromUtf8("k"))
        self.gridLayout.addWidget(self.k, 2, 0, 1, 1)
        self.hkl_label = QtGui.QLabel(self.layoutWidget)
        self.hkl_label.setObjectName(_fromUtf8("hkl_label"))
        self.gridLayout.addWidget(self.hkl_label, 0, 0, 1, 1)
        self.pushButton_to_uvw = QtGui.QPushButton(self.layoutWidget)
        self.pushButton_to_uvw.setObjectName(_fromUtf8("pushButton_to_uvw"))
        self.gridLayout.addWidget(self.pushButton_to_uvw, 4, 0, 1, 1)
        self.pushButton_copy2 = QtGui.QPushButton(self.layoutWidget)
        self.pushButton_copy2.setObjectName(_fromUtf8("pushButton_copy2"))
        self.gridLayout.addWidget(self.pushButton_copy2, 14, 0, 1, 1)

        self.retranslateUi(hkl_uvw)
        QtCore.QMetaObject.connectSlotsByName(hkl_uvw)
        hkl_uvw.setTabOrder(self.h, self.k)
        hkl_uvw.setTabOrder(self.k, self.l)
        hkl_uvw.setTabOrder(self.l, self.pushButton_to_uvw)
        hkl_uvw.setTabOrder(self.pushButton_to_uvw, self.uvw_entry)
        hkl_uvw.setTabOrder(self.uvw_entry, self.pushButton_copy1)
        hkl_uvw.setTabOrder(self.pushButton_copy1, self.u)
        hkl_uvw.setTabOrder(self.u, self.v)
        hkl_uvw.setTabOrder(self.v, self.w)
        hkl_uvw.setTabOrder(self.w, self.pushButton_to_hkl)
        hkl_uvw.setTabOrder(self.pushButton_to_hkl, self.hkl_entry)
        hkl_uvw.setTabOrder(self.hkl_entry, self.pushButton_copy2)

    def retranslateUi(self, hkl_uvw):
        hkl_uvw.setWindowTitle(_translate("hkl_uvw", "hkl-uvw", None))
        self.pushButton_to_hkl.setText(_translate("hkl_uvw", "To hkl", None))
        self.pushButton_copy1.setText(_translate("hkl_uvw", "Copy to pole", None))
        self.uvw_label.setText(_translate("hkl_uvw", "uvw", None))
        self.hkl_label.setText(_translate("hkl_uvw", "hkl", None))
        self.pushButton_to_uvw.setText(_translate("hkl_uvw", "To uvw", None))
        self.pushButton_copy2.setText(_translate("hkl_uvw", "Copy to pole", None))

