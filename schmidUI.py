# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'schmid.ui'
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

class Ui_Schmid(object):
    def setupUi(self, Schmid):
        Schmid.setObjectName(_fromUtf8("Schmid"))
        Schmid.resize(352, 326)
        self.layoutWidget = QtGui.QWidget(Schmid)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 10, 318, 298))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.b1 = QtGui.QLineEdit(self.layoutWidget)
        self.b1.setObjectName(_fromUtf8("b1"))
        self.gridLayout.addWidget(self.b1, 2, 0, 1, 1)
        self.schmid_factor_label = QtGui.QLabel(self.layoutWidget)
        self.schmid_factor_label.setText(_fromUtf8(""))
        self.schmid_factor_label.setObjectName(_fromUtf8("schmid_factor_label"))
        self.gridLayout.addWidget(self.schmid_factor_label, 4, 0, 1, 3)
        self.buttonBox = QtGui.QDialogButtonBox(self.layoutWidget)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout.addWidget(self.buttonBox, 6, 0, 1, 2)
        self.n_label = QtGui.QLabel(self.layoutWidget)
        self.n_label.setObjectName(_fromUtf8("n_label"))
        self.gridLayout.addWidget(self.n_label, 0, 1, 1, 1)
        self.n2 = QtGui.QLineEdit(self.layoutWidget)
        self.n2.setObjectName(_fromUtf8("n2"))
        self.gridLayout.addWidget(self.n2, 3, 1, 1, 1)
        self.T_label = QtGui.QLabel(self.layoutWidget)
        self.T_label.setObjectName(_fromUtf8("T_label"))
        self.gridLayout.addWidget(self.T_label, 0, 2, 1, 1)
        self.n1 = QtGui.QLineEdit(self.layoutWidget)
        self.n1.setObjectName(_fromUtf8("n1"))
        self.gridLayout.addWidget(self.n1, 2, 1, 1, 1)
        self.b2 = QtGui.QLineEdit(self.layoutWidget)
        self.b2.setObjectName(_fromUtf8("b2"))
        self.gridLayout.addWidget(self.b2, 3, 0, 1, 1)
        self.n0 = QtGui.QLineEdit(self.layoutWidget)
        self.n0.setObjectName(_fromUtf8("n0"))
        self.gridLayout.addWidget(self.n0, 1, 1, 1, 1)
        self.T0 = QtGui.QLineEdit(self.layoutWidget)
        self.T0.setObjectName(_fromUtf8("T0"))
        self.gridLayout.addWidget(self.T0, 1, 2, 1, 1)
        self.T1 = QtGui.QLineEdit(self.layoutWidget)
        self.T1.setObjectName(_fromUtf8("T1"))
        self.gridLayout.addWidget(self.T1, 2, 2, 1, 1)
        self.T2 = QtGui.QLineEdit(self.layoutWidget)
        self.T2.setObjectName(_fromUtf8("T2"))
        self.gridLayout.addWidget(self.T2, 3, 2, 1, 1)
        self.schmid_text = QtGui.QTextEdit(self.layoutWidget)
        self.schmid_text.setObjectName(_fromUtf8("schmid_text"))
        self.gridLayout.addWidget(self.schmid_text, 5, 0, 1, 3)
        self.b_label = QtGui.QLabel(self.layoutWidget)
        self.b_label.setObjectName(_fromUtf8("b_label"))
        self.gridLayout.addWidget(self.b_label, 0, 0, 1, 1)
        self.b0 = QtGui.QLineEdit(self.layoutWidget)
        self.b0.setObjectName(_fromUtf8("b0"))
        self.gridLayout.addWidget(self.b0, 1, 0, 1, 1)

        self.retranslateUi(Schmid)
        QtCore.QMetaObject.connectSlotsByName(Schmid)
        Schmid.setTabOrder(self.b0, self.b1)
        Schmid.setTabOrder(self.b1, self.b2)
        Schmid.setTabOrder(self.b2, self.n0)
        Schmid.setTabOrder(self.n0, self.n1)
        Schmid.setTabOrder(self.n1, self.n2)
        Schmid.setTabOrder(self.n2, self.buttonBox)

    def retranslateUi(self, Schmid):
        Schmid.setWindowTitle(_translate("Schmid", "Schmid Factor", None))
        self.n_label.setText(_translate("Schmid", "n", None))
        self.T_label.setText(_translate("Schmid", "T", None))
        self.b_label.setText(_translate("Schmid", "b", None))

