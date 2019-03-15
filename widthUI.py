# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'width.ui'
#
# Created by: PyQt4 UI code generator 4.12.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import pyplot as plt

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

class NavigationToolbar(NavigationToolbar):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Pan', 'Zoom')]
    def set_message(self, msg):
        pass
        
class Ui_Width(object):
    def setupUi(self, Width):
        Width.setObjectName(_fromUtf8("Width"))
        Width.resize(682, 615)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, None)
        
        self.layoutWidget = QtGui.QWidget(Width)
        self.layoutWidget.setGeometry(QtCore.QRect(9, 9, 661, 601))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.n0 = QtGui.QLineEdit(self.layoutWidget)
        self.n0.setObjectName(_fromUtf8("n0"))
        self.gridLayout.addWidget(self.n0, 1, 0, 1, 1)
        self.widget_1 = QtGui.QWidget(self.layoutWidget)
        self.widget_1.setObjectName(_fromUtf8("widget_1"))
        self.gridLayout.addWidget(self.canvas, 0, 0, 1, 6)
        self.n1 = QtGui.QLineEdit(self.layoutWidget)
        self.n1.setObjectName(_fromUtf8("n1"))
        self.gridLayout.addWidget(self.n1, 1, 1, 1, 1)
        self.n2 = QtGui.QLineEdit(self.layoutWidget)
        self.n2.setObjectName(_fromUtf8("n2"))
        self.gridLayout.addWidget(self.n2, 1, 2, 1, 1)
        self.widget_2 = QtGui.QWidget(self.layoutWidget)
        self.widget_2.setObjectName(_fromUtf8("widget_2"))
        self.gridLayout.addWidget(self.toolbar, 3, 0, 1, 2)
        self.buttonBox = QtGui.QDialogButtonBox(self.layoutWidget)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.gridLayout.addWidget(self.buttonBox, 3, 2, 1, 4)
        self.thickness = QtGui.QLineEdit(self.layoutWidget)
        self.thickness.setObjectName(_fromUtf8("thickness"))
        self.gridLayout.addWidget(self.thickness, 2, 0, 1, 1)
        self.thickness_label = QtGui.QLabel(self.layoutWidget)
        self.thickness_label.setObjectName(_fromUtf8("thickness_label"))
        self.gridLayout.addWidget(self.thickness_label, 2, 1, 1, 1)
        self.trace_radio_button = QtGui.QRadioButton(self.layoutWidget)
        self.trace_radio_button.setObjectName(_fromUtf8("trace_radio_button"))
        self.gridLayout.addWidget(self.trace_radio_button, 1, 3, 1, 1)
        self.thickness_checkBox = QtGui.QCheckBox(self.layoutWidget)
        self.thickness_checkBox.setObjectName(_fromUtf8("thickness_checkBox"))
        self.gridLayout.addWidget(self.thickness_checkBox, 2, 2, 1, 1)

        self.retranslateUi(Width)
        QtCore.QMetaObject.connectSlotsByName(Width)
        Width.setTabOrder(self.n0, self.n1)
        Width.setTabOrder(self.n1, self.n2)
        Width.setTabOrder(self.n2, self.trace_radio_button)
        Width.setTabOrder(self.trace_radio_button, self.thickness)
        Width.setTabOrder(self.thickness, self.thickness_checkBox)
        Width.setTabOrder(self.thickness_checkBox, self.buttonBox)

    def retranslateUi(self, Width):
        Width.setWindowTitle(_translate("Width", "Apparent Width", None))
        self.thickness_label.setText(_translate("Width", "thickness (nm)", None))
        self.trace_radio_button.setText(_translate("Width", "trace dir.", None))
        self.thickness_checkBox.setText(_translate("Width", "w (nm)", None))

