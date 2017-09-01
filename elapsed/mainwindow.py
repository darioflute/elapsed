#!/usr/bin/env python

# Essential libraries
import sys,os
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Matplotlib
import matplotlib
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1


# Connection with PyQt5
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QWidget, QMainWindow, QMessageBox,QToolBar,QAction,QStatusBar,QSizePolicy,
                             QHBoxLayout, QVBoxLayout, QApplication, QListWidget,QSplitter,QMenu,QTabWidget)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt, QTimer, QThread, pyqtSignal, QObject


class MplCanvas(FigureCanvas):
    """ Basic matplotlib canvas class """

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.compute_initial_figure()

    def compute_initial_figure(self):
        pass

class ImageCanvas(MplCanvas):
    """ Canvas to plot an image """
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)

    def compute_initial_figure(self):
        self.axes1 = self.fig.add_subplot(111)
        self.axes1.set_xlim([0,1])
        self.axes1.set_ylim([0,10])
        pass


class ProfileCanvas(MplCanvas):
    """ Canvas to plot the growth profile """
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)


    def compute_initial_figure(self):
        pass

class SedCanvas(MplCanvas):
    """ Canvas to plot the SED """

    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)

    def compute_initial_figure(self):
        pass


# Style sheet
# Background color (FFF7C0 is buttermilk, DCAE1D is honey, F2D388 is butter)
# Colors from https://designschool.canva.com/blog/website-color-schemes/
# To work on MAC-OSX for QToolBar, one has to set the border to force the style on the system
# https://bugreports.qt.io/browse/QTBUG-12717
# Qt documentation: http://doc.qt.io/qt-5/stylesheet-examples.html

style = """
        QMainWindow {
        background-color: QLinearGradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 LemonChiffon, stop: 1 #F2D388);
        }
        QMenu {
        background-color: #D8AB4E;
        background: '#FFF6BA';
        color: 'black';
        }
        QMenuBar {
        background-color: QLinearGradient(x1:0, y1:0, x2:0, y2:1, stop:0 #FFF6BA, stop:1 #F2D388);
        background: #F2D388;
        color: 'black';
        }
        QMenuBar::item {
        background: transparent;
        spacing: 3px;
        padding: 1px 4px;
        border-radius: 4px;
        }
        QMenuBar::item:selected { /* when selected using mouse or keyboard */
        background: #FFF6BA;
        }
        QMenuBar::item:pressed {
        background: #DCAE1D;
        }
        QStatusBar {
        background-color: #FFF6BA;
        border: 1px solid black;
        border-radius: 3px;
        }
        QToolBar#tb1, QToolBar#tb2, QToolBar#tb {
        background-color: transparent;
        border: 1px transparent;
        }
        QToolBar::separator{
        background-color: transparent;
        }
        QToolButton:pressed {
        background-color: LemonChiffon;
        border-radius: 3px;
        }
        QToolButton:hover {
        background-color: LemonChiffon;
        border-radius: 3px;
        }
        QToolButton:focused {
        background-color: LemonChiffon;
        border-radius: 3px;
        }
        QToolButton:checked {
        background-color: LemonChiffon;
        border-radius: 3px;
        }
        QToolTip {
        border: 1px solid black;
        padding: 2px;
        border-radius: 3px;
        opacity: 200;
        background-color: LemonChiffon;
        }
        QTabBar::tab-bar{
        alignment: right;
        }
"""

    

class ApplicationWindow(QMainWindow):
    '''
    Layout and methods of the main window
    '''
            
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(Qt.WA_DeleteOnClose)
        # Get the path of the package
        path0, file0 = os.path.split(__file__)

        # Define style
        self.setStyleSheet(style)

        # Menu
        self.file_menu = self.menuBar().addMenu('&File')
        self.quit_program = QAction('Quit',self,shortcut='Ctrl+q',triggered=self.fileQuit)
        self.file_menu.addAction(self.quit_program)
        
        self.help_menu = self.menuBar().addMenu('&Help')
        self.about_code = QAction('About',self,shortcut='Ctrl+h',triggered=self.about)
        self.help_menu.addAction(self.about_code)

        # Especially for MAC OS (do not put menu on the top bar)
        self.menuBar().setNativeMenuBar(False)

        # Define main widget
        self.main_widget = QWidget(self)

        # Define sub widgets
        self.ic1 = ImageCanvas(self.main_widget, width=9, height=9, dpi=100)
        self.ic2 = ImageCanvas(self.main_widget, width=9, height=9, dpi=100)
        self.ic3 = ImageCanvas(self.main_widget, width=9, height=9, dpi=100)
        self.ic4 = ImageCanvas(self.main_widget, width=9, height=9, dpi=100)
        self.ic5 = ImageCanvas(self.main_widget, width=9, height=9, dpi=100)
        #self.mpl_toolbar1 = NavigationToolbar(self.ic1, self)
        #self.mpl_toolbar1.pan('on')
        #self.mpl_toolbar1.setObjectName('tb1')
        #self.tb2 = NavigationToolbar(self.ic2, self)
        #self.tb2.pan('on')
        #self.tb2.setObjectName('tb2')

        self.pc = ProfileCanvas(self.main_widget, width=2, height=4, dpi=100)
        self.sc = SedCanvas(self.main_widget, width=2, height=4, dpi=100)


        # Status Bar
        self.sb = QStatusBar()
        self.sb.showMessage("Welcome to ElApSED !", 10000)
        
        # Layout
        mainLayout = QHBoxLayout(self.main_widget)
        splitter1 = QSplitter(Qt.Horizontal)
        splitter2 = QSplitter(Qt.Vertical)

        # the image
        #imageWidget = QWidget()
        #imageWidget.layout = QVBoxLayout()

        tabs = QTabWidget()
        tabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        tabs.layout = QVBoxLayout()
        tabs.layout.setSpacing(0)
        tabs.layout.setContentsMargins(0,0,0,0)        
        #tabs.resize(450,450)

        tab1 = QWidget()
        tab2 = QWidget()
        tab3 = QWidget()
        tab4 = QWidget()
        tab5 = QWidget()
        
        tablist = [tab1,tab2,tab3,tab4,tab5]
        iclist  = [self.ic1,self.ic2,self.ic3,self.ic4,self.ic5]
        bands = ['u','g','r','i','z']

        for t,ic,b in zip(tablist,iclist,bands):            
            t.layout = QVBoxLayout()
            t.layout.addWidget(ic)
            tabs.addTab(t,b)

        #imageWidget.layout.addWidget(tabs)
        #imageWidget.layout.addWidget(self.sb)
        
        # imageLayout = QVBoxLayout()
        # imageWidget = QWidget()
        # imageWidget.setLayout(imageLayout)
        # imageLayout.addWidget(self.sb)

        sedLayout = QVBoxLayout()
        sedWidget = QWidget()
        sedWidget.setLayout(sedLayout)
        sedLayout.addWidget(self.sc)

        profileLayout = QVBoxLayout()
        profileWidget = QWidget()
        profileWidget.setLayout(profileLayout)
        profileLayout.addWidget(self.pc)

        
        splitter2.addWidget(profileWidget)
        splitter2.addWidget(sedWidget)
        splitter2.addWidget(self.sb)
        
        splitter1.addWidget(tabs)
        splitter1.addWidget(splitter2)

        mainLayout.addWidget(splitter1)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

    def fileQuit(self):
        self.close()
       
    def about(self):
        # Get path of the package
        path0,file0 = os.path.split(__file__)
        file=open(path0+"/copyright.txt","r")
        message=file.read()
        QMessageBox.about(self, "About", message)

""" Main code """
        
def main():
    app = QApplication(sys.argv)
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    aw = ApplicationWindow()
    aw.setGeometry(100, 100, width*0.9, width*0.35)
    progname = 'Elliptical Aperture photometry SED (ElApSED)'
    aw.setWindowTitle("%s" % progname)
    aw.show()
    app.exec_()
 
