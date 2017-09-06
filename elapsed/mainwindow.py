#!/usr/bin/env python

# Essential libraries
import sys,os,glob
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# Matplotlib
import matplotlib
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1

from matplotlib.widgets import Slider

# FITS
from astropy.io import fits
from astropy.wcs import WCS


# Connection with PyQt5
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QWidget, QMainWindow, QMessageBox,QToolBar,QAction,QStatusBar,QSizePolicy,
                             QHBoxLayout, QVBoxLayout, QStackedLayout, QApplication, QListWidget,QSplitter,QMenu,QTabWidget)
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

    def compute_initial_figure(self, image=None, wcs=None):
        if image == None:
            pass
        else:
            self.axes = self.fig.add_subplot(111, projection = wcs)
            self.axes.coords[0].set_major_formatter('hh:mm:ss')
            self.image = self.axes.imshow(image, cmap='gist_heat_r', origin='lower', interpolation='none')
            # Colorbar
            cbaxes = self.fig.add_axes([0.9,0.1,0.02,0.8])
            self.fig.colorbar(self.image, cax=cbaxes)
            # Sliders to adjust intensity
            self.ax_cmin = self.figure.add_axes([0.1, 0.01, 0.8, 0.01])
            self.ax_cmax = self.figure.add_axes([0.1, 0.04, 0.8, 0.01])
            self.ax_cmin.clear()
            self.ax_cmax.clear()
            vmin0=np.nanmin(image); vmax0=np.nanmax(image)
            d0 = (vmax0-vmin0)/20.
            self.s_cmin = Slider(self.ax_cmin, 'low', vmin0-d0, vmax0+d0, valinit=vmin0, facecolor='goldenrod')
            self.s_cmax = Slider(self.ax_cmax, 'high', vmin0-d0, vmax0+d0, valinit=vmax0, facecolor='goldenrod')
            self.s_cmin.valtext.set_visible(False)
            self.s_cmax.valtext.set_visible(False)
            self.slider1=self.s_cmin.on_changed(self.updateScale)
            self.slider2=self.s_cmax.on_changed(self.updateScale)
        #self.axes.set_xlim([0,1])
        #self.axes.set_ylim([0,10])

    def updateScale(self,val):
        _cmin = self.s_cmin.val
        _cmax = self.s_cmax.val
        self.image.set_clim([_cmin, _cmax])
        self.fig.canvas.draw_idle()




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
        QToolBar {
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
        QTabWidget::tab-bar{
        alignment: left;
        }
        QTabWidget{
        background-color: transparent;
        }
"""

def strip(text):
    ''' strips whitespaces out of field'''
    try:
        return text.strip()
    except AttributeError:
        return text


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


        self.pc = ProfileCanvas(self.main_widget, width=3, height=1.5, dpi=100)
        self.sc = SedCanvas(self.main_widget, width=3, height=2, dpi=100)


        # Status Bar
        self.sb = QStatusBar()
        self.sb.showMessage("Welcome to ElApSED !", 10000)
        
        # Layout
        mainLayout = QHBoxLayout(self.main_widget)
        splitter1 = QSplitter(Qt.Horizontal)
        splitter2 = QSplitter(Qt.Vertical)

        # Check the file with centers
        centers = pd.read_csv('centers.csv',  names=['source','ra','dec'],skiprows=1, header=None)
        print ('source: ', centers.head())
        print ('RA: ', centers['ra'])
        ras = centers['ra'].values
        print("first RA is: ", ras[0])
        
        tabs = QTabWidget()
        tabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)

        # Check images for the source
        sources = centers['source'].values
        source = sources[0]

        path = './'+'{:d}'.format(source)+'/'
        files = sorted(glob.glob(path+'*.fits'))
        files = np.array(files)
        print ('files are: ', files)

        
        bands = np.array([os.path.splitext(os.path.basename(f))[0] for f in files])
        print ('bands are: ', bands)

        # Reorder bands according to wavelength
        filters = pd.read_csv('filters.csv',  names=['filter','wvl','zp','f0','delta'],skiprows=1,
                              header=None, index_col='filter',
                              converters = {'filter': strip}, skipinitialspace=True
        )

        w = np.array([filters.loc[b]['wvl'] for b in bands])
        bands = bands[np.argsort(w)].tolist()

        # Here I should check if other images (not included in the filter list )are present and discard them

        # Open tabs
        tabi = []
        ici = []
        for b in bands:
            t,ic = self.addImage(b, tabs)
            tabi.append(t)
            ici.append(ic)
            
        # Plot images
        for ima in bands:
            print ('image ', ima)
            image,wcs = self.readFits('/Users/dfadda/Python/A85/2/'+ima+'.fits')
            ici[bands.index(ima)].compute_initial_figure(image,wcs)

            

        #self.mpl_toolbar1 = NavigationToolbar(self.ic1, self)
        #self.mpl_toolbar1.pan('on')
        #self.mpl_toolbar1.setObjectName('tb1')

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


    def addImage(self,b, tabs):
        ''' Add a tab with an image '''
        t = QWidget()
        t.layout = QVBoxLayout()
        tabs.addTab(t, b)
        ic = ImageCanvas(t, width=4, height=4, dpi=100)
        ntb = NavigationToolbar(ic, self)
        ntb.pan('on')
        t.layout.addWidget(ic)
        t.layout.addWidget(ntb)
        t.setLayout(t.layout)
        return t,ic

    def readFits(self, infile):
        ''' Read a fits file '''
        hdl = fits.open(infile)
        header = hdl[0].header
        image = hdl[0].data
        hdl.close()
        wcs = WCS(header)
        return image, wcs
       
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
    aw.setGeometry(100, 100, width*0.7, width*0.45)
    progname = 'Elliptical Aperture photometry for Spectral Energy Distributions (ElApSED)'
    aw.setWindowTitle("%s" % progname)
    aw.show()
    app.exec_()
 
