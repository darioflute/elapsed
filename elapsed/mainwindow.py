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

        
    def compute_initial_figure(self, image=None, wcs=None, center=None):
        if wcs == None:
            pass
        else:
            self.wcs = wcs
            self.axes = self.fig.add_subplot(111, projection = wcs)
            #wcs.wcs.print_contents()
            # This works only if the x-axis is R.A. ?
            self.axes.coords[0].set_major_formatter('hh:mm:ss')
            self.image = self.axes.imshow(image, cmap='gist_heat_r', origin='lower', interpolation='none')
            self.axes.grid(color='black', ls='dashed')
            #self.axes.set_xlabel('R.A.')
            #self.axes.set_ylabel('Dec')
            # Mark center
            xc,yc = wcs.wcs_world2pix(center[0],center[1],1)
            self.axes.scatter(x=xc,y=yc,marker='+',s=400,c='green')
            # Colorbar
            cbaxes = self.fig.add_axes([0.9,0.1,0.02,0.8])
            self.fig.colorbar(self.image, cax=cbaxes)
            # Sliders to adjust intensity
            self.ax_cmin = self.figure.add_axes([0.1, 0.01, 0.8, 0.01])
            self.ax_cmax = self.figure.add_axes([0.1, 0.04, 0.8, 0.01])
            self.ax_cmin.clear()
            self.ax_cmax.clear()
            vmed0=np.nanmedian(image)
            d0 = np.nanstd(image)
            self.s_cmin = Slider(self.ax_cmin, 'low', vmed0-2*d0, vmed0+5*d0, valinit=vmed0-d0, facecolor='goldenrod')
            self.s_cmax = Slider(self.ax_cmax, 'high', vmed0-2*d0, vmed0+5*d0, valinit=vmed0+4*d0, facecolor='goldenrod')
            self.image.set_clim([vmed0-d0,vmed0+4*d0])
            self.s_cmin.valtext.set_visible(False)
            self.s_cmax.valtext.set_visible(False)
            self.slider1=self.s_cmin.on_changed(self.updateScale)
            self.slider2=self.s_cmax.on_changed(self.updateScale)

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
        with open(path0+'/yellow.stylesheet',"r") as fh:
            print ('reading stylesheet')
            self.setStyleSheet(fh.read())

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
        alphas = centers['ra'].values
        deltas = centers['dec'].values
        
        self.tabs = QTabWidget()
        self.tabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)

        # Check images for the source
        sources = centers['source'].values
        source = sources[0]

        path = './'+'{:d}'.format(source)+'/'
        files = sorted(glob.glob(path+'*.fits'))
        files = np.array(files)
        
        bands = [os.path.splitext(os.path.basename(f))[0] for f in files]

        # Reorder bands according to wavelength
        filters = pd.read_csv('filters.csv',  names=['filter','wvl','zp','f0','delta'],skiprows=1,
                              header=None, index_col='filter',
                              converters = {'filter': strip}, skipinitialspace=True
        )

        w = []
        newbands = bands.copy()
        for b in bands:
            try:
                w.append(filters.loc[b]['wvl'])
            except:
                print (b,' removed because is not in the filter list')
                newbands.remove(b)
        bands = np.array(newbands)
        bands = bands[np.argsort(w)].tolist()

        # Open tabs
        self.tabi = []
        self.ici = []
        for b in bands:
            t,ic = self.addImage(b)
            self.tabi.append(t)
            self.ici.append(ic)
            
        # Plot images
        for ima in bands:
            print ('image ', ima)
            image,wcs = self.readFits(path+ima+'.fits')
            self.ici[bands.index(ima)].compute_initial_figure(image=image,wcs=wcs,center=(alphas[0],deltas[0]))

            

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
        
        splitter1.addWidget(self.tabs)
        splitter1.addWidget(splitter2)

        mainLayout.addWidget(splitter1)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        

    def fileQuit(self):
        self.close()


    def addImage(self,b):
        ''' Add a tab with an image '''
        t = QWidget()
        t.layout = QVBoxLayout()
        self.tabs.addTab(t, b)
        ic = ImageCanvas(t, width=4, height=4, dpi=100)
        ic.toolbar = NavigationToolbar(ic, self)
        ic.toolbar.pan('on')
        t.layout.addWidget(ic)
        t.layout.addWidget(ic.toolbar)
        t.setLayout(t.layout)
        # connect image
        ic.mpl_connect('draw_event', self.onDraw)
        return t,ic

    def onDraw(self, event):
        ''' New limits for all images if zoom is activated '''
        itab = self.tabs.currentIndex()
        ic = self.ici[itab]
        if ic.toolbar._active == 'ZOOM':
            x = ic.axes.get_xlim()
            y = ic.axes.get_ylim()
            ra,dec = ic.wcs.wcs_pix2world(x,y,1)
            self.sb.showMessage("Resizing figures .... ", 1000)
            for ima in self.ici:
                x,y = ima.wcs.wcs_world2pix(ra,dec,1)
                ima.axes.set_xlim(x)
                ima.axes.set_ylim(y)
                ima.fig.canvas.draw_idle()
            self.sb.showMessage("Figures have been resized ", 1000)
            ic.toolbar.pan('on')

        
        
    def readFits(self, infile):
        import astropy.units as u
        import math
        ''' Read a fits file '''
        hdl = fits.open(infile)
        header = hdl[0].header
        image = hdl[0].data
        hdl.close()
        wcs = WCS(header)

        if hasattr(wcs.wcs,'cd'):
            CD11 = wcs.wcs.cd[0,0]
            CD12 = wcs.wcs.cd[0,1]
            CD21 = wcs.wcs.cd[1,0]
            CD22 = wcs.wcs.cd[1,1]
            if (abs(CD21) > abs(CD22)) and (CD21 >= 0): 
                North = "Right"
                positionAngle = 270.*u.deg + math.degrees(math.atan(CD22/CD21))*u.deg
            elif (abs(CD21) > abs(CD22)) and (CD21 < 0):
                North = "Left"
                positionAngle = 90.*u.deg + math.degrees(math.atan(CD22/CD21))*u.deg
            elif (abs(CD21) < abs(CD22)) and (CD22 >= 0):
                North = "Up"
                positionAngle = 0.*u.deg + math.degrees(math.atan(CD21/CD22))*u.deg
            elif (abs(CD21) < abs(CD22)) and (CD22 < 0):
                North = "Down"
                positionAngle = 180.*u.deg + math.degrees(math.atan(CD21/CD22))*u.deg
            if (abs(CD11) > abs(CD12)) and (CD11 > 0): East = "Right"
            if (abs(CD11) > abs(CD12)) and (CD11 < 0): East = "Left"
            if (abs(CD11) < abs(CD12)) and (CD12 > 0): East = "Up"
            if (abs(CD11) < abs(CD12)) and (CD12 < 0): East = "Down"
            if North == "Up" and East == "Left": imageFlipped = False
            if North == "Up" and East == "Right": imageFlipped = True
            if North == "Down" and East == "Left": imageFlipped = True
            if North == "Down" and East == "Right": imageFlipped = False
            if North == "Right" and East == "Up": imageFlipped = False
            if North == "Right" and East == "Down": imageFlipped = True
            if North == "Left" and East == "Up": imageFlipped = True
            if North == "Left" and East == "Down": imageFlipped = False
            print ('image flipped is ', imageFlipped)
        
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
 
