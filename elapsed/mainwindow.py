#!/usr/bin/env python

# Essential libraries
import sys,os,glob
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
#import time

# Matplotlib
import matplotlib
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1

# FITS
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales as pixscales

# Connection with PyQt5
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QWidget, QMainWindow, QMessageBox,QToolBar,QAction,QStatusBar,QSizePolicy,
                             QHBoxLayout, QVBoxLayout, QApplication, QSplitter, QTabWidget)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt, QThread, QTimer#, QObject

from elapsed.canvas import ImageCanvas, ImageHistoCanvas, ProfileCanvas, SedCanvas
from elapsed.canvas import sourceDialog

def strip(text):
    ''' strips whitespaces out of field'''
    try:
        return text.strip()
    except AttributeError:
        return text


class AlignImagesThread(QThread):
    ''' Thread to align images '''
    
    def __init__(self, ic, ici, parent=None):
        super(AlignImagesThread, self).__init__(parent)
        self.ic = ic
        self.ici = ici
        
    def run(self):
        x = self.ic.axes.get_xlim()
        y = self.ic.axes.get_ylim()
        ra,dec = self.ic.wcsn.all_pix2world(x,y,1)
        for ima in self.ici:
            x,y = ima.wcsn.all_world2pix(ra,dec,1)
            ima.axes.set_xlim(x)
            ima.axes.set_ylim(y)
            ima.fig.canvas.draw_idle()

class ApplicationWindow(QMainWindow):
    '''
    Layout and methods of the main window
    '''
            
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(Qt.WA_DeleteOnClose)
        
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480
        
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Get the path of the package
        self.path0, file0 = os.path.split(__file__)

        # Define style
        with open(self.path0+'/yellow.stylesheet',"r") as fh:
            self.setStyleSheet(fh.read())
        
        self.createMenu()

        # Define main widget
        self.main_widget = QWidget(self)

        # Define main plot canvases
        self.pc = ProfileCanvas(self.main_widget, width=4.5, height=1, dpi=100)
        self.sc = SedCanvas(self.main_widget, width=4.5, height=3, dpi=100)

        # Status Bar
        self.sb = QStatusBar()
        self.sb.showMessage("Welcome to ElApSED !", 10000)
                
        # Tabs with images
        self.readCentersFilters()
        # Check the file with centers
        self.nSource = 0
        
        self.createToolbar()

        self.tabs = QTabWidget()
        self.tabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        self.tabs.currentChanged.connect(self.onChange)  # things to do when changing tab

        # Check images for the source (this will be put in a separate list later)
        # Select default source number for start
        self.createTabs()
        
        # Layout
        mainLayout = QHBoxLayout(self.main_widget)
        splitter1 = QSplitter(Qt.Horizontal)
        splitter2 = QSplitter(Qt.Vertical)

        sedLayout = QVBoxLayout()
        sedWidget = QWidget()
        sedWidget.setLayout(sedLayout)
        sedLayout.addWidget(self.sc)

        profileLayout = QVBoxLayout()
        profileWidget = QWidget()
        profileWidget.setLayout(profileLayout)
        profileLayout.addWidget(self.pc)

        
        plotLayout = QVBoxLayout()
        plotWidget = QWidget()
        plotWidget.setLayout(plotLayout)
        plotLayout.addWidget(sedWidget)
        plotLayout.addWidget(self.tb)
        plotLayout.addWidget(self.sb)
        
        
        splitter2.addWidget(profileWidget)
        splitter2.addWidget(plotWidget)
        
        splitter1.addWidget(self.tabs)
        splitter1.addWidget(splitter2)

        mainLayout.addWidget(splitter1)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        # Timer for periodical events
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.blinkTab)
        
    def createMenu(self):
        # Menu
        self.file_menu = self.menuBar().addMenu('&File')
        self.quit_program = QAction('Quit',self,shortcut='Ctrl+q',triggered=self.fileQuit)
        self.file_menu.addAction(self.quit_program)
        
        self.help_menu = self.menuBar().addMenu('&Help')
        self.about_code = QAction('About',self,shortcut='Ctrl+h',triggered=self.about)
        self.help_menu.addAction(self.about_code)

        # Especially for MAC OS (do not put menu on the top bar)
        self.menuBar().setNativeMenuBar(False)
        
    def createToolbar(self):
        # Actions
        alignAction = self.createAction(self.path0+'/icons/align.png','Align images','Ctrl+A',self.alignImages)
        self.blink = 'off'
        blinkAction = self.createAction(self.path0+'/icons/blink.png','Blink between 2 images','Ctrl+B',self.blinkImages)        
        levelsAction = self.createAction(self.path0+'/icons/levels.png','Adjust image levels','Ctrl+L',self.changeVisibility)        
        sourceAction = self.createAction(self.path0 + '/icons/importimage.png', 'Select a source file', 'Ctrl+S', self.updateSource)
        sources = self.centers['source'].values
        self.sourceList = sources
        # import this
        self.selectSource = sourceDialog(self.sourceList, self.nSource)
        
        # Commented out to prevent 'instantaneous' changing of images. Too laggy.
        # self.selectSource.slist.currentRowChanged.connect(self.updateSource)
        
        # Toolbar
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.setObjectName('toolbar')
        self.tb.addAction(alignAction)
        self.tb.addAction(blinkAction)
        self.tb.addAction(levelsAction)
        self.tb.addAction(sourceAction)
     
    def updateSource(self, event):
        sourceDialog(self.sourceList, self.nSource)
        self.selectSource.exec_()
        if (self.selectSource.launchTabs == True) :
            self.createTabs()
        else :
            pass

    def readCentersFilters(self):
        self.centers = pd.read_csv('centers.csv',  names=['source','ra','dec'],skiprows=1, header=None)
        self.filters = pd.read_csv('filters.csv',  names=['filter','wvl','zp','f0','delta'],skiprows=1,
                              header=None, index_col='filter',
                              converters = {'filter': strip}, skipinitialspace=True
        )
    
    def createTabs(self):
        
        # removes any tabs that already exist from display; does not delete them.
        self.tabs.clear()
        
        sources = self.centers['source'].values
        source = sources[self.nSource]
        
        try:
            source = '{:d}'.format(source)
        except:
            source = strip(source)
            
        #path = './'+'{:d}'.format(source)+'/'
        path = './'+source+'/'
        files = sorted(glob.glob(path+'*.fits'))
        files = np.array(files)
        
        bands = [os.path.splitext(os.path.basename(f))[0] for f in files]
        w = []
        newbands = bands.copy()
        for b in bands:
            try:
                w.append(self.filters.loc[b]['wvl'])
            except:
                print (b,' removed because is not in the filter list')
                newbands.remove(b)
        bands = np.array(newbands)
        self.bands = bands[np.argsort(w)].tolist()
        
        self.openTabs()
            
        # Plot images
        alphas = self.centers['ra'].values
        deltas = self.centers['dec'].values
        for ima in self.bands:
             print ('image ', ima)
             image,wcs = self.readFits(path+ima+'.fits')
             ic = self.ici[self.bands.index(ima)]
             ic.compute_initial_figure(image=image,wcs=wcs,
                                       center=(alphas[self.nSource],deltas[self.nSource]),title=source)
             # Callback to propagate axes limit changes among images
             ic.cid = ic.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.zoomAll)
             ih = self.ihi[self.bands.index(ima)]
             clim = ic.image.get_clim()
             ih.compute_initial_figure(image=image,xmin=clim[0],xmax=clim[1])
            
    def openTabs(self):
        # Open tabs
        self.tabi = []
        self.ici = []
        self.ihi = []
        for b in self.bands:
            t,ic,ih = self.addImage(b)
            self.tabi.append(t)
            self.ici.append(ic)
            self.ihi.append(ih)

    def changeVisibility(self):
        itab = self.tabs.currentIndex()
        ih = self.ihi[itab]
        state = ih.isVisible()
        ih.setVisible(not state)

    def onChange(self, itab):
        ''' When tab changes check if latest update of ellipse are implemented '''
        #print("current index is ", itab)
        if itab < len(self.ici):
            ima = self.ici[itab]
            if ima.changed:
                canvas = ima.arcell[0].figure.canvas
                canvas.draw_idle()
                ima.changed = False
            if self.blink == 'select':
                # Select 2nd tab and start blinking until blink status changes ...
                self.btab[1] = itab
                self.blink = 'on'
                self.timer.start(1000)

    def blinkTab(self):
        ''' keep switching between two tabs until blink changes state '''
        itab = self.tabs.currentIndex()
        if itab == self.btab[0]:
            i = 1
        else:
            i = 0
        self.tabs.setCurrentIndex(self.btab[i])
        
    def blinkImages(self, event):
        ''' Blink between two images in different tabs or stop blinking'''
        if self.blink == 'off':
            self.btab = [self.tabs.currentIndex(),0]
            self.sb.showMessage("Select another tab to blink / click again to stop blinking", 2000)
            self.blink = 'select'
        else:
            self.blink = 'off'
            self.timer.stop()
                        
    def fileQuit(self):
        self.close()

    def createAction(self,icon,text,shortcut,action):
        act = QAction(QIcon(icon),text, self)
        act.setShortcut(shortcut)
        act.triggered.connect(action)
        return act

        
    def addImage(self,b):
        ''' Add a tab with an image '''
        t = QWidget()
        t.layout = QVBoxLayout()
        self.tabs.addTab(t, b)
        ic = ImageCanvas(t, width=6, height=5.5, dpi=100)
        ih = ImageHistoCanvas(t, width=6, height=0.25, dpi=100)
        ih.mySignal.connect(self.onChangeIntensity)
        #ih.setVisible(False)
        ic.toolbar = NavigationToolbar(ic, self)
        #ic.toolbar.pan('on')
        t.layout.addWidget(ic)
        t.layout.addWidget(ih)
        t.layout.addWidget(ic.toolbar)
        t.setLayout(t.layout)
        # connect image to draw events
        ic.mpl_connect('button_release_event', self.onDraw)
        #ih.mpl_connect('button_release_event', self.onChangeIntensity)
        ic.mpl_connect('scroll_event',self.onWheel)
        return t,ic,ih

    def onChangeIntensity(self, event):
        itab = self.tabs.currentIndex()
        ic = self.ici[itab]
        ih = self.ihi[itab]
        # apply intensity limits to the relative figure
        ic.image.set_clim(ih.limits)
        ic.fig.canvas.draw_idle()

    
    def onDraw(self, event):
        ''' propagate ellipse changes to other figures '''
        itab = self.tabs.currentIndex()
        ic = self.ici[itab]
        if ic.drrEllipse.changed:
            ic.drrEllipse.changed=False
            x,y = ic.arcell[0].center
            ra,dec = ic.wcsn.all_pix2world(x,y,1)
            delta = pixscales(ic.wcsn)
            w = ic.arcell[0].width*delta[0] # assuming pixels are square
            h = ic.arcell[0].height*delta[0]
            a = ic.arcell[0].angle
            if ic.flip:
                a = -a
            ici = self.ici.copy()
            ici.remove(ic)
            for ima in ici:
                x,y = ima.wcsn.all_world2pix(ra,dec,1)
                delta = pixscales(ima.wcsn)
                w_ = w/delta[0]
                h_ = h/delta[0]
                if ima.flip:
                    a_ = -a
                else:
                    a_ = a
                for arc in ima.arcell:
                    arc.width = w_
                    arc.height= h_
                    arc.angle = a_
                    arc.center = (x,y)
                ima.changed = True # To redraw next time tab is open

    def zoomAll(self, event):
        ''' propagate limit changes to all images '''
        itab = self.tabs.currentIndex()
        ic = self.ici[itab]
        if ic.axes == event: # only consider axes on screen (not other tabs)
            if ic.toolbar._active == 'ZOOM':
                ic.toolbar.zoom()  # turn off zoom
            x = ic.axes.get_xlim()
            y = ic.axes.get_ylim()
            #print (x,y)
            ra,dec = ic.wcsn.all_pix2world(x,y,1)
            #self.sb.showMessage("Resizing figures .... ", 1000)
            ici = self.ici.copy()
            ici.remove(ic)
            #print(ra,dec)
            for ima in ici:
                x,y = ima.wcsn.all_world2pix(ra,dec,1)
                ima.axes.set_xlim(x)
                ima.axes.set_ylim(y)
                ima.changed = True

    def onWheel(self,event):
        ''' enable zoom with mouse wheel and propagate changes to other tabs '''
        eb = event.button
        itab = self.tabs.currentIndex()
        ic = self.ici[itab]

        curr_xlim = ic.axes.get_xlim()
        curr_ylim = ic.axes.get_ylim()
        curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
        curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
        if eb == 'up':
            factor=0.9
        elif eb == 'down':
            factor=1.1
        new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
        new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
        x = [curr_x0-new_width,curr_x0+new_width]
        y = [curr_y0-new_height,curr_y0+new_height]
        ic.axes.set_xlim(x)
        ic.axes.set_ylim(y)
        ic.fig.canvas.draw_idle()
        ici = self.ici.copy()
        ici.remove(ic)
        ra,dec = ic.wcsn.all_pix2world(x,y,1)
        for ima in ici:
            x,y = ima.wcsn.all_world2pix(ra,dec,1)
            ima.axes.set_xlim(x)
            ima.axes.set_ylim(y)
            ima.changed = True

            
    def alignImages(self, event):
        ''' Align all images using the same RA-Dec limits as current image '''

        itab = self.tabs.currentIndex()
        ic = self.ici[itab]
        if ic.toolbar._active == 'ZOOM':
            ic.toolbar.zoom() # turn off zoom
        ici = self.ici.copy()
        ici.remove(ic)
        x = ic.axes.get_xlim()
        y = ic.axes.get_ylim()
        ra,dec = ic.wcsn.all_pix2world(x,y,1)
        for ima in self.ici:
            x,y = ima.wcsn.all_world2pix(ra,dec,1)
            ima.axes.set_xlim(x)
            ima.axes.set_ylim(y)
            ima.changed = True

    def readFits(self, infile):
        ''' Read a fits file '''
        hdl = fits.open(infile)
        header = hdl[0].header
        image = hdl[0].data
        hdl.close()
        wcs = WCS(header)
        
        # Center the reference pixel
        nx, ny = header["naxis1"], header["naxis2"]
        i0, j0 = (float(nx) + 1.0)/2, (float(ny) + 1.0)/2
        [ra0, dec0], = wcs.all_pix2world([[i0, j0]], 1)
        header.update(crpix1=i0, crpix2=j0, crval1=ra0, crval2=dec0)
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
    aw.setGeometry(width*0.025, 0, width*0.95, width*0.5)
    progname = 'Elliptical Aperture photometry for Spectral Energy Distributions (ElApSED)'
    aw.setWindowTitle("%s" % progname)
    aw.show()
    app.exec_()
 
