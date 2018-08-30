import matplotlib
import numpy as np
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QSizePolicy
from PyQt5.QtCore import QSize, pyqtSignal#, QObject

# FITS
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales as pixscales

#from matplotlib.widgets import Slider 
from matplotlib.widgets import SpanSelector
# from matplotlib.patches import Ellipse,Arc
# from elapsed.apertures import EllipseInteractor
# from matplotlib.widgets import EllipseSelector


class MplCanvas(FigureCanvas):
    """ Basic matplotlib canvas class """

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        #FigureCanvas.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
        FigureCanvas.setSizePolicy(self,QSizePolicy.MinimumExpanding,QSizePolicy.MinimumExpanding)
        FigureCanvas.updateGeometry(self)
        self.compute_initial_figure()

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(1,1)
    
    def compute_initial_figure(self):
        pass

    
class ImageCanvas(MplCanvas):
    """ Canvas to plot an image """
    
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)

            
    def compute_initial_figure(self, image=None, wcs=None, center=None, title=None):
        if wcs == None:
            ''' initial definition when images are not yet read '''
            pass
        else:
            h1 = wcs.to_header()
            #if hasattr(wcs.wcs,'cd'):
            try:
                pc11=h1["PC1_1"]
                pc12=h1["PC1_2"]
                pc21=h1["PC2_1"]
                pc22=h1["PC2_2"]
                pc1 = -np.hypot(pc11,pc12)
                if h1["CRVAL2"] < 0:    # if object in Southern Emisphere
                    print('crval2 < 0')
                    pc2 = -np.hypot(pc21,pc22)
                else:
                    pc2 = np.hypot(pc21,pc22)
                print ('pc2 is ', pc2)
                h1.update(
                    pc1_1=pc1,
                    pc1_2=0.0, 
                    pc2_1=0.0, 
                    pc2_2=pc2,
                    orientat=0.0)
            except:
                pass
            self.wcsn = WCS(h1)
            self.wcs = wcs

            self.axes = self.fig.add_axes([0.1,0.1,.8,.8], projection = self.wcsn)
            self.image = self.axes.imshow(image, cmap='gist_heat_r',interpolation='none',
                                          transform=self.axes.get_transform(self.wcs))
            self.axes.coords[0].set_major_formatter('hh:mm:ss')
            self.axes.grid(color='black', ls='dashed')
            self.axes.set_xlabel('R.A.')
            self.axes.set_ylabel('Dec')

            # Plot with North up (corners are clockwise from left-bottom)
            corners = self.wcsn.calc_footprint()
            self.flip = False
            if corners[0,1] < 0:
                if corners[0,1] > corners[1,1]:
                    ylim = self.axes.get_ylim()
                    self.axes.set_ylim([ylim[1],ylim[0]])
                    self.flip = True

            # Add title
            self.fig.suptitle('Source: '+title)

            # Colorbar
            cbaxes = self.fig.add_axes([0.9,0.1,0.02,0.85])
            self.fig.colorbar(self.image, cax=cbaxes)

            # Sliders to adjust intensity
            #self.ax_cmin = self.figure.add_axes([0.1, 0.01, 0.8, 0.01])
            #self.ax_cmax = self.figure.add_axes([0.1, 0.04, 0.8, 0.01])
            #self.ax_cmin.clear()
            #self.ax_cmax.clear()
            vmed0=np.nanmedian(image)
            d0 = np.nanstd(image)
            #self.s_cmin = Slider(self.ax_cmin, 'low', vmed0-2*d0, vmed0+5*d0, valinit=vmed0-d0, facecolor='goldenrod')
            #self.s_cmax = Slider(self.ax_cmax, 'high', vmed0-2*d0, vmed0+5*d0, valinit=vmed0+4*d0, facecolor='goldenrod')
            self.image.set_clim([vmed0-d0,vmed0+4*d0])
            #self.s_cmin.valtext.set_visible(False)
            #self.s_cmax.valtext.set_visible(False)
            #self.slider1=self.s_cmin.on_changed(self.updateScale)
            #self.slider2=self.s_cmax.on_changed(self.updateScale)

            # Mark center
            xc,yc = self.wcsn.all_world2pix(center[0],center[1],1)
            self.axes.scatter(x=xc,y=yc,marker='+',s=400,c='green')
            self.pixscale = pixscales(self.wcsn)[0]*3600.
            self.changed = False
            # Draw canvas
            #canvas = self.axes.figure.canvas
            #canvas.draw()            
            self.apertures=[]

    def updateScale(self,val):
        _cmin = self.s_cmin.val
        _cmax = self.s_cmax.val
        self.image.set_clim([_cmin, _cmax])
        self.fig.canvas.draw_idle()


class ImageHistoCanvas(MplCanvas):
    """ Canvas to plot the histogram of image intensity """
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)

    mySignal = pyqtSignal(str)
        
    def compute_initial_figure(self, image=None,xmin=None,xmax=None):
        if image is None:
            ''' initial definition when images are not yet read '''
            pass
        else:
            self.axes = self.fig.add_axes([0.0,0.4,1.,1.])
            self.axes.yaxis.set_major_formatter(plt.NullFormatter())
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            self.axes.spines['left'].set_visible(False)
            # Print the histogram of finite values
            ima = image.ravel()
            mask = np.isfinite(ima)
            ima = ima[mask]
            ima = np.sort(ima)
            s = np.size(ima)
            smax = min(int(s*0.9995),s-1)
            nbins=256
            n, self.bins, patches = self.axes.hist(ima, bins=nbins, 
                                                   range=(np.nanmin(ima), ima[smax]), fc='k', ec='k')
            # Define the interval containing 99% of the values

            self.x = np.arange(s)
            if xmin == None:
                xmin = ima[int(s*0.01)]
            if xmax == None:
                xmax = ima[int(s*0.99)-1]
            self.onSelect(xmin,xmax)
            # Start a span selector
            self.span = SpanSelector(self.axes, self.onSelect, 'horizontal', useblit=True,
                                     rectprops=dict(alpha=0.5, facecolor='LightSalmon'),button=1)
            # Define the way to draw/read a shaded interval (check from sospex)
            # Communicate values to the image (from the tab we know the canvas, but I have to communicate changes to the main window)


    def onSelect(self,xmin, xmax):
        indmin, indmax = np.searchsorted(self.bins, (xmin, xmax))
        indmax = min(len(self.bins) - 1, indmax)
        self.limits = [self.bins[indmin],self.bins[indmax]]
        try:
            self.shade.remove()
        except:
            pass
        self.mySignal.emit('limits changed')
        self.shade = self.axes.axvspan(self.limits[0],self.limits[1],facecolor='Lavender',
                                       alpha=0.5,linewidth=0)
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

from PyQt5.QtWidgets import QDialog, QLabel, QListWidget, QListWidgetItem, QPushButton, QVBoxLayout
import os
from PyQt5.QtGui import QIcon

class sourceDialog(QDialog):

    dirSignal = pyqtSignal(str)
    
    def __init__(self,stlist,  currentST, parent=None):
        super().__init__()
        # self.currentRow = currentST
        path0, file0 = os.path.split(__file__)
        self.setWindowTitle('Title')
        layout = QVBoxLayout()
        self.launchTabs = False
        self.list = stlist
        label = QLabel("Source Files")        
        self.slist = QListWidget(self)
        self.slist.setSizePolicy(QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding))
        self.slist.setMaximumSize(QSize(160,150))
        for st in stlist:
            st = str(st)
            QListWidgetItem(QIcon(path0+"/icons/"+st+"_.png"),st,self.slist)
        # n = stlist.index(currentST)
        # self.slist.setCurrentRow(self.currentRow)
        
        # Button with OK to close dialog
        b2 = QPushButton("OK",self)
        b2.clicked.connect(self.accept)

        # Layout
        layout.addWidget(label)
        layout.addWidget(self.slist)
        layout.addWidget(b2)
        self.setLayout(layout)

    def end(self):
        self.close()
    
    def accept(self):
        self.launchTabs = True
        self.end()