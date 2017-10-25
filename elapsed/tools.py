import numpy as np

class DragResizeRotateEllipse:
    #from PyQt5.QtCore import pyqtSignal    
    #updateEllipses = pyqtSignal('QString')

    lock = None
    def __init__(self, arcell, border_tol=0.6, changed=False):
        self.arcell = arcell
        self.ellipse = arcell[0]
        self.arc1 = arcell[1]
        self.arc2 = arcell[2]
        self.arc3 = arcell[3]
        self.arc4 = arcell[4]
        for arc in arcell:
            arc.set_animated(True)
        self.border_tol = border_tol
        self.lock = None
        self.press = None
        self.background = None
        self.changed = changed
        self.axes = arcell[0].axes
        canvas = self.axes.figure.canvas        
        # connect to events
        self.ciddraw = canvas.mpl_connect('draw_event', self.on_draw)
        self.cidpress = canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas = canvas
        self.spot = 0
        
    def on_press(self, event):
        '''on button press it stores some data if mouse is over it'''
        if event.inaxes != self.ellipse.axes: return
        if DragResizeRotateEllipse.lock is not None: return
        contains, attrd = self.ellipse.contains(event)
        if not contains: return

        x0, y0 = self.ellipse.center
        w0, h0 = self.ellipse.width, self.ellipse.height
        theta0 = self.ellipse.angle
        self.lock = "pressed"
        self.press = x0, y0, w0, h0, theta0, event.xdata, event.ydata
        DragResizeRotateEllipse.lock = self
        # draw everything
        self.canvas.draw_idle()
        # register changes
        self.changed = True

    def on_draw(self, event):
        ''' Drawing the ellipse '''
        # first capture the background
        self.background = self.canvas.copy_from_bbox(self.axes.bbox)
        # then draw the ellipse
        for arc in self.arcell:
            self.axes.draw_artist(arc)

    def on_motion(self, event):
        '''on motion it will act on the ellipse if the mouse is over it'''
        if DragResizeRotateEllipse.lock is not self: return
        if not self.ellipse.contains(event): return

        # store original position, compute coordinates increment, update ellipse
        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        self.dx = event.xdata - xpress
        self.dy = event.ydata - ypress
        self.update_ellipse()

        # restore the background region
        # redraw just the current ellipse
        self.canvas.restore_region(self.background)
        for arc in self.arcell:
            self.axes.draw_artist(arc)
        # blit just the redrawn area
        #self.canvas.blit(self.axes.bbox)
        self.canvas.update()
        self.canvas.flush_events()

        
    def on_release(self, event):
        '''on release it resets the press data'''
        if DragResizeRotateEllipse.lock is not self:
            return

        self.press = None
        DragResizeRotateEllipse.lock = None
        self.lock = "released"
        self.spot = 0

        # turn off the animation property and reset the background
        self.background = None
        
    def disconnect(self):
        'disconnect all the stored connection ids'
        self.canvas.mpl_disconnect(self.ciddraw)
        self.canvas.mpl_disconnect(self.cidpress)
        self.canvas.mpl_disconnect(self.cidrelease)
        self.canvas.mpl_disconnect(self.cidmotion)

    def update_ellipse(self):
        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        dx, dy = self.dx, self.dy
        bt = self.border_tol
        # Normalized point (to the circle)
        xnorm, ynorm = self.ellipse.get_patch_transform().inverted().transform_point((xpress, ypress))

        # lock into a mode
        if self.lock == "pressed":
            print('Pressed !')
            rnorm = np.sqrt(xnorm*xnorm+ynorm*ynorm)
            if rnorm > bt:
                anorm = np.arctan2(ynorm,xnorm)*180./np.pi
                self.lock = "resizerotate"
                th0 = theta0/180.*np.pi
                c, s = np.cos(th0), np.sin(th0)
                R = np.matrix('{} {}; {} {}'.format(c, s, -s, c))
                (dx_,dy_), = np.array(np.dot(R,np.array([dx,dy])))
                if abs(dx_) > 1.2*abs(dy_) and (abs(anorm) < 30.):
                    self.spot = 1
                elif abs(dx_) > 1.2*abs(dy_) and (abs(anorm) > 150.):
                    self.spot = 3
                elif abs(dy_) > 1.2*abs(dx_) and (anorm > 60. and anorm < 120.):
                    self.spot = 2
                elif abs(dy_) > 1.2*abs(dx_) and (anorm < -60. and anorm > -120.):
                    self.spot = 4
                else:
                    self.spot = 0
            else:
                self.lock = "move"
                
        elif self.lock == "move":
            print('Move !!!')
            xn = x0+dx; yn =y0+dy
            if xn < 0: xn = x0
            if yn < 0: yn = y0
            for arc in self.arcell:
                arc.center = (xn,yn)
        elif self.lock == "resizerotate":
            print('Resize rotate !!!')
            dtheta = np.arctan2(ypress+dy-y0,xpress+dx-x0)-np.arctan2(ypress-y0,xpress-x0)
            dtheta *= 180./np.pi
            theta_ = theta0+dtheta

            anorm = np.arctan2(ynorm,xnorm)*180./np.pi
            th0 = theta0/180.*np.pi
            c, s = np.cos(th0), np.sin(th0)
            R = np.matrix('{} {}; {} {}'.format(c, s, -s, c))
            (dx_,dy_), = np.array(np.dot(R,np.array([dx,dy])))

            if self.spot == 1:
                #print('arc1')
                w_ = w0+2*dx_  if (w0+2*dx_) > 0 else w0 # Avoid flipping
                h_ = h0
            elif self.spot ==3:
                #print('arc3')
                w_ = w0-2*dx_  if (w0-2*dx_) > 0 else w0
                h_ = h0
            elif self.spot == 2:
                #print('arc2')
                h_ = h0+2*dy_  if (h0+2*dx_) > 0 else h0
                w_ = w0
            elif self.spot == 4:
                #print('arc4')
                h_ = h0-2*dy_  if (h0-2*dx_) > 0 else h0
                w_ = w0
            else:
                self.lock = "released"

            if self.lock != "released":
                for arc in self.arcell:
                    arc.width = w_
                    arc.height = h_
                    arc.angle = theta_

from matplotlib.patches import Rectangle
class RectangleSelect:
    ''' Add rectangle to the figure to crop the spectral cube'''
    def __init__(self, parent):
        self.frame = parent
        self.rect = Rectangle((0,0), 0.1, 0.1, facecolor='None', edgecolor='None')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.rectpatch = self.frame.axes.add_patch(self.rect)
        self.pressed = False
        self.connect()
        
    def connect(self):
        ''' connect to events '''
        self.cidpress = self.rect.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def disconnect(self):
        '''disconnect all the stored connection ids'''
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)
        self.rectpatch.remove()
        self.frame.canvas.draw_idle()
        
    def on_press(self, event):
        ''' Callback to handle the mouse being clicked and held over the canvas'''
        # Check the mouse press was actually on the canvas 
        if event.xdata is not None and event.ydata is not None:
            # Upon initial press of the mouse record the origin and record the mouse as pressed
            self.pressed = True
            self.rect.set_linestyle('dashed')
            self.rect.set_edgecolor('blue')
            self.x0 = event.xdata
            self.y0 = event.ydata
            print ("pressed ",self.x0,self.y0)
            
    def on_motion(self, event):
        '''Callback to handle the motion event created by the mouse moving over the canvas'''

        # If the mouse has been pressed draw an updated rectangle when the mouse is moved so 
        # the user can see what the current selection is
        if self.pressed:
            # Check the mouse was released on the canvas, and if it wasn't then just leave the width and 
            # height as the last values set by the motion event
            if event.xdata is not None and event.ydata is not None:
                self.x1 = event.xdata
                self.y1 = event.ydata
                
            # Set the width and height and draw the rectangle
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.frame.canvas.draw()

    def get_xlim(self):
        if self.x1 > self.x0:
            return (self.x0,self.x1)
        else:
            return (self.x1,self.x0)
        
    def get_ylim(self):
        if self.y1 > self.y0:
            return (self.y0,self.y1)
        else:
            return (self.y1,self.y0)

    def get_rect(self):
        return self.rect
 
    def on_release(self, event):
        '''Callback to handle the mouse being released over the canvas'''
        
        # Check that the mouse was actually pressed on the canvas to begin with and this isn't a rouge mouse 
        # release event that started somewhere else
        if self.pressed:
            # Upon release draw the rectangle as a solid rectangle
            self.pressed = False
            self.rect.set_linestyle('solid')

            # Check the mouse was released on the canvas, and if it wasn't then just leave the width and 
            # height as the last values set by the motion event
            if event.xdata is not None and event.ydata is not None:
                self.x1 = event.xdata
                self.y1 = event.ydata

            # Set the width and height and origin of the bounding rectangle
            self.boundingRectWidth =  self.x1 - self.x0
            self.boundingRectHeight =  self.y1 - self.y0
            self.bouningRectOrigin = (self.x0, self.y0)

            # Draw the bounding rectangle
            self.rect.set_width(self.boundingRectWidth)
            self.rect.set_height(self.boundingRectHeight)
            self.rect.set_xy((self.x0, self.y0))
            self.frame.canvas.draw()

