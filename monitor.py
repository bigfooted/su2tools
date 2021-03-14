# monitor: displays and updates the residuals of the su2 history file in a window.
import sys
import os.path
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout,\
                            QHBoxLayout, QGridLayout, QSplitter, QFrame, QTextEdit, QFileDialog, QGroupBox,\
                            QButtonGroup, QLabel, QCheckBox,QLabel
from PyQt5 import QtCore
import pandas as pd
import matplotlib
from matplotlib import figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


# we first open the file and check if it is a legacy file format
# for legacy format we have to skip a couple of rows in the header
def skipRowsInHeader(filename):
    checkheader = open(filename, "r")
    line = checkheader.readline()
    if line[0:5] == "TITLE":
        #print("legacy file")
        #isLegacy=True
        skipRows = [0,2]
    else:
        #print("no legacy file")
        #isLegacy=False
        skipRows = []
    checkheader.close()    
    
    return skipRows

# ############################################################### #
# ############################################################### #
# ############################################################### #
class MatplotlibFigure(QWidget):
    """ setup of layout with canvas and toolbar for matplotlib figure
    """
    # this is the class for the actual figure containing the residual data
    
    # constructor
    def __init__(self):
        super().__init__()
        
        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        matplotlib.rcParams['lines.linewidth'] = 2
        matplotlib.rcParams['figure.edgecolor'] = 'blue'
        matplotlib.rcParams['figure.frameon'] = True


        #matplotlib.rcParams["figure.subplot.left"]=0.07
        #matplotlib.rcParams["figure.subplot.right"]=0.75
        #matplotlib.rcParams["figure.subplot.top"]=0.96

        matplotlib.rcParams["figure.subplot.left"]=0.10
        matplotlib.rcParams["figure.subplot.right"]=0.6
        matplotlib.rcParams["figure.subplot.top"]=0.9
        matplotlib.rcParams["figure.subplot.bottom"]=0.1
        
        matplotlib.rcParams["figure.subplot.wspace"]=0.1
        matplotlib.rcParams["figure.subplot.hspace"]=0.1

        # We create a vertical layout with the toolbar at the top and the canvas (that contains the actual figure)
        # at the bottom. We then group everything so we can create a bordered layout
        layout = QGridLayout()
        self.setLayout(layout)

        groupbox = QGroupBox("Convergence")
        layout.addWidget(groupbox)
        
        mpllayout = QVBoxLayout()
        groupbox.setLayout(mpllayout)
        mpllayout.addWidget(self.toolbar)
        mpllayout.addWidget(self.canvas)

        # original
        #mpllayout = QVBoxLayout(self)
        #mpllayout.addWidget(self.toolbar)
        #mpllayout.addWidget(self.canvas)
        

    def plot(self,listonoff,filename,normalize):
        # this will plot the actual residuals in a (x,y) plot inside the canvas
        if (filename != ""):  
            
            skipNrRows = skipRowsInHeader(filename)
            
            self.figure.clf()
            # pandas data 
            #self.dataframe = pd.read_csv(filename,skiprows=[skipNrRows])
            self.dataframe = pd.read_csv(filename,skiprows=skipNrRows)
            # get rid of quotation marks in the column names
            self.dataframe.columns = self.dataframe.columns.str.replace('"','')
            # remove legacy string
            #self.dataframe.columns = self.dataframe.columns.str.replace('VARIABLES=','')

            # get rid of spaces in the column names
            self.dataframe.columns = self.dataframe.columns.str.replace(' ','')
            
            # limit the columns to the ones containing the strings rms and Res
            self.dfrms = self.dataframe.filter(regex='rms|Res')
    
            x = [i for i in range(len(self.dfrms))]
            y=[]
            for c in range(len(self.dfrms.columns)):
                y.append(self.dfrms.iloc[:,c].tolist())
            
            #self.figure.tight_layout()
            #self.figure.adjust(right=0.5)
            # one row,3 columns, and the figure spans the first and second column
            #ax = self.figure.add_subplot(1,4,(1,3))
            ax = self.figure.add_subplot(111)
            #ax = self.figure.subplots(nrows=1,ncols=1)

            ax.set_xlabel('iterations')
            ax.set_ylabel('residuals')
            
            #ax.legend(self.dfrms.columns,bbox_to_anchor=(1.01,1),loc='upper left',prop={'size': 6,'weight':'bold'},borderaxespad=None)
            #ax.legend(self.dfrms.columns)
                
            plotlist=[]
            for l in range(len(listonoff)):
                if (listonoff[l]==True):
                    norm = 1.0
                    yn = y[l]
                    if normalize==True:
                        if len(x)>1:
                            if abs(y[l][1]) > 1e-6:
                                norm = y[l][1]
                                yn = [1+i-norm for i in y[l]]
                    
                    #print("x=",x[1],", y=",norm)

                    p, = ax.plot(x, yn,label=self.dfrms.columns[l])
                    #print("xrange=",ax.xaxis.get_data_interval())
                    #print("yrange=",ax.yaxis.get_data_interval())                    
                    plotlist.append(p)

            
            #ax.legend(self.dfrms.columns,bbox_to_anchor=(1.01,1),loc='upper left',prop={'size': 6,'weight':'bold'},borderaxespad=None)
            ax.legend(bbox_to_anchor=(1.01,1),loc='upper left',prop={'size': 6,'weight':'bold'},borderaxespad=None)
            self.figure.tight_layout()
            
            self.canvas.draw_idle()
            #ax.legend(handles=[p1, p2], title='title', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
            #ax.legend(handles=plotlist,bbox_to_anchor=(1.05, 1), loc='upper left')
            #ax.set_yscale('log')


# ############################################################### #
# ############################################################### #
# ############################################################### #
class MainApp(QMainWindow):
    """ Create the main window
    """
    def __init__(self):

        super().__init__()
        self.title = 'SU2 convergence monitor'
        
        # start at top left corner of the screen
        self.left = 0
        self.top = 0

        # set initial window height and width
        self.width = 1280
        self.height = 780

        # set the styles of fonts etc.
        self.setStyleSheet('font-size: 12pt')

        # user interface
        self.makeUI()

    def makeUI(self):
        # set window geometry
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # generate the tab layout
        self.table_widget = MyTableWidget()
        self.setCentralWidget(self.table_widget)

        self.show()

# ############################################################### #
# ############################################################### #
# ############################################################### #
class MyTableWidget(QWidget):
    """ setup of the tab containing all the buttons and data
    """

    # ################################### #
    def __init__(self):

        super().__init__()
        self.layout = QVBoxLayout(self)
        
        # initialization of the filename containing the residual data
        self.filename=""
        
        # initialize the tab screens
        self.tabs = QTabWidget()

        # add the tab screens to the tab widget
        self.tabs.resize(600, 600)

        # normalize the plot with the value at the first iteration
        self.normalize = False
        
        # make the layout of the tabs
        self.make_tab()

        # initialize the tabs
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        # setup a timer for plot update every 'setInterval'  ms
        self.timer = QtCore.QTimer()
        self.timer.setInterval(500)
        self.timer.timeout.connect(self.plot_data)
        self.timer.start()


        

    # ################################### #
    def delete_tab(self):    
        for i in range(len(self.chkbx)):           
            self.btn_layout.removeWidget(self.chkbx[i])
            self.chkbx[i].setParent(None)
            self.chkbx[i].deleteLater()
        self.btn_layout.update()
        
        
    # ################################### #
    # checkboxes for all the individual residual lines of the variables    
    def update_chkbx(self):
        self.chkbx = []

        if (self.filename != ""):  
            
            skipNrRows = skipRowsInHeader(self.filename)
            
            # read file into pandas dataframe, specific to su2 
            self.dataframe = pd.read_csv(self.filename,skiprows=skipNrRows)
            # get rid of quotation marks in the column names
            self.dataframe.columns = self.dataframe.columns.str.replace('"','')
            # remove legacy string
            #self.dataframe.columns = self.dataframe.columns.str.replace('VARIABLES=','')
            # get rid of spaces in the column names
            self.dataframe.columns = self.dataframe.columns.str.replace(' ','')
            # limit the columns to the ones containing the string 'rms'
            self.dfrms = self.dataframe.filter(regex='rms|Res')
            
            for name in self.dfrms:
                self.chkbx.append(QCheckBox(name))
            for i in range(len(self.chkbx)):    
                self.chkbx[i].toggle()
                self.chkbx[i].toggled.connect(self.changedata)
    
    # ################################### #
    def update_tab(self):    

        self.update_chkbx()

        '''Button layout section'''
        for i in range(len(self.chkbx)):           
            self.btn_layout.addWidget(self.chkbx[i])
        self.btn_layout.addStretch()
        self.btn_layout.update()
    
    # ################################### #
    # update plot based on normalization checkbox
    def update_normalize(self,state):
        if state == QtCore.Qt.Checked:
            self.normalize=True
            #print('Checked')
        else:
            self.normalize=False
            #print('Unchecked')
            
            
    
    # ################################### #
    def make_tab(self):
        
        self.figure = MatplotlibFigure()

        self.update_chkbx()
        
        '''file dialog section'''
        self.btn_openfile = QPushButton("Select File")
        self.btn_openfile.clicked.connect(self.getfile)
        self.show_filename = QLabel("No file selected")
        self.btn_openfile.setToolTip('click to select the file containing the residual data')
        
        '''button section'''
        # temporarily stop the realtime update
        self.btn_plot_data = QPushButton('stop')
        self.btn_plot_data.setToolTip('Click to stop real time data gathering')
        self.btn_plot_data.resize(50, 50) 
        self.btn_plot_data.clicked.connect(self.timer_startstop)

        # normalize checkbox    
        self.btn_normalize = QCheckBox('Normalize')
        self.btn_normalize.setToolTip('Normalize all residuals with the value at iteration 1 ')
        self.btn_normalize.stateChanged.connect(self.update_normalize)
        
        # now we add all the buttons in a vertical list
        # 
        self.btn_layout = QVBoxLayout()
        self.btn_layout.addWidget(self.btn_openfile)        
        self.btn_layout.addWidget(self.show_filename)
        self.btn_layout.addWidget(self.btn_plot_data)
        self.btn_layout.addWidget(self.btn_normalize)
        
        for i in range(len(self.chkbx)):           
            #print("i=",i)
            self.btn_layout.addWidget(self.chkbx[i])

        
        #self.btn_normalize.addWidget(self.btn_normalize)
            
        # this adds stretch, so the button is always top-aligned (in a vbox)
        self.btn_layout.addStretch()

        # those were the individual widgets, we now place the widgets in a frame and separate the frame from the plot using a splitter

        '''left-side layout'''
        left = QFrame()
        left.setLayout(self.btn_layout)


        # combine the buttons and the canvas in a splitter layout
        splitter1 = QSplitter(QtCore.Qt.Horizontal)
        splitter1.addWidget(left)
        splitter1.addWidget(self.figure)
        splitter1.setSizes([75, 700])

        '''master layout of tab'''
        self.tab = QWidget()
        
        # create tab with horizontal layout 
        self.tab.layout = QHBoxLayout()

        # add the last splitter to the layout
        self.tab.layout.addWidget(splitter1)
        self.tab.setLayout(self.tab.layout)
        
        self.tabs.addTab(self.tab, 'Convergence plot')

        checklist=[]
        for c in self.chkbx:
            checklist.append(c.isChecked())
            
        self.figure.plot(checklist,self.filename,self.normalize)       
        
        
    # ################################### #
    # file dialog for monitor input
    def getfile(self):

        # now we need to clear the list of widgets
        self.delete_tab()

        #fname = QFileDialog.getOpenFileName(self, 'Open file', 'history.csv',"csv files (*.txt *.csv)")
        fname = QFileDialog.getOpenFileName(self, 'Open file', 'history.csv',filter='csv files (*.csv)\n'
                                                                                    'dat files (*.dat)\n'
                                                                                    'all files (*)')
        self.filename=fname[0]
        # split the filename 
        fn = os.path.basename(self.filename)
        fb = os.path.basename(os.path.dirname(self.filename))
        self.show_filename.setText(os.path.join(fb,fn))
        self.show_filename.setToolTip(self.filename)

        # update tab
        self.update_tab()    
        
 
    # ################################### #        
    # start or stop the realtime update using the start/stop button
    def timer_startstop(self):
        if (self.timer.isActive()==True):
            #print("timer is active. it will be stopped")
            self.timer.stop()
            self.btn_plot_data.setText('start')
            self.btn_plot_data.setToolTip('Click to start real time data gathering')
            self.btn_plot_data.update()
            self.btn_layout.update()
        else:
            self.timer.start()
            #print("timer is inactive. it will be started")
            self.btn_plot_data.setText('stop')
            self.btn_plot_data.setToolTip('Click to stop real time data gathering')
            self.btn_plot_data.update()
            self.btn_layout.update()
            
    # ################################### #
    def plot_data(self):
        checklist=[]
        for c in self.chkbx:
            checklist.append(c.isChecked())
        self.figure.plot(checklist,self.filename,self.normalize)

    # ################################### #
    def changedata(self):
        checklist=[]
        for c in self.chkbx:
            checklist.append(c.isChecked())
        self.figure.plot(checklist,self.filename,self.normalize)

        
# ############################################################### #
# ############################################################### #
# ############################################################### #
if __name__ == '__main__':

    def run_app():
        app = QApplication(sys.argv)
        mainWin = MainApp()
        mainWin.show()
        app.exec_()

    run_app()
