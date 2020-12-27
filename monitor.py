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
        

    def plot(self,listonoff,filename):
        # this will plot the actual residuals in a (x,y) plot inside the canvas
        if (filename != ""):  
            self.figure.clf()
            # pandas data 
            self.dataframe = pd.read_csv(filename)
            # get rid of quotation marks in the column names
            self.dataframe.columns = self.dataframe.columns.str.replace('"','')
            # get rid of spaces in the column names
            self.dataframe.columns = self.dataframe.columns.str.replace(' ','')
            # limit the columns to the ones containing the string 'rms'
            self.dfrms = self.dataframe.filter(regex='rms')
    
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
                    p, = ax.plot(x, y[l],label=self.dfrms.columns[l])
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
        
        
    def update_chkbx(self):
        self.chkbx = []

        if (self.filename != ""):  
            # pandas data, specific to su2 
            self.dataframe = pd.read_csv(self.filename)
            # get rid of quotation marks in the column names
            self.dataframe.columns = self.dataframe.columns.str.replace('"','')
            # get rid of spaces in the column names
            self.dataframe.columns = self.dataframe.columns.str.replace(' ','')
            # limit the columns to the ones containing the string 'rms'
            self.dfrms = self.dataframe.filter(regex='rms')
            
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
    def make_tab(self):
        
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

         # setting size of button (width,height)
        self.btn_plot_data.resize(50, 50) 
        self.btn_plot_data.clicked.connect(self.timer_startstop)

        # this is a one-column list with the start/stop button at the top
        # followed by the list of lines
        self.btn_layout = QVBoxLayout()
        self.btn_layout.addWidget(self.btn_openfile)        
        self.btn_layout.addWidget(self.show_filename)
        self.btn_layout.addWidget(self.btn_plot_data)
        for i in range(len(self.chkbx)):           
            print("i=",i)
            self.btn_layout.addWidget(self.chkbx[i])

        # this adds stretch, so the button is always top-aligned (in a vbox)
        self.btn_layout.addStretch()


        '''left-side layout'''

        left = QFrame()
        left.setLayout(self.btn_layout)

        self.figure = MatplotlibFigure()

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
        self.figure.plot(checklist,self.filename)       
        
        
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
        self.figure.plot(checklist,self.filename)

    # ################################### #
    def changedata(self):
        checklist=[]
        for c in self.chkbx:
            checklist.append(c.isChecked())
        self.figure.plot(checklist,self.filename)

        
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
