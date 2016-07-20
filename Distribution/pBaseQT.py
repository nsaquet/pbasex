# -*- coding: utf-8 -*-
# Created: Mon Sep 30 10:02:05 2013
#      by: pyside-uic 0.2.13 running on PySide 1.1.1
#
# GUI interface for pBase program with a Qt interface
# Use some of the function build in pBaseCore 
# Use  matplotlib for displaying images

import os
from PySide import QtCore, QtGui
import matplotlib
matplotlib.use('Qt4Agg')
import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg
from matplotlib.figure import Figure
from matplotlib.pyplot import figure,setp
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.ticker import NullFormatter

from pBaseCore import Datas,theta_f,symmetrize

waitCondition = QtCore.QWaitCondition()
mutex = QtCore.QMutex()

class pBaseForm(QtGui.QMainWindow):
    def setupUi(self, MainWindow):
        #Define datas
        self.plotsettings=PlotSettings()
        self.workflow=Datas()
        self.file_path=''
        self.file_name=''
        self.path_to_Basis=""
        #Create Main Window
        MainWindow.setObjectName("MainWindow")
        
        #Fit Main window with a Widget
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        #Create the frame for the graph and add it to the main layout
        self.Wplot = QtGui.QWidget()
        self.LPlot= QtGui.QVBoxLayout()
        self.LPlot.setContentsMargins(5, 5, 5, 5)
        self.fig = Figure((6., 6.6), dpi=150,frameon=False)
        self.canvas = FigCanvas(self.fig)
        self.canvas.setParent(self.Wplot)
        self.canvas.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        self.gs=gridspec.GridSpec(2,1,height_ratios=[4,1],hspace=0.1,bottom=0.05,top=0.995)
        self.axes=self.fig.add_subplot(self.gs[0])
        self.axesPES=self.fig.add_subplot(self.gs[1])
        self.gs.tight_layout(self.fig,pad=0.1)
        self.Wplot.setObjectName("Image")
        #mpl_toolbar = NavigationToolbar(self.canvas, self.Wplot)
        #self.LPlot.addWidget(mpl_toolbar)
        self.LPlot.addWidget(self.canvas)
        self.Wplot.setLayout(self.LPlot)
        
        # Bind the 'pick' event for clicking on one of the bars
        self.canvas.mpl_connect('button_press_event', self.on_press)
        
        #Create the Tool Box on the side
        self.VerticalWidget = QtGui.QWidget()
        self.VerticalWidget.setFixedWidth(350)
        self.VerticalBox = QtGui.QVBoxLayout()
        self.VerticalBox.setContentsMargins(5, 5, 5, 5)
        self.VerticalBox.setObjectName("verticalLayout")
        #Design the external Group box
        self.groupBox = QtGui.QGroupBox(self.VerticalWidget)
        self.groupBox.setObjectName("ToolBox")
        
        self.verticalLayout = QtGui.QVBoxLayout(self.groupBox)
        self.verticalLayout.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout.setObjectName("verticalLayout")
        
        #Design a choice in display
        self.DisplayCLayout = QtGui.QHBoxLayout()
        self.DisplayCLayout.setContentsMargins(5, 5, 5, 5)
        self.DisplayCLayout.setObjectName("DisplayChoiceLayout")
        self.label22 = QtGui.QLabel(self.groupBox)
        self.label22.setObjectName("label22")
        self.DisplayCLayout.addWidget(self.label22)
        self.DisplayChoiceBox = QtGui.QComboBox(self.groupBox)
        self.DisplayChoiceBox.setObjectName("DisplayChoiceBox")
        self.DisplayChoiceBox.addItem("Original")
        self.DisplayChoiceBox.addItem("Inverted")
        self.DisplayCLayout.addWidget(self.DisplayChoiceBox)
        self.DisplayCLayout.insertSpacing(2,15) 
        self.verticalLayout.addLayout(self.DisplayCLayout)  
        self.DisplayChoiceBox.activated[str].connect(self.changeDisplay) 
        
        
        #Design the polynom group box
        self.PolyGBox = QtGui.QGroupBox()
        self.PolyGBox.setFlat(True)
        self.PolyGBox.setObjectName("PolynomeGBox")
        self.PolyLayout = QtGui.QHBoxLayout(self.PolyGBox)
        self.PolyLayout.setSpacing(5)
        self.PolyLayout.setContentsMargins(5, 5, 5, 5)
        self.PolyLayout.setObjectName("PolynomeBoxLayout")
        self.label = QtGui.QLabel(self.PolyGBox)
        self.label.setObjectName("label")
        self.PolyLayout.addWidget(self.label)
        self.PolyBox = QtGui.QSpinBox(self.PolyGBox)
        self.PolyBox.setWrapping(False)
        self.PolyBox.setFrame(True)
        self.PolyBox.setButtonSymbols(QtGui.QAbstractSpinBox.PlusMinus)
        self.PolyBox.setSpecialValueText("")
        self.PolyBox.setAccelerated(False)
        self.PolyBox.setSuffix("")
        self.PolyBox.setMinimum(0)
        self.PolyBox.setMaximum(100)
        self.PolyBox.setProperty("value", 2)
        self.PolyBox.setObjectName("PolyBox")
        self.PolyBox.valueChanged.connect(self.ChangeLmax)
        self.PolyLayout.addWidget(self.PolyBox)
        self.PolyLayout.insertSpacing(2,15)

        self.OddBox = QtGui.QCheckBox(self.PolyGBox)
        self.OddBox.setObjectName("OddBox")
        self.OddBox.stateChanged.connect(self.ChangeOdd)
        self.PolyLayout.addWidget(self.OddBox)
        self.PolyGBox.setLayout(self.PolyLayout)
        self.verticalLayout.addWidget(self.PolyGBox)
        
                
        self.DataChoiceGBox = QtGui.QGroupBox()
        self.DataChoiceGBox.setFlat(True)
        self.DataChoiceGBox.setObjectName("DataChoiceBox")
        self.DCLayout = QtGui.QHBoxLayout(self.DataChoiceGBox)
        self.DCLayout.setContentsMargins(5, 5, 5, 5)
        self.DCLayout.setObjectName("DataChoiceLayout")
        self.label2 = QtGui.QLabel(self.DataChoiceGBox)
        self.label2.setObjectName("label2")
        self.DCLayout.addWidget(self.label2)
        self.DataChoiceBox = QtGui.QComboBox(self.DataChoiceGBox)
        self.DataChoiceBox.setObjectName("DataChoiceBox")
        self.DataChoiceBox.addItem("")
        self.DataChoiceBox.addItem("")
        self.DCLayout.addWidget(self.DataChoiceBox)
        self.DCLayout.insertSpacing(2,15) 
        self.DataChoiceGBox.setLayout(self.DCLayout)
        self.verticalLayout.addWidget(self.DataChoiceGBox)   
        
        #Set ellipticity
        self.EllipticityGBox = QtGui.QGroupBox()
        self.EllipticityGBox.setFlat(True)
        self.EllipticityGBox.setObjectName("EllipticityBox")
        self.ElLayout = QtGui.QHBoxLayout(self.EllipticityGBox)
        self.ElLayout.setContentsMargins(5, 5, 5, 5)
        self.ElLayout.setObjectName("EllipticityLayout")
        self.ESlider = QtGui.QSlider(self.EllipticityGBox)
        self.ESlider.setOrientation(QtCore.Qt.Horizontal)
        self.ESlider.setTickInterval(5)
        self.ESlider.setTickPosition(QtGui.QSlider.TicksBelow)
        self.ESlider.setValue(100)
        self.ESlider.setPageStep(1)
        self.ESlider.setRange(90,110)
        self.ESlider.valueChanged[int].connect(self.changeEValue)
        self.ESlider.setObjectName("ESlider")
        self.ElLayout.addWidget(self.ESlider)
        self.ElLayout.insertSpacing(2,15) 
        self.EllipticityGBox.setLayout(self.ElLayout)
        self.verticalLayout.addWidget(self.EllipticityGBox) 
        
        #Set invert option
        self.InvertOptGBox = QtGui.QGroupBox()
        self.InvertOptGBox.setFlat(True)
        self.InvertOptGBox.setObjectName("InvertOptionBox")
        self.IOptLayout = QtGui.QHBoxLayout(self.InvertOptGBox)
        self.IOptLayout.setContentsMargins(5, 5, 5, 5)
        self.IOptLayout.setObjectName("InvertOptionLayout")
        self.InvHalfBox = QtGui.QCheckBox(self.InvertOptGBox)
        self.InvHalfBox.stateChanged.connect(self.ChangeHalfInvert)
        self.InvHalfBox.setObjectName("InvHalfBox")
        self.IOptLayout.addWidget(self.InvHalfBox)
        self.WhichHalfBox = QtGui.QComboBox(self.InvertOptGBox)
        self.WhichHalfBox.setObjectName("WhichHalfBox")
        self.WhichHalfBox.addItem("Left")
        self.WhichHalfBox.addItem("Right")
        self.IOptLayout.addWidget(self.WhichHalfBox)
        self.IOptLayout.insertSpacing(2,15) 
        self.InvertOptGBox.setLayout(self.IOptLayout)
        self.verticalLayout.addWidget(self.InvertOptGBox)   
        
        #Define the Centering Tool Box
        self.CenterBox = QtGui.QGroupBox(self.groupBox)
        self.CenterBox.setFlat(True)
        self.CenterBox.setObjectName("CenterBox")
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setContentsMargins(5, 5, 5, 5)
        self.gridLayout.setVerticalSpacing(8)
        self.gridLayout.setHorizontalSpacing(8)
        self.gridLayout.setObjectName("CenterBoxLayout")
        self.XSlider = QtGui.QSlider(self.CenterBox)
        self.XSlider.setOrientation(QtCore.Qt.Horizontal)
        self.XSlider.setPageStep(1)
        self.XSlider.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.XSlider.valueChanged[int].connect(self.changeXValue)
        self.XSlider.setObjectName("XSlider")
        self.gridLayout.addWidget(self.XSlider, 0, 1, 1, 1)
        self.YSlider = QtGui.QSlider(self.CenterBox)
        self.YSlider.setOrientation(QtCore.Qt.Horizontal)
        self.YSlider.setPageStep(1)
        self.YSlider.setFocusPolicy(QtCore.Qt.StrongFocus)   
        self.YSlider.valueChanged[int].connect(self.changeYValue)
        self.YSlider.setObjectName("YSlider")
        self.gridLayout.addWidget(self.YSlider, 1, 1, 1, 1)
        self.TransposeButton = QtGui.QPushButton(self.CenterBox)
        self.TransposeButton.setObjectName("TransposeButton")
        self.TransposeButton.clicked.connect(self.TransposeFn)
        self.gridLayout.addWidget(self.TransposeButton, 2, 0, 1, 1)
        
        self.FixCenterBox = QtGui.QCheckBox(self.CenterBox)
        self.FixCenterBox.stateChanged.connect(self.ChangeFixCenter)
        self.FixCenterBox.setObjectName("FixCenterBox")
        self.gridLayout.addWidget(self.FixCenterBox, 2, 1, 1, 1)
        self.AutoButton = QtGui.QPushButton(self.CenterBox)
        self.AutoButton.setFlat(False)
        self.AutoButton.clicked.connect(self.AutoCenterFn)
        self.AutoButton.setObjectName("AutoButton")
        self.gridLayout.addWidget(self.AutoButton, 0, 0, 1, 1)
        self.SymButton = QtGui.QPushButton(self.CenterBox)
        self.SymButton.setObjectName("SymButton")
        self.SymButton.clicked.connect(self.SymmetrizeFn)
        self.gridLayout.addWidget(self.SymButton, 1, 0, 1, 1)
        self.ImageCenterBtn = QtGui.QPushButton(self.CenterBox)
        self.ImageCenterBtn.setObjectName("ImageCenterBtn")
        self.ImageCenterBtn.clicked.connect(self.ICenterFn)
        self.gridLayout.addWidget(self.ImageCenterBtn, 3, 0, 1, 1)
        self.gridLayout.setColumnStretch(0, 3)
        self.gridLayout.setColumnStretch(1, 4)
        self.CenterBox.setLayout(self.gridLayout)
        self.verticalLayout.addWidget(self.CenterBox)
        
        #Define the plot toolbox
        self.ColorBox = QtGui.QGroupBox(self.groupBox)
        self.ColorBox.setFlat(True)
        self.ColorBox.setCheckable(False)
        self.ColorBox.setObjectName("ColorBox")
        self.ColorLayout = QtGui.QHBoxLayout()
        self.ColorLayout.setContentsMargins(5, 5, 5, 5)
        self.ColorLayout.setObjectName("ColorLayout")
        self.ColorMapBox = QtGui.QComboBox(self.ColorBox)
        self.ColorMapBox.setObjectName("ColorMapBox")
        self.ColorMapBox.addItem("")
        self.ColorMapBox.addItem("")
        self.ColorMapBox.addItem("")
        self.ColorMapBox.addItem("")
        self.ColorMapBox.addItem("")
        self.ColorMapBox.addItem("")
        self.ColorMapBox.activated[str].connect(self.OnChooseCM)
        self.ColorLayout.addWidget(self.ColorMapBox)
        self.InvColorBox = QtGui.QCheckBox(self.ColorBox)
        self.InvColorBox.setObjectName("InvColorBox")
        self.InvColorBox.stateChanged.connect(self.InvertCM)
        self.ColorLayout.addWidget(self.InvColorBox)
        self.SqrtColorBox = QtGui.QCheckBox(self.ColorBox)
        self.SqrtColorBox.setObjectName("SqrtColorBox")
        self.SqrtColorBox.stateChanged.connect(self.SqrtCM)
        self.ColorLayout.addWidget(self.SqrtColorBox)
        self.ColorLayout.setStretch(0, 3)
        self.ColorLayout.setStretch(1, 2)
        self.ColorLayout.setStretch(2, 1)
        self.ColorBox.setLayout(self.ColorLayout)
        self.verticalLayout.addWidget(self.ColorBox)
        self.verticalLayout.setStretch(0, 1)
        self.verticalLayout.setStretch(1, 1)
        self.verticalLayout.setStretch(2, 1)
        self.verticalLayout.setStretch(3, 1.8)
        self.verticalLayout.setStretch(4, 1)
        
        #Define the Main Tool Box: Invert|Save|Close buttons
        self.HBox = QtGui.QHBoxLayout()
        self.HBox.setContentsMargins(2, 2, 2, 2)
        self.HBox.setObjectName("HBox")
        self.InvertButton = QtGui.QPushButton(self.groupBox)
        font = QtGui.QFont()
        font.setWeight(75)
        font.setStrikeOut(False)
        font.setBold(True)
        self.InvertButton.setFont(font)
        self.InvertButton.setObjectName("InvertButton")
        self.InvertButton.clicked.connect(self.InvertFn)
        self.HBox.addWidget(self.InvertButton)
        self.SaveButton = QtGui.QPushButton(self.groupBox)
        self.SaveButton.setObjectName("SaveButton")
        self.SaveButton.clicked.connect(self.openSave)
        self.HBox.addWidget(self.SaveButton)
        self.CloseButton = QtGui.QPushButton(self.groupBox)
        self.CloseButton.setObjectName("CloseButton")
        self.HBox.addWidget(self.CloseButton)
        
        #Add the layout to the central widget and then to the window
        self.VerticalBox.addWidget(self.groupBox)
        self.VerticalBox.addLayout(self.HBox)
        self.VerticalWidget.setLayout(self.VerticalBox)
        
        #Design the main Layout Frame | Control box: Horizontal
        self.MLayout = QtGui.QHBoxLayout()
        self.MLayout.setSpacing(5)
        self.MLayout.setContentsMargins(5, 5, 5, 5)
        self.MLayout.setObjectName("HorizontalBox")
        self.MLayout.addWidget(self.Wplot)
        self.MLayout.addWidget(self.VerticalWidget)
        self.MLayout.setAlignment(self.Wplot,QtCore.Qt.AlignCenter)
        
        #Define the Main window menu bar and statusbar
        self.menubar = QtGui.QMenuBar()
        self.menubar.setGeometry(QtCore.QRect(0, 0, 830, 22))
        self.menubar.setObjectName("menubar")
        self.menuFiles = QtGui.QMenu(self.menubar)
        self.menuFiles.setObjectName("menuFiles")
        MainWindow.setMenuBar(self.menubar)
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setObjectName("actionOpen")
        self.actionOpen.setStatusTip("Open new File")
        self.actionSave = QtGui.QAction(MainWindow)
        self.actionSave.setObjectName("actionSave")
        self.actionSave.setStatusTip("Save current Analyse")
        self.menuFiles.addAction(self.actionOpen)
        self.menuFiles.addAction(self.actionSave)
        self.menubar.addAction(self.menuFiles.menuAction())
        #Status bar
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.statlabel = QtGui.QLabel()
        self.statlabel.setText("Ready !")
        self.statusbar.addWidget(self.statlabel,1)
        self.statcoordinates = QtGui.QLabel()
        self.statcoordinates.setText("Center: x= , y= ")
        self.statusbar.addWidget(self.statcoordinates,2)
        self.progressBar = QtGui.QProgressBar()
        self.statusbar.addWidget(self.progressBar,2)
        self.progressBar.setValue(0)
        self.progressBar.setRange(0,20)
        self.progressBar.setVisible(False)

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.CloseButton, QtCore.SIGNAL("released()"), MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        self.actionOpen.triggered.connect(self.openFile)
        self.actionSave.triggered.connect(self.openSave)
        
        self.centralwidget.setLayout(self.MLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.display()
		
    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "pBaseQt", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setTitle(QtGui.QApplication.translate("MainWindow", "Treatments", None, QtGui.QApplication.UnicodeUTF8))
        self.PolyGBox.setTitle(QtGui.QApplication.translate("MainWindow", "Polynoms", None, QtGui.QApplication.UnicodeUTF8))
        self.DataChoiceGBox.setTitle(QtGui.QApplication.translate("MainWindow", "Choice of Data for the inversion", None, QtGui.QApplication.UnicodeUTF8))
        self.EllipticityGBox.setTitle(QtGui.QApplication.translate("MainWindow", "Image Ellipticity", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Max Legendre Polynoms:", None, QtGui.QApplication.UnicodeUTF8))
        self.label2.setText(QtGui.QApplication.translate("MainWindow", "Datas to be inverted:", None, QtGui.QApplication.UnicodeUTF8))
        self.label22.setText(QtGui.QApplication.translate("MainWindow", "Datas displayed:", None, QtGui.QApplication.UnicodeUTF8))
        self.OddBox.setText(QtGui.QApplication.translate("MainWindow", "Odd", None, QtGui.QApplication.UnicodeUTF8))
        self.DataChoiceBox.setItemText(0, QtGui.QApplication.translate("MainWindow", "Original", None, QtGui.QApplication.UnicodeUTF8))
        self.DataChoiceBox.setItemText(1, QtGui.QApplication.translate("MainWindow", "Current", None, QtGui.QApplication.UnicodeUTF8))
        self.CenterBox.setTitle(QtGui.QApplication.translate("MainWindow", "Centering", None, QtGui.QApplication.UnicodeUTF8))
        self.TransposeButton.setText(QtGui.QApplication.translate("MainWindow", "Transpose", None, QtGui.QApplication.UnicodeUTF8))
        self.FixCenterBox.setText(QtGui.QApplication.translate("MainWindow", "Fixed Center", None, QtGui.QApplication.UnicodeUTF8))
        self.AutoButton.setText(QtGui.QApplication.translate("MainWindow", "Auto center", None, QtGui.QApplication.UnicodeUTF8))
        self.SymButton.setText(QtGui.QApplication.translate("MainWindow", "Symmetrize V", None, QtGui.QApplication.UnicodeUTF8))
        self.ImageCenterBtn.setText(QtGui.QApplication.translate("MainWindow", "Center of the image", None, QtGui.QApplication.UnicodeUTF8))
        self.ColorBox.setTitle(QtGui.QApplication.translate("MainWindow", "Colormaps", None, QtGui.QApplication.UnicodeUTF8))
        self.ColorMapBox.setItemText(0, QtGui.QApplication.translate("MainWindow", "Igor", None, QtGui.QApplication.UnicodeUTF8))
        self.ColorMapBox.setItemText(1, QtGui.QApplication.translate("MainWindow", "Jet", None, QtGui.QApplication.UnicodeUTF8))
        self.ColorMapBox.setItemText(2, QtGui.QApplication.translate("MainWindow", "Gnuplot", None, QtGui.QApplication.UnicodeUTF8))
        self.ColorMapBox.setItemText(3, QtGui.QApplication.translate("MainWindow", "Gray", None, QtGui.QApplication.UnicodeUTF8))
        self.ColorMapBox.setItemText(4, QtGui.QApplication.translate("MainWindow", "Hot", None, QtGui.QApplication.UnicodeUTF8))
        self.ColorMapBox.setItemText(5, QtGui.QApplication.translate("MainWindow", "Inferno", None, QtGui.QApplication.UnicodeUTF8))        
        self.InvColorBox.setText(QtGui.QApplication.translate("MainWindow", "Inverse Colors", None, QtGui.QApplication.UnicodeUTF8))
        self.SqrtColorBox.setText(QtGui.QApplication.translate("MainWindow", "Sqrt", None, QtGui.QApplication.UnicodeUTF8))
        self.InvertButton.setText(QtGui.QApplication.translate("MainWindow", "Invert !", None, QtGui.QApplication.UnicodeUTF8))
        self.SaveButton.setText(QtGui.QApplication.translate("MainWindow", "Save", None, QtGui.QApplication.UnicodeUTF8))
        self.CloseButton.setText(QtGui.QApplication.translate("MainWindow", "Close", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFiles.setTitle(QtGui.QApplication.translate("MainWindow", "Files", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setText(QtGui.QApplication.translate("MainWindow", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+O", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave.setText(QtGui.QApplication.translate("MainWindow", "Save", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+S", None, QtGui.QApplication.UnicodeUTF8))
        self.InvertOptGBox.setTitle(QtGui.QApplication.translate("MainWindow", "Invert Options", None, QtGui.QApplication.UnicodeUTF8))
        self.InvHalfBox.setText(QtGui.QApplication.translate("MainWindow", "1/2 image only", None, QtGui.QApplication.UnicodeUTF8))
        
    def openFile(self):
        fname, _ = QtGui.QFileDialog.getOpenFileName(self,self.tr("Open data file"),"~/",self.tr("Fit Files (*.fit *.fits);;Image Files (*.tiff *.jpg *.bmp);;Txt (*.txt *.dat)"))
        if fname:
        	file=QtCore.QFileInfo(fname)
        	self.file_path=file.path()
        	self.file_name=file.fileName()
        	self.statlabel.setText("Opening File %s" %fname)
        	self.workflow.OpenFile(fname)
        	self.statlabel.setText("Go !")
        	self.workflow.lmax=self.PolyBox.value()
        	self.workflow.odd=self.OddBox.isChecked()
        else: self.statlabel.setText("Failed to open File")
        self.display()
        
    def openSave(self):
    	
    	msgbox=QtGui.QMessageBox(self)
    	msgbox.setText("Which file format do you wish to use?")
    	fitsbutton=msgbox.addButton("FITS",QtGui.QMessageBox.ActionRole)
    	datbutton=msgbox.addButton("Dat",QtGui.QMessageBox.ActionRole)
    	pdfbutton=msgbox.addButton("PDF",QtGui.QMessageBox.ActionRole)
    	msgbox.setDefaultButton(fitsbutton)
    	msgbox.exec_()
    	if msgbox.clickedButton()==fitsbutton:
    		outputname=self.file_name[:-4]+'_output.fit'
    		suggestedname=os.path.join(self.file_path,outputname)
        	fname, _ = QtGui.QFileDialog.getSaveFileName(self,self.tr("Save data file"),suggestedname,self.tr("Fits Files (*.fit)"))
        	if fname:
        		self.statlabel.setText("Saving File %s" %fname)
        		self.workflow.SaveFileFits(fname)
        	else: self.statlabel.setText("Failed to save File")
    	elif msgbox.clickedButton()==datbutton:
    		if self.file_name[-5:]=='.fits': outputname=self.file_name[:-5]+'_output.dat'
    		else: outputname=self.file_name[:-4]+'_output.dat'
    		suggestedname=os.path.join(self.file_path,outputname)
        	fname, _ = QtGui.QFileDialog.getSaveFileName(self,self.tr("Save data file"),suggestedname,self.tr("Dat Files (*.dat)"))
        	if fname:
        		self.statlabel.setText("Saving File %s" %fname)
        		self.workflow.SaveFileDat(fname)
        	else: self.statlabel.setText("Failed to save File")
    	elif msgbox.clickedButton()==pdfbutton:
    		outputname=self.file_name[:-4]+'_output.dat'
    		suggestedname=os.path.join(self.file_path,outputname)
        	fname, _ = QtGui.QFileDialog.getSaveFileName(self,self.tr("Save data file"),suggestedname,self.tr("PDF Files (*.pdf)"))
        	if fname:
        		fname_pes=fname[:-4]+'_PES.pdf'
        		fname_img=fname[:-4]+'_ImgInv.pdf'
        		self.statlabel.setText("Saving File %s" %fname_img)
        		
        		fig=figure()
        		axes=fig.add_subplot(111)
        		axes.plot(self.workflow.normed_pes,'k')
        		axes.set_yticks([0,0.5,1.])
        		axes.set_xlim([0,self.workflow.r])
        		fig.savefig(fname_pes,format='pdf')
        		del fig,axes
        		
        		fig=figure(figsize=(6.,6.),dpi=150)
        		axes=fig.add_subplot(1,1,1)
        		palette=self.plotsettings.cmapdic[self.plotsettings.palettename]
        		if self.plotsettings.IsSqrt and self.plotsettings.IsR: palette=self.plotsettings.cmapdic_r_sqrt[self.plotsettings.palettename]	
        		axes.imshow(self.workflow.datas,origin='lower',cmap=palette,vmax=0.8*self.workflow.datas.max())
        		nullfmt = NullFormatter()
        		axes.yaxis.set_major_formatter(nullfmt)
        		axes.xaxis.set_major_formatter(nullfmt)
        		fig.savefig(fname_img,format='pdf')
        		del fig,axes,palette
        		
        	else: self.statlabel.setText("Failed to save File")
        
    def display(self):
        """
            Deal with displaying the picture
        """
        self.axes.clear()
        self.axesPES.clear()
        #Deal with image
        xmax,ymax=self.workflow.datas.shape
        self.XSlider.setMaximum(xmax)
        self.YSlider.setMaximum(ymax)
        self.XSlider.setValue(self.workflow.center[0])
        self.YSlider.setValue(self.workflow.center[1])
        
        palette=self.plotsettings.cmapdic[self.plotsettings.palettename]
        if self.plotsettings.IsSqrt and self.plotsettings.IsR: palette=self.plotsettings.cmapdic_r_sqrt[self.plotsettings.palettename]	
        elif self.plotsettings.IsSqrt: palette=self.plotsettings.cmapdic_sqrt[self.plotsettings.palettename]
        elif self.plotsettings.IsR: palette=self.plotsettings.cmapdic_r[self.plotsettings.palettename]
        
        self.axes.imshow(self.workflow.datas,extent=[0,ymax,0,xmax],origin='lower',cmap=palette,vmax=0.8*self.workflow.datas.max())
        
        if self.workflow.r!=0.:
            xc,yc= paint_circle(self.workflow.center,self.workflow.r,self.workflow.scale.ellipticity)
            scalarMap = cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=256), cmap=palette)
            self.axes.plot(xc,yc,color=scalarMap.to_rgba(256),lw=1)
        
        self.axes.axis([0,ymax,0,xmax])
        # no label
        nullfmt = NullFormatter()
        self.axes.yaxis.set_major_formatter(nullfmt)
        self.axes.xaxis.set_major_formatter(nullfmt)
        self.axes.tick_params(bottom='off',top='off',left='off',right='off')
        #Deal with PES
        self.axesPES.plot(self.workflow.radial,self.workflow.normed_pes,'k')
        self.axesPES.set_yticks([0,0.5,1.])
        self.axesPES.set_xlim([0,self.workflow.r])
        setp(self.axesPES.get_xticklabels(),fontsize=10)
        setp(self.axesPES.get_yticklabels(),fontsize=10)
        self.axesPES.patch.set_alpha(0.1)
        #Deal with display
        self.canvas.draw()
        del palette
        self.statcoordinates.setText(u"Center: x={0[0]} , y={0[1]}.\t Radius: {1} \t Ellipticity: {2}".format(self.workflow.center,int(self.workflow.r),self.workflow.scale.ellipticity))
        
    def OnChooseCM(self,text):
        """
            Change cmap settings
        """
        self.plotsettings.palettename=text
        self.display()
        
    def InvertCM(self,state):
        """
            Inverse cmap scale 
        """
        if state == QtCore.Qt.Checked: self.plotsettings.IsR=True
        else: self.plotsettings.IsR=False
        self.display()
    
    def changeDisplay(self,state):
    	if state=='Inverted':
    		self.workflow.datas=self.workflow.output
    	else:
    		self.workflow.datas=self.workflow.raw
    	self.display()
        
    def SqrtCM(self,state):
        """
            Sqrt cmap scale 
        """
        if (state == QtCore.Qt.Checked): self.plotsettings.IsSqrt=True
        else: self.plotsettings.IsSqrt=False
        self.display() 

    def ChangeOdd(self,state):
        if state == QtCore.Qt.Checked: self.workflow.odd=1
        else: self.workflow.odd=0
    
    def ChangeHalfInvert(self,state):
        if state == QtCore.Qt.Checked: self.workflow.half=1
        else: self.workflow.half=0
    
    def ChangeLmax(self,number):
        self.workflow.lmax=number
    
    def ChangeFixCenter(self,state):
    	"""
    		If pressed, if fixes the centre of the image to the most recent value and 
    		prevents the user from changing it by clicking the mouse on the image window. 
    		It also forces the panel to give the distance from the centre instead of the 
    		cartesian coordinates.
    	"""
    	if state == QtCore.Qt.Checked: self.plotsettings.IsFixed=True
    	else: self.plotsettings.IsFixed=False
    
    def TransposeFn(self):
        self.workflow.datas=self.workflow.datas.T
        self.workflow.center=(self.workflow.center[1],self.workflow.center[0])
        self.statlabel.setText("Transposed")
        self.display()
        
    def SymmetrizeFn(self):
    	"""
    		Symmetrise a 2_D circular selection vertically (ie about a horizontal axis). 
    		Assume that the centre is mid-pixel (x0,y0) rather than at lower left corner
    		of pixel x0,y0. Note that no symmetrisation is needed 
    		horizontally since the Legendre polynomials are already 
    		symmetric along the vertical axis. (Vertically being the polarisation 
    		axis of the ligth, or the direction of propagation in the case of cpl).
    	"""
        self.workflow.datas=symmetrize(self.workflow.datas,self.workflow.center,self.workflow.r)
        self.statlabel.setText("Symmetrized")
        self.display()
    
    def AutoCenterFn(self):
        if self.workflow.r==0.: self.display()
        else: 
            self.progressBar.reset()
            self.progressBar.setVisible(True)
            if not self.plotsettings.IsFixed: 
                self.statlabel.setText("Centering")
                self.RunProcess()
                while not self.process.isFinished(): QtCore.QCoreApplication.processEvents()
            del self.process 
            self.statlabel.setText("Centered")
            self.progressBar.setVisible(False)
            self.display()
            
    def RunProcess(self):
        self.process=CenterProcesser(self)
        QtCore.QObject.connect(self.process,QtCore.SIGNAL("progress(int)"),self.progressBar,QtCore.SLOT("setValue(int)"))
        if not self.process.isRunning():
            self.process.exiting = False
            self.process.start() 
                
    def ICenterFn(self):
        if not self.plotsettings.IsFixed: self.workflow.get_com()
        self.display()
    
    def InvertFn(self):
    	#Path to the basis file
    	if not self.path_to_Basis:
    		self.path_to_Basis= QtGui.QFileDialog.getExistingDirectory(self,self.tr("Basis Files Directory"),dir="~/")
    	self.progressBar.reset()
        self.progressBar.setVisible(True)
    	self.process=InvertProcesser(self,self.path_to_Basis)
        QtCore.QObject.connect(self.process,QtCore.SIGNAL("progress(int)"),self.progressBar,QtCore.SLOT("setValue(int)"))
        if not self.process.isRunning():
            self.process.exiting = False
            self.process.start() 
        while not self.process.isFinished(): QtCore.QCoreApplication.processEvents()
        del self.process
        self.progressBar.setVisible(False)
    	self.workflow.datas=self.workflow.output
    	self.DisplayChoiceBox.setCurrentIndex(1)
    	self.display()
    
    def changeXValue(self,value):
    	if not self.plotsettings.IsFixed: self.workflow.center=(value,self.workflow.center[1])
    	self.display()
    	
    def changeYValue(self,value):
    	if not self.plotsettings.IsFixed: self.workflow.center=(self.workflow.center[0],value)
    	self.display()
    
    def changeEValue(self,value):
    	self.workflow.scale.ellipticity=value/100.
    	self.display()
    	
    def on_press(self,event):
        if event.inaxes and not self.plotsettings.IsFixed: 
            self.workflow.center=(event.xdata,event.ydata)
        self.canvas.mpl_connect('button_release_event', self.on_release)
    
    def on_release(self,event):
        if event.inaxes:
            x=event.xdata
            y=event.ydata
            if self.workflow.center==(event.xdata,event.ydata): self.workflow.r=10
            else:
	            #Limit the maximum r available to within the image.
    	        smalldim=int(min(self.workflow.raw.shape-np.array([self.workflow.center[1],self.workflow.center[0]])))
    	        smalldim=int(min(smalldim,min(self.workflow.center)))
    	        self.workflow.r=min(np.sqrt((x-self.workflow.center[0])**2+(y-self.workflow.center[1])**2),smalldim-1)/max(self.workflow.scale.ellipticity,1./self.workflow.scale.ellipticity)
            self.display() 
            
        
def cmap_xmap(function,cmap):
    """ Applies function, on the indices of colormap cmap. Beware, function
    should map the [0, 1] segment to itself, or you are in for surprises.
    """
    cdict = cmap._segmentdata.copy()
    function_to_map = lambda x : (function(x[0]), x[1], x[2])
    for key in ('red','green','blue'):         
        cdict[key] = map(function_to_map, cdict[key])
        cdict[key].sort()
        assert (cdict[key][0]<0 or cdict[key][-1]>1), "Resulting indices extend out of the [0, 1] segment."
    return colors.LinearSegmentedColormap('name_sqrt',cdict,2048)

def paint_circle(center,radius,ellipticity):
        theta=np.linspace(-np.pi,np.pi,1001)
        return center[0]+radius*np.cos(theta),center[1]+radius*np.sin(theta)*ellipticity

class CenterProcesser(QtCore.QThread):
    __errorHappened=False
    def __init__(self,gui,parent=None):
        QtCore.QThread.__init__(self,parent)
        self.workflow=gui.workflow
        self.gui=gui
        self.exiting=False
        
    def run(self):
        Cmax=0
        center,Cn=self.workflow.Newcenter(10)
        for i in np.arange(20):
            if Cn>Cmax:
        	self.workflow.center=center
        	self.gui.display()
        	self.emit(QtCore.SIGNAL("progress(int)"),i)
        	Cmax=Cn
        	center,Cn=self.workflow.Newcenter(10)
            else: 
                self.emit(QtCore.SIGNAL("progress(int)"),20)
                break
                
class InvertProcesser(QtCore.QThread):
    __errorHappened=False
    def __init__(self,gui,path,parent=None):
        QtCore.QThread.__init__(self,parent)
        self.workflow=gui.workflow
        self.gui=gui
        self.exiting=False
        self.path=path
        
    def run(self):
    	self.gui.statlabel.setText("Start the inversion procedure")
    	base=self.workflow.LoadBasis(self.path)
    	if len(base.shape)<2: 
    		QtGui.QMessageBox.warning(self,"No Basis","Basis file don't exist yet !!! \n Please build it first. :(")
    		return 0
    	self.gui.statlabel.setText("Basis Loaded")
    	self.emit(QtCore.SIGNAL("progress(int)"),5)
    	
    	if self.gui.DataChoiceBox.currentText()=="Original":
    		self.gui.statlabel.setText("Original") 
    		dat=np.copy(self.workflow.raw)
    	else:
    		self.gui.statlabel.setText("Current")
    		dat=dat=np.copy(self.workflow.datas)
    	
    	if not self.workflow.half:
    		dat=symmetrize(dat,self.workflow.center,self.workflow.r)
    		polar=self.workflow.to_polar(dat)
    	elif self.gui.WhichHalfBox.currentText()=="Left":
    		polar=self.workflow.to_polar(dat)
    	elif self.gui.WhichHalfBox.currentText()=="Right":
    		center=np.copy(self.workflow.center)
    		self.workflow.center=(dat.shape[0]-center[0],center[1])
    		polar=self.workflow.to_polar(np.fliplr(dat))
    		self.workflow.center=(center[0],center[1])
    	else: polar=self.workflow.to_polar(dat) 
    		
        self.gui.statlabel.setText("Polar Image")
    	self.emit(QtCore.SIGNAL("progress(int)"),10)
    	
        self.workflow.Invert(polar,base)
        self.gui.statlabel.setText("Fitted")
    	self.emit(QtCore.SIGNAL("progress(int)"),15)
    	del polar,base
    	
        self.workflow.image_for_display()
        self.gui.statlabel.setText("Image inverted")
    	self.emit(QtCore.SIGNAL("progress(int)"),20)

class PlotSettings():
    def __init__(self):
        self.palettename='Igor'
        self.IsSqrt=False
        self.IsR=False
        self.IsFixed=False
        cdict_coolheat ={'red'  :  ((0., 0., 0.), (0.25,0.,0.), (0.5,1.,1.), (0.75,1.0,1.0),  (1., 1., 1.)),
    					 'green':  ((0., 0., 0.), (0.25,0.,0.), (0.5,0.,0.), (0.75,1.0,1.0),  (1., 1., 1.)),
    					 'blue' :  ((0., 0., 0.), (0.25,1.,1.), (0.5,0.,0.), (0.75,0.0,0.0),  (1., 1., 1.))
    					}
    	cdict_coolheat_r ={'red'  :  ((0., 1., 1.), (0.25,1.,1.), (0.5,1,1.), (0.75,0,0),  (1., 0., 0.)),
    					   'green':  ((0., 1., 1.), (0.25,1.,1.), (0.5,0.,0.), (0.75,0,0),  (1., 0., 0.)),
    					   'blue' :  ((0., 1., 1.), (0.25,0.,0.), (0.5,0.,0.), (0.75,1,1),  (1., 0., 0.))
    					  }
    	coolheat = colors.LinearSegmentedColormap('mycoolheat', cdict_coolheat,2048)
    	coolheat_r = colors.LinearSegmentedColormap('mycoolheat_r', cdict_coolheat_r,2048)
    	
    	cdict_gray ={'red'  :  ((0., 0., 0.),(0.25, 0.5, 0.4), (1., 1., 1.)),
    				 'green':  ((0., 0., 0.),(0.25, 0.5, 0.4), (1., 1., 1.)),
    				 'blue' :  ((0., 0., 0.),(0.25, 0.5, 0.4), (1., 1., 1.))
    					}
    	cdict_gray_r ={'red'  :  ((0., 1., 1.),(0.25, 0.5, 0.6), (1., 0., 0.)),
    				   'green':  ((0., 1., 1.),(0.25, 0.5, 0.6), (1., 0., 0.)),
    				   'blue' :  ((0., 1., 1.),(0.25, 0.5, 0.6), (1., 0., 0.))
    				  }
    	gray = colors.LinearSegmentedColormap('mygray', cdict_gray,2048)
    	gray_r = colors.LinearSegmentedColormap('mygray_r', cdict_gray_r,2048)
    	
    	cdict_terrain ={'red': ((0.0, 0., 0.),(0.17, 1.0, 1.0),(0.3, 1.0, 1.0),(0.46, 0.0, 0.0),(0.68, 0.0, 0.0),(0.87, 1.0, 1.0),(0.89, 1.0, 1.0),(0.97, 0.0, 0.0),(1.0, 0.0, 0.0)),
 			   'green': ((0.0, 0., 0.),(0.31, 1.0, 1.0),(0.51, 1.0, 1.0),(0.74, 0.35, 0.35),(0.87, 1.0, 1.0),(0.9, 1.0, 1.0),(1., 0., 0.)),
 			   'blue': ((0.0, 0., 0.),(0.5, 0.0, 0.0),(0.87, 1.0, 1.0),(1.0, 1.0, 1.0))
    					}
    	cdict_terrain_r ={'red': ((0.0, 0., 0.),(0.17, 0.0, 0.0),(0.3, 1.0, 1.0),(0.46, 1.0, 1.0),(0.68, 0.0, 0.0),(0.87, 0.0, 0.0),(0.89, 1.0, 1.0),(0.97, 1.0, 1.0),(1.0, 0.0, 0.0)),
 						  'blue': ((0.0, 1., 1.),(0.5, 1.0, 1.0),(0.87, 0.0, 0.0),(1.0, 0.0, 0.0)),
 						  'green': ((0.0, 1., 1.),(0.25, 0.36, 0.36),(0.5, 1.0, 1.0),(0.75, 0.8, 0.8),(0.85, 0.6, 0.6),(1.0, 0.2, 0.2))
    					 }
    	terrain = colors.LinearSegmentedColormap('myterrain', cdict_terrain,2048)
    	terrain_r = colors.LinearSegmentedColormap('myterrain_r', cdict_terrain_r,2048)
    	
    	cdict_hot =cm.hot._segmentdata.copy()
    	cdict_hot_r =cm.hot_r._segmentdata.copy()
    	hot = colors.LinearSegmentedColormap('myhot', cdict_hot,2048)
    	hot_r = colors.LinearSegmentedColormap('myhot_r', cdict_hot_r,2048)
    	cdict_inferno =[[0.001462, 0.000466, 0.013866],
                 [0.002267, 0.001270, 0.018570],
                 [0.003299, 0.002249, 0.024239],
                 [0.004547, 0.003392, 0.030909],
                 [0.006006, 0.004692, 0.038558],
                 [0.007676, 0.006136, 0.046836],
                 [0.009561, 0.007713, 0.055143],
                 [0.011663, 0.009417, 0.063460],
                 [0.013995, 0.011225, 0.071862],
                 [0.016561, 0.013136, 0.080282],
                 [0.019373, 0.015133, 0.088767],
                 [0.022447, 0.017199, 0.097327],
                 [0.025793, 0.019331, 0.105930],
                 [0.029432, 0.021503, 0.114621],
                 [0.033385, 0.023702, 0.123397],
                 [0.037668, 0.025921, 0.132232],
                 [0.042253, 0.028139, 0.141141],
                 [0.046915, 0.030324, 0.150164],
                 [0.051644, 0.032474, 0.159254],
                 [0.056449, 0.034569, 0.168414],
                 [0.061340, 0.036590, 0.177642],
                 [0.066331, 0.038504, 0.186962],
                 [0.071429, 0.040294, 0.196354],
                 [0.076637, 0.041905, 0.205799],
                 [0.081962, 0.043328, 0.215289],
                 [0.087411, 0.044556, 0.224813],
                 [0.092990, 0.045583, 0.234358],
                 [0.098702, 0.046402, 0.243904],
                 [0.104551, 0.047008, 0.253430],
                 [0.110536, 0.047399, 0.262912],
                 [0.116656, 0.047574, 0.272321],
                 [0.122908, 0.047536, 0.281624],
                 [0.129285, 0.047293, 0.290788],
                 [0.135778, 0.046856, 0.299776],
                 [0.142378, 0.046242, 0.308553],
                 [0.149073, 0.045468, 0.317085],
                 [0.155850, 0.044559, 0.325338],
                 [0.162689, 0.043554, 0.333277],
                 [0.169575, 0.042489, 0.340874],
                 [0.176493, 0.041402, 0.348111],
                 [0.183429, 0.040329, 0.354971],
                 [0.190367, 0.039309, 0.361447],
                 [0.197297, 0.038400, 0.367535],
                 [0.204209, 0.037632, 0.373238],
                 [0.211095, 0.037030, 0.378563],
                 [0.217949, 0.036615, 0.383522],
                 [0.224763, 0.036405, 0.388129],
                 [0.231538, 0.036405, 0.392400],
                 [0.238273, 0.036621, 0.396353],
                 [0.244967, 0.037055, 0.400007],
                 [0.251620, 0.037705, 0.403378],
                 [0.258234, 0.038571, 0.406485],
                 [0.264810, 0.039647, 0.409345],
                 [0.271347, 0.040922, 0.411976],
                 [0.277850, 0.042353, 0.414392],
                 [0.284321, 0.043933, 0.416608],
                 [0.290763, 0.045644, 0.418637],
                 [0.297178, 0.047470, 0.420491],
                 [0.303568, 0.049396, 0.422182],
                 [0.309935, 0.051407, 0.423721],
                 [0.316282, 0.053490, 0.425116],
                 [0.322610, 0.055634, 0.426377],
                 [0.328921, 0.057827, 0.427511],
                 [0.335217, 0.060060, 0.428524],
                 [0.341500, 0.062325, 0.429425],
                 [0.347771, 0.064616, 0.430217],
                 [0.354032, 0.066925, 0.430906],
                 [0.360284, 0.069247, 0.431497],
                 [0.366529, 0.071579, 0.431994],
                 [0.372768, 0.073915, 0.432400],
                 [0.379001, 0.076253, 0.432719],
                 [0.385228, 0.078591, 0.432955],
                 [0.391453, 0.080927, 0.433109],
                 [0.397674, 0.083257, 0.433183],
                 [0.403894, 0.085580, 0.433179],
                 [0.410113, 0.087896, 0.433098],
                 [0.416331, 0.090203, 0.432943],
                 [0.422549, 0.092501, 0.432714],
                 [0.428768, 0.094790, 0.432412],
                 [0.434987, 0.097069, 0.432039],
                 [0.441207, 0.099338, 0.431594],
                 [0.447428, 0.101597, 0.431080],
                 [0.453651, 0.103848, 0.430498],
                 [0.459875, 0.106089, 0.429846],
                 [0.466100, 0.108322, 0.429125],
                 [0.472328, 0.110547, 0.428334],
                 [0.478558, 0.112764, 0.427475],
                 [0.484789, 0.114974, 0.426548],
                 [0.491022, 0.117179, 0.425552],
                 [0.497257, 0.119379, 0.424488],
                 [0.503493, 0.121575, 0.423356],
                 [0.509730, 0.123769, 0.422156],
                 [0.515967, 0.125960, 0.420887],
                 [0.522206, 0.128150, 0.419549],
                 [0.528444, 0.130341, 0.418142],
                 [0.534683, 0.132534, 0.416667],
                 [0.540920, 0.134729, 0.415123],
                 [0.547157, 0.136929, 0.413511],
                 [0.553392, 0.139134, 0.411829],
                 [0.559624, 0.141346, 0.410078],
                 [0.565854, 0.143567, 0.408258],
                 [0.572081, 0.145797, 0.406369],
                 [0.578304, 0.148039, 0.404411],
                 [0.584521, 0.150294, 0.402385],
                 [0.590734, 0.152563, 0.400290],
                 [0.596940, 0.154848, 0.398125],
                 [0.603139, 0.157151, 0.395891],
                 [0.609330, 0.159474, 0.393589],
                 [0.615513, 0.161817, 0.391219],
                 [0.621685, 0.164184, 0.388781],
                 [0.627847, 0.166575, 0.386276],
                 [0.633998, 0.168992, 0.383704],
                 [0.640135, 0.171438, 0.381065],
                 [0.646260, 0.173914, 0.378359],
                 [0.652369, 0.176421, 0.375586],
                 [0.658463, 0.178962, 0.372748],
                 [0.664540, 0.181539, 0.369846],
                 [0.670599, 0.184153, 0.366879],
                 [0.676638, 0.186807, 0.363849],
                 [0.682656, 0.189501, 0.360757],
                 [0.688653, 0.192239, 0.357603],
                 [0.694627, 0.195021, 0.354388],
                 [0.700576, 0.197851, 0.351113],
                 [0.706500, 0.200728, 0.347777],
                 [0.712396, 0.203656, 0.344383],
                 [0.718264, 0.206636, 0.340931],
                 [0.724103, 0.209670, 0.337424],
                 [0.729909, 0.212759, 0.333861],
                 [0.735683, 0.215906, 0.330245],
                 [0.741423, 0.219112, 0.326576],
                 [0.747127, 0.222378, 0.322856],
                 [0.752794, 0.225706, 0.319085],
                 [0.758422, 0.229097, 0.315266],
                 [0.764010, 0.232554, 0.311399],
                 [0.769556, 0.236077, 0.307485],
                 [0.775059, 0.239667, 0.303526],
                 [0.780517, 0.243327, 0.299523],
                 [0.785929, 0.247056, 0.295477],
                 [0.791293, 0.250856, 0.291390],
                 [0.796607, 0.254728, 0.287264],
                 [0.801871, 0.258674, 0.283099],
                 [0.807082, 0.262692, 0.278898],
                 [0.812239, 0.266786, 0.274661],
                 [0.817341, 0.270954, 0.270390],
                 [0.822386, 0.275197, 0.266085],
                 [0.827372, 0.279517, 0.261750],
                 [0.832299, 0.283913, 0.257383],
                 [0.837165, 0.288385, 0.252988],
                 [0.841969, 0.292933, 0.248564],
                 [0.846709, 0.297559, 0.244113],
                 [0.851384, 0.302260, 0.239636],
                 [0.855992, 0.307038, 0.235133],
                 [0.860533, 0.311892, 0.230606],
                 [0.865006, 0.316822, 0.226055],
                 [0.869409, 0.321827, 0.221482],
                 [0.873741, 0.326906, 0.216886],
                 [0.878001, 0.332060, 0.212268],
                 [0.882188, 0.337287, 0.207628],
                 [0.886302, 0.342586, 0.202968],
                 [0.890341, 0.347957, 0.198286],
                 [0.894305, 0.353399, 0.193584],
                 [0.898192, 0.358911, 0.188860],
                 [0.902003, 0.364492, 0.184116],
                 [0.905735, 0.370140, 0.179350],
                 [0.909390, 0.375856, 0.174563],
                 [0.912966, 0.381636, 0.169755],
                 [0.916462, 0.387481, 0.164924],
                 [0.919879, 0.393389, 0.160070],
                 [0.923215, 0.399359, 0.155193],
                 [0.926470, 0.405389, 0.150292],
                 [0.929644, 0.411479, 0.145367],
                 [0.932737, 0.417627, 0.140417],
                 [0.935747, 0.423831, 0.135440],
                 [0.938675, 0.430091, 0.130438],
                 [0.941521, 0.436405, 0.125409],
                 [0.944285, 0.442772, 0.120354],
                 [0.946965, 0.449191, 0.115272],
                 [0.949562, 0.455660, 0.110164],
                 [0.952075, 0.462178, 0.105031],
                 [0.954506, 0.468744, 0.099874],
                 [0.956852, 0.475356, 0.094695],
                 [0.959114, 0.482014, 0.089499],
                 [0.961293, 0.488716, 0.084289],
                 [0.963387, 0.495462, 0.079073],
                 [0.965397, 0.502249, 0.073859],
                 [0.967322, 0.509078, 0.068659],
                 [0.969163, 0.515946, 0.063488],
                 [0.970919, 0.522853, 0.058367],
                 [0.972590, 0.529798, 0.053324],
                 [0.974176, 0.536780, 0.048392],
                 [0.975677, 0.543798, 0.043618],
                 [0.977092, 0.550850, 0.039050],
                 [0.978422, 0.557937, 0.034931],
                 [0.979666, 0.565057, 0.031409],
                 [0.980824, 0.572209, 0.028508],
                 [0.981895, 0.579392, 0.026250],
                 [0.982881, 0.586606, 0.024661],
                 [0.983779, 0.593849, 0.023770],
                 [0.984591, 0.601122, 0.023606],
                 [0.985315, 0.608422, 0.024202],
                 [0.985952, 0.615750, 0.025592],
                 [0.986502, 0.623105, 0.027814],
                 [0.986964, 0.630485, 0.030908],
                 [0.987337, 0.637890, 0.034916],
                 [0.987622, 0.645320, 0.039886],
                 [0.987819, 0.652773, 0.045581],
                 [0.987926, 0.660250, 0.051750],
                 [0.987945, 0.667748, 0.058329],
                 [0.987874, 0.675267, 0.065257],
                 [0.987714, 0.682807, 0.072489],
                 [0.987464, 0.690366, 0.079990],
                 [0.987124, 0.697944, 0.087731],
                 [0.986694, 0.705540, 0.095694],
                 [0.986175, 0.713153, 0.103863],
                 [0.985566, 0.720782, 0.112229],
                 [0.984865, 0.728427, 0.120785],
                 [0.984075, 0.736087, 0.129527],
                 [0.983196, 0.743758, 0.138453],
                 [0.982228, 0.751442, 0.147565],
                 [0.981173, 0.759135, 0.156863],
                 [0.980032, 0.766837, 0.166353],
                 [0.978806, 0.774545, 0.176037],
                 [0.977497, 0.782258, 0.185923],
                 [0.976108, 0.789974, 0.196018],
                 [0.974638, 0.797692, 0.206332],
                 [0.973088, 0.805409, 0.216877],
                 [0.971468, 0.813122, 0.227658],
                 [0.969783, 0.820825, 0.238686],
                 [0.968041, 0.828515, 0.249972],
                 [0.966243, 0.836191, 0.261534],
                 [0.964394, 0.843848, 0.273391],
                 [0.962517, 0.851476, 0.285546],
                 [0.960626, 0.859069, 0.298010],
                 [0.958720, 0.866624, 0.310820],
                 [0.956834, 0.874129, 0.323974],
                 [0.954997, 0.881569, 0.337475],
                 [0.953215, 0.888942, 0.351369],
                 [0.951546, 0.896226, 0.365627],
                 [0.950018, 0.903409, 0.380271],
                 [0.948683, 0.910473, 0.395289],
                 [0.947594, 0.917399, 0.410665],
                 [0.946809, 0.924168, 0.426373],
                 [0.946392, 0.930761, 0.442367],
                 [0.946403, 0.937159, 0.458592],
                 [0.946903, 0.943348, 0.474970],
                 [0.947937, 0.949318, 0.491426],
                 [0.949545, 0.955063, 0.507860],
                 [0.951740, 0.960587, 0.524203],
                 [0.954529, 0.965896, 0.540361],
                 [0.957896, 0.971003, 0.556275],
                 [0.961812, 0.975924, 0.571925],
                 [0.966249, 0.980678, 0.587206],
                 [0.971162, 0.985282, 0.602154],
                 [0.976511, 0.989753, 0.616760],
                 [0.982257, 0.994109, 0.631017],
                 [0.988362, 0.998364, 0.644924]]
 	cdict_inferno_r = cdict_inferno[::-1]
    	inferno = colors.LinearSegmentedColormap.from_list('myinferno', cdict_inferno,2048)
    	inferno_r = colors.LinearSegmentedColormap.from_list('myinferno_r', cdict_inferno_r,2048)
    	
    	cdict_jet ={'blue': ((0.0, 0.5, 0.75),(0.11, 1, 1),(0.34, 1, 1),(0.65, 0, 0),(1, 0, 0.5)),
 					'green': ((0.0, 0, 0.5),(0.125, 0, 0.5),(0.375, 1, 1),(0.64, 1, 1),(0.91, 0, 0.5),(1, 0, 0.5)),
 					'red': ((0.0, 0, 0.5), (0.35, 0, 0.5), (0.66, 1, 1), (0.89, 1, 1), (1, 0.5, 0.75))
 					}
 	cdict_jet_r ={'blue': ((0.0, 0., 0.5),(0.35, 0, 0),(0.76, 1, 1),(0.89, 1, 1),(1, 0.5, 0.75)),
 					'green': ((0.0, 0, 0.5),(0.125, 0, 0.5),(0.36, 1, 1),(0.64, 1, 1),(0.91, 0, 0.5),(1, 0, 0.5)),
 					'red': ((0.0, 0.5, 0.75), (0.11, 1, 1), (0.44, 1, 1), (0.65, 0, 0.5), (1, 0., 0.5))
 					}
 	jet = colors.LinearSegmentedColormap('myjet', cdict_jet,2048)
    	jet_r = colors.LinearSegmentedColormap('myjet_r', cdict_jet_r,2048)
    	self.cmapdic={'Jet':jet,'Hot':hot,'Gray':gray,'Gnuplot':coolheat,'Igor':terrain,'Inferno':inferno}
    	self.cmapdic_sqrt={'Jet':cmap_xmap(lambda x: x**2,jet),'Hot':cmap_xmap(lambda x: x**2,hot),'Gray':cmap_xmap(lambda x: x**2,gray),'Gnuplot':cmap_xmap(lambda x: x**2,coolheat),'Igor':cmap_xmap(lambda x: x**2,terrain),'Inferno':cmap_xmap(lambda x: x**2,inferno)}
    	self.cmapdic_r={'Jet':jet_r,'Hot':hot_r,'Gray':gray_r,'Gnuplot':coolheat_r,'Igor':terrain_r,'Inferno':inferno_r}
    	self.cmapdic_r_sqrt={'Jet':cmap_xmap(lambda x: x**2,jet_r),'Hot':cmap_xmap(lambda x: x**2,hot_r),'Gray':cmap_xmap(lambda x: x**2,gray_r),'Gnuplot':cmap_xmap(lambda x: x**2,coolheat_r),'Igor':cmap_xmap(lambda x: x**2,terrain_r),'Inferno':cmap_xmap(lambda x: x**2,inferno_r)}