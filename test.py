import sys
from pBaseQT import pBaseForm
from PySide import QtCore
from PySide.QtGui import QMainWindow,QApplication

class MainWindow(QMainWindow):
  def __init__(self, parent=None):
    QMainWindow.__init__(self,parent)
    self.ui =  pBaseForm()
    self.ui.setupUi(self)

if __name__ == "__main__":
    app=QApplication(sys.argv)
    pBase = MainWindow()
    pBase.show()
    sys.exit(app.exec_())
    
    
    
    """
    compile library for 32bits and 34 bits arch by adding -arch i386 -arch x86_64
    """