import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

import joblib
import numpy as np
from plug import SLEder

class DlgMain(QDialog):
    def __init__(self):
        super().__init__()
        self.resize(300, 300)
        self.setWindowTitle("SLEPLUG v.1.0")
        #
        # pal = QPalette()
        # role = QPalette.Background
        # pal.setColor(role, QColor(200, 200, 200))
        # self.setPalette(pal)
        #
        self.main_text = QtWidgets.QLabel(self)
        self.main_text.setText("WELCOME TO SLEPLUG")
        self.main_text.move(70, 100)
        self.main_text.adjustSize()
        #
        self.btn = QPushButton('Open file', self)
        self.btn.move(105, 130)
        self.btn.clicked.connect(self.evt_btn_clicked)

    def evt_btn_clicked(self):
        path_to_file = QFileDialog.getOpenFileName(self,
                                          "Find the patient's file",
                                          '~/SLEder', 'TXT File (*.txt)')[0]
        name, answer, confidence = SLEder("best_model_all_fech_lr_29_08_22.pkl",
                                          "scaling_mean.txt",
                                          "scaling_var.txt",
                                          path_to_file)
        QMessageBox.information(self, "Conclusion",
                                'Name: '+name+"\n""Status: "+answer+' | confidence: '+str(confidence)+'%')

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dlgMain = DlgMain()
    dlgMain.show()
    sys.exit(app.exec_())
