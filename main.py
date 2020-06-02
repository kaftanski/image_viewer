import sys

from PyQt5 import QtWidgets

from project_ui_prototype import ProjectMainWindow


# def except_hook(cls, exception, traceback):
#     sys.__excepthook__(cls, exception, traceback)


if __name__ == "__main__":
    #sys.excepthook = except_hook

    app = QtWidgets.QApplication(sys.argv)
    window = ProjectMainWindow()
    sys.exit(app.exec_())
