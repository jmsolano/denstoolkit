#ifndef DTKMAINWINDOW_H
#define DTKMAINWINDOW_H

#include <QMainWindow>
#include <QAction>

namespace Ui {
   class DTKMainWindow;
}

class DTKMainWindow : public QMainWindow
{
   Q_OBJECT

public:
   explicit DTKMainWindow(QWidget *parent = 0);
   ~DTKMainWindow();

public slots:
   void updateStatusBar(void);
   void loadMolecule(void);
   void loadTestMolecule(void);
   void clearViewPort(void);
   void exportViewPortImage(void);

private slots:
   void on_resetPushButton_clicked();

   void on_viewAtLblsCheckBox_clicked();

private:
   Ui::DTKMainWindow *ui;
   // ++++++++++++++++++++ ACTIONS ++++++++++++++++++++
   QAction *loadMoleculeAction;
   QAction *loadTestMoleculeAction;
   QAction *clearViewPortAction;
   QAction *exportViewPortImageAction;
   // +++++++++++++++++++++ MENUS +++++++++++++++++++++
   QMenu *fileMenu;
   // +++++++++++++++++++++++++++++++++++++++++++++++++
   void createMenus(void);
   void createActions(void);
};

#endif // DTKMAINWINDOW_H
