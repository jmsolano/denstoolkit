#ifndef DTKMAINWINDOW_H
#define DTKMAINWINDOW_H

#include <QMainWindow>
#include <QAction>
#include <QToolBar>

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
   void setViewAtomLabels(void);

private slots:
   void on_resetPushButton_clicked();

   void on_viewAtLblsCheckBox_clicked();

private:
   Ui::DTKMainWindow *ui;
   // ++++++++++++++++++++ ACTIONS ++++++++++++++++++++
   // When adding new actions, do not forget to delete
   // them within ~DTKMainWindow().
   QAction *loadMoleculeAction;
   QAction *loadTestMoleculeAction;
   QAction *clearViewPortAction;
   QAction *exportViewPortImageAction;
   QAction *viewAtomLabelsAction;
   // +++++++++++++++++++++ MENUS +++++++++++++++++++++
   // When adding new menus, do not forget to delete
   // them within ~DTKMainWindow().
   QMenu *fileMenu;
   QMenu *viewMenu;
   // +++++++++++++++++++++++++++++++++++++++++++++++++
   // +++++++++++++++++++ TOOLBARS ++++++++++++++++++++
   // When adding new toolbars, do not forget to delete
   // them within ~DTKMainWindow().
   QToolBar *fileToolBar;
   // +++++++++++++++++++++++++++++++++++++++++++++++++
   void createMenus(void);
   void createActions(void);
   void setupMainToolbar(void);
};

#endif // DTKMAINWINDOW_H
