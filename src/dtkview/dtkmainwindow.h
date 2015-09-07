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
   void setViewBondGradientPaths(void);
   void setViewRingGradientPaths(void);
   void setViewCageGradientPaths(void);
   void showAboutDTK(void);

private slots:
   void on_resetPushButton_clicked();
   void on_viewAtLblsCheckBox_clicked();
   void on_viewCGPsCheckBox_clicked();
   void on_viewRGPsCheckBox_clicked();
   void on_viewBGPsCheckBox_clicked();

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
   QAction *viewBondGradientPathsAction;
   QAction *viewRingGradientPathsAction;
   QAction *viewCageGradientPathsAction;
   QAction *showAboutDTKAction;
   // +++++++++++++++++++++ MENUS +++++++++++++++++++++
   // When adding new menus, do not forget to delete
   // them within ~DTKMainWindow().
   QMenu *fileMenu;
   QMenu *viewMenu;
   QMenu *helpMenu;
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
