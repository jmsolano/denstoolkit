/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.2.0
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
 */

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
   void setTransparentAtomsAndLinks(void);
   void setViewAtoms(void);
   void setViewAtomLabels(void);
   void setViewRegularBonds(void);
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
   void on_viewRegularBondsCheckBox_clicked();
   void on_viewAtomsCheckBox_clicked();
   void on_setTransparentCheckBox_clicked();

private:
   Ui::DTKMainWindow *ui;
   // ++++++++++++++++++++ ACTIONS ++++++++++++++++++++
   // When adding new actions, do not forget to delete
   // them within ~DTKMainWindow().
   QAction *loadMoleculeAction;
   QAction *loadTestMoleculeAction;
   QAction *clearViewPortAction;
   QAction *exportViewPortImageAction;

   QAction *viewAtomsAction;
   QAction *viewAtomLabelsAction;
   QAction *viewRegularBondsAction;
   QAction *setTransparentAtomsAndLinksAction;

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
