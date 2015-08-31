#include "dtkmainwindow.h"
#include <QApplication>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

int main(int argc, char *argv[])
{
   glutInit(&argc,argv);
   QApplication a(argc, argv);
   DTKMainWindow w;
   w.setWindowIcon(QIcon(":/images/dtk32x32.png"));
   w.show();

   return a.exec();
}
