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
   w.show();

   return a.exec();
}
