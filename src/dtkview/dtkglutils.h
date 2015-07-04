#ifndef DTKGLUTILS_H
#define DTKGLUTILS_H
#include <QVector3D>
namespace dtkglutils{
void getRotationVectorAndAngle(const QVector3D &v1,const QVector3D &v2,\
                                  QVector3D &vres,float &ares);
}
#endif // DTKGLUTILS_H

