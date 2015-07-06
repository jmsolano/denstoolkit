#ifndef DTKGLUTILS_H
#define DTKGLUTILS_H
#include <QVector3D>
#include <QString>
namespace dtkglutils{
void getRotationVectorAndAngle(const QVector3D &v1,const QVector3D &v2,\
                                  QVector3D &vres,float &ares);
QString getWFNFileNameFromCPX(QString cpxname);
}
#endif // DTKGLUTILS_H

