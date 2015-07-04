#include "dtkglutils.h"
#include <cmath>


void dtkglutils::getRotationVectorAndAngle(const QVector3D &v1, const QVector3D &v2, QVector3D &vres, float &ares)
{
    const float oeoPI = 180.0f/acos(-1.0f);
    static const QVector3D zAxis(0.0f, 0.0f, 1.0f);
    QVector3D diff=v2-v1;
    vres=QVector3D::crossProduct(zAxis,diff);
    float radians=acos(QVector3D::dotProduct(zAxis,diff)/(diff.length()));
    ares = radians * oeoPI;
}
