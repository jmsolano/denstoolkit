/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
          Copyright (c) 2013-2020, Juan Manuel Solano-Altamirano
                                   <jmsolanoalt@gmail.com>
  
   -------------------------------------------------------------------
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
   ---------------------------------------------------------------------
  
   If you want to redistribute modifications of the suite, please
   consider to include your modifications in our official release.
   We will be pleased to consider the inclusion of your code
   within the official distribution. Please keep in mind that
   scientific software is very special, and version control is 
   crucial for tracing bugs. If in despite of this you distribute
   your modified version, please do not call it DensToolKit.
  
   If you find DensToolKit useful, we humbly ask that you cite
   the paper(s) on the package --- you can find them on the top
   README file.
*/

#include "dtkglutils.h"
#include <QFileInfo>
#include <cmath>
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include "../common/iofuncts-cpx.h"


void dtkglutils::getRotationVectorAndAngle(const QVector3D &v1, const QVector3D &v2, QVector3D &vres, float &ares)
{
    const float oeoPI = 180.0f/acos(-1.0f);
    static const QVector3D zAxis(0.0f, 0.0f, 1.0f);
    QVector3D diff=v2-v1;
    vres=QVector3D::crossProduct(zAxis,diff);
    float radians=acos(QVector3D::dotProduct(zAxis,diff)/(diff.length()));
    ares = radians * oeoPI;
}


QString dtkglutils::getWFNFileNameFromCPX(QString cpxname)
{
   QFileInfo fi(cpxname);
   string str=cpxname.toStdString();
   ifstream ifil(str.c_str());
   str=cpxGetWFXFileName(ifil);
   ifil.close();
   QString res=fi.absolutePath()+QString("/")+QString::fromStdString(str);
   return res;
}
