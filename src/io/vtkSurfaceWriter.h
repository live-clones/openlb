/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef IO_VTK_SURFACE_WRITER_H
#define IO_VTK_SURFACE_WRITER_H

#ifdef FEATURE_VTK

#include "stlReader.h"
#include "fileName.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkXMLUnstructuredGridWriter.h>

namespace olb {

template<typename T>
class vtkSurfaceWriter {
private:
  STLreader<T>& _surfaceI;
  const std::string _fileName;
  std::vector<FunctorPtr<AnalyticalF3D<T,T>>> _f;

public:
  vtkSurfaceWriter(STLreader<T>& surfaceI,
                   const std::string& fileName):
    _surfaceI(surfaceI),
    _fileName(fileName)
  { }

  void addFunctor(FunctorPtr<AnalyticalF3D<T,T>>&& f) {
    _f.emplace_back(std::move(f));
  }

  void write(int iT);

};

template<typename T>
void vtkSurfaceWriter<T>::write(int iT)
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

  std::size_t pointId = 0;
  for (const STLtriangle<T>& triangle : _surfaceI.getMesh().getTriangles()) {
    for (const auto& point : triangle.point) {
      points->InsertNextPoint(point.coords[0], point.coords[1], point.coords[2]);
    }

    vtkSmartPointer<vtkTriangle> cell = vtkSmartPointer<vtkTriangle>::New();
    cell->GetPointIds()->SetId(0, pointId++);
    cell->GetPointIds()->SetId(1, pointId++);
    cell->GetPointIds()->SetId(2, pointId++);
    cells->InsertNextCell(cell);
  }
  grid->SetPoints(points);
  grid->SetCells(VTK_TRIANGLE, cells);

  for (auto& f : _f) {
    vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetName(f->getName().c_str());
    data->SetNumberOfComponents(f->getTargetDim());

    for (const STLtriangle<T>& triangle : _surfaceI.getMesh().getTriangles()) {
      for (const auto& point : triangle.point) {
        auto physR = point.coords;
        std::vector<T> result(f->getTargetDim(), 0);
        f(result.data(), physR.data());
        data->InsertNextTuple(result.data());
      }
    }
    grid->GetPointData()->AddArray(data);
  }

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  std::string filePath = singleton::directories().getVtkOutDir()
                       + "data/"
                       + createFileName(_fileName, iT)
                       + ".vtu";
  writer->SetFileName(filePath.c_str());
  writer->SetInputData(grid);
  writer->Write();
}

}

#endif

#endif
