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
#include <vtkXMLPUnstructuredGridWriter.h>

namespace olb {

template<typename T>
class vtkSurfaceWriter {
private:
  STLreader<T>& _surfaceI;
  CuboidGeometry3D<T>& _cuboidGeometry;
  LoadBalancer<T>& _loadBalancer;

  const std::string _fileName;
  std::vector<FunctorPtr<AnalyticalF3D<T,T>>> _f;

public:
  vtkSurfaceWriter(STLreader<T>& surfaceI,
                   CuboidGeometry3D<T>& cuboidGeometry,
                   LoadBalancer<T>& loadBalancer,
                   const std::string& fileName):
    _surfaceI(surfaceI),
    _cuboidGeometry(cuboidGeometry),
    _loadBalancer(loadBalancer),
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

  std::vector<const STLtriangle<T>*> localTriangles;
  std::size_t pointId = 0;
  for (const STLtriangle<T>& triangle : _surfaceI.getMesh().getTriangles()) {
    int globC{};
    if (_cuboidGeometry.getC(triangle.getCenter(), globC)) {
      if (_loadBalancer.isLocal(globC)) {
        localTriangles.emplace_back(&triangle);
        for (const auto& point : triangle.point) {
          points->InsertNextPoint(point.coords[0], point.coords[1], point.coords[2]);
        }

        vtkSmartPointer<vtkTriangle> cell = vtkSmartPointer<vtkTriangle>::New();
        cell->GetPointIds()->SetId(0, pointId++);
        cell->GetPointIds()->SetId(1, pointId++);
        cell->GetPointIds()->SetId(2, pointId++);
        cells->InsertNextCell(cell);
      }
    }
  }

  grid->SetPoints(points);
  grid->SetCells(VTK_TRIANGLE, cells);

  std::vector<vtkSmartPointer<vtkFloatArray>> data(_f.size());
  for (std::size_t iF=0; iF < _f.size(); ++iF) {
    data[iF] = vtkSmartPointer<vtkFloatArray>::New();
    data[iF]->SetName(_f[iF]->getName().c_str());
    data[iF]->SetNumberOfComponents(_f[iF]->getTargetDim());
    grid->GetPointData()->AddArray(data[iF]);
  }

  std::string filePath = singleton::directories().getVtkOutDir()
                       + createFileName(_fileName, iT)
                       + ".pvtu";

  if (singleton::mpi().isMainProcessor()) {
    vtkSmartPointer<vtkXMLPUnstructuredGridWriter> multiWriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
    multiWriter->SetInputData(grid);
    multiWriter->SetNumberOfPieces(singleton::mpi().getSize());
    multiWriter->SetFileName(filePath.c_str());
    multiWriter->SetUseSubdirectory(true);
    multiWriter->SetStartPiece(0);
    multiWriter->SetEndPiece(singleton::mpi().getSize()-1);
    multiWriter->SetWriteSummaryFile(1);
    multiWriter->Write();
  }

  for (const STLtriangle<T>* triangle : localTriangles) {
    for (const auto& point : triangle->point) {
      for (std::size_t iF=0; iF < _f.size(); ++iF) {
        auto physR = point.coords;
        std::vector<T> result(_f[iF]->getTargetDim(), 0);
        _f[iF](result.data(), physR.data());
        data[iF]->InsertNextTuple(result.data());
      }
    }
  }

  vtkSmartPointer<vtkXMLPUnstructuredGridWriter> multiWriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
  multiWriter->SetWriteSummaryFile(0);
  multiWriter->SetInputData(grid);
  multiWriter->SetNumberOfPieces(singleton::mpi().getSize());
  multiWriter->SetFileName(filePath.c_str());
  multiWriter->SetUseSubdirectory(true);
  multiWriter->SetStartPiece(singleton::mpi().getRank());
  multiWriter->SetEndPiece(singleton::mpi().getRank());
  multiWriter->Write();
}

}

#endif

#endif
