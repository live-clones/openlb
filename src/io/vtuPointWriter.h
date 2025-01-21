#ifndef VTU_POINT_WRITER_H
#define VTU_POINT_WRITER_H


namespace olb {

struct position{
  double x,y;
  double z = 0.0;
};
template<typename T>
class vtuPointWriter{
protected:
  int _dim = 3;
  std::string _name;
  bool _binary;
  bool _haveMaster;
  mutable OstreamManager clout;




  //  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  //  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);
  ///  performes <VTKFile ...>, <ImageData ...> and <PieceExtent ...>
  void preambleVTU(const std::string& fullName, int num);
  ///  performes </ImageData> and </VTKFile>
  void closeVTU(const std::string& fullNamePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVD(int iT, int i, const std::string& fullNamePVD,
               const std::string& namePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVDmaster(int iT, int i, const std::string& fullNamePVDMaster,
                     const std::string& namePiece);
  ///  writes functors stored at pointerVec
  void dataArray(const std::string& fullName, std::vector<position> pos, std::vector<double> val);

  void writePosition( std::ofstream& fout, std::vector<position> pos );
  void writeConst( std::ofstream& fout , std::vector<olb::position> pos, std::vector<double> val );
public:
  vtuPointWriter(  std::string const name,
                            bool binary = true);
  vtuPointWriter( const vtuPointWriter<T>& rhs);
  vtuPointWriter( const vtuPointWriter<T>&& rhs);

  void createMasterFile();
  void write(std::size_t iT, std::vector<position> pos, std::vector<double> val);

};

template<typename T>
vtuPointWriter<T>::vtuPointWriter( std::string const name,
                                                 bool binary )
  : _name( name ), _binary( binary ),
    _haveMaster(false), clout(std::cout, "vtuPointWriter")
{}

template<typename T>
vtuPointWriter<T>::vtuPointWriter(const vtuPointWriter<T>& rhs)
  : _name(rhs._name), _binary(rhs._binary),
    _haveMaster(rhs._haveMaster), clout(std::cout, "vtuPointWriter"), _dim( rhs._dim)
{}

template<typename T>
vtuPointWriter<T>::vtuPointWriter(const vtuPointWriter<T>&& rhs)
  : _name(rhs._name), _binary(rhs._binary),
    _haveMaster(rhs._haveMaster), clout(std::cout, "vtuPointWriter"), _dim( rhs._dim )
{}



template<typename T>
void vtuPointWriter<T>::write( std::size_t iT, std::vector<position> pos, std::vector<double> val )
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  if (rank == 0) { // master only
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName(_name) + "_master.pvd";
    std::string fullNamePVD = singleton::directories().getVtkOutDir() + "data/"
                              + createFileName(_name, iT) + ".pvd";
    preamblePVD( fullNamePVD );         // timestep
    std::string namePiece =  "data/" + createFileName(_name, iT, 0) + ".vtu";
    // puts name of .vti piece to a .pvd file [fullNamePVD]
    dataPVD(iT, 1, fullNamePVD, namePiece);
    // adds a namePiece to master.pvd file.
    // To do so we overwrite closePVD() and add new entry.
    dataPVDmaster(iT, 1, fullNamePVDmaster, namePiece);
    closePVD(fullNamePVD);            // timestep
  } // master only
  if ( rank == 0){ //master only or different systems each
    std::string fullNameVTU = singleton::directories().getVtkOutDir()
                            + "data/" + createFileName(_name, iT, rank) + ".vtu";
    preambleVTU(fullNameVTU, pos.size());

    this->dataArray(fullNameVTU, pos, val);
    closeVTU(fullNameVTU);
}
}


template<typename T>
void vtuPointWriter<T>::createMasterFile()
{
  std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                  + createFileName(_name) + "_master.pvd";
  preamblePVD(fullNamePVDmaster);
  closePVD(fullNamePVDmaster);
  _haveMaster = true;

}

template<typename T>
void vtuPointWriter<T>::preamblePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n" << "<Collection>\n";
  fout.close();
}


template<typename T>
void vtuPointWriter<T>::closePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}


template<typename T>
void vtuPointWriter<T>::preambleVTU(
  const std::string& fullName, int num)
{
  std::ofstream fout(fullName.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>" << std::endl << std::flush;
  fout
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
      << std::endl;
  fout << "<UnstructuredGrid>" << std::endl;
  fout << "<Piece NumberOfPoints=\"" << num
       << "\" NumberOfCells=\"" << num << "\">"
       << std::endl;
  fout << "<PointData Vectors=\"Particles\">" << std::endl;
  fout.close();
}

template<typename T>
void vtuPointWriter<T>::closeVTU(
  const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</UnstructuredGrid>\n";
  fout << "</VTKFile>\n";
  fout.close();
}


template<typename T>
void vtuPointWriter<T>::dataPVD(int iT, int i,
    const std::string& fullNamePVD, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" " << i
       << "\" " << "file=\"" << namePiece << "\"/>\n";
  fout.close();
}

template<typename T>
void vtuPointWriter<T>::dataPVDmaster(int iT, int i,
    const std::string& fullNamePVDMaster, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVDMaster.c_str(),
                     std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25, std::ios::end); // jump -25 form the end of file to overwrite closePVD
    fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" "
         << i << "\" " << "file=\"" << namePiece << "\"/>\n";
    fout.close();
    closePVD(fullNamePVDMaster);
  } else {
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
}


template<typename T>
void vtuPointWriter<T>::dataArray(
    const std::string& fullName, std::vector<position> pos, std::vector<double> val)
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  writeConst(fout, pos, val);
  fout << "</PointData>" << std::endl;
  fout << "<CellData /> " << std::endl;
  fout << "<Cells>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
       << std::endl;
  for ( int i=0; i < 1; ++i){
    for ( unsigned long j=0; j<pos.size(); ++j){
      fout << j << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
       << std::endl;
  for ( int i=0; i < 1; ++i){
    for ( unsigned long j=1; j <= pos.size(); ++j){
      fout << j << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
       << std::endl;
    for ( unsigned long j=0; j < pos.size(); ++j){
      fout << 1 << " ";
    }
  fout << "</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
  fout << "<Points>" << std::endl;
  writePosition( fout, pos );
  fout << "</Points>" << std::endl;
  fout << "</Piece>" << std::endl;

  fout.close();
}


template<typename T>
void vtuPointWriter<T>::writePosition( std::ofstream& fout , std::vector<position> pos ){
  fout << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\""<< _dim
       << "\">" << std::endl;

    int num = pos.size();
    for ( int j=0; j < num; ++j){
      fout << pos[j].x << " "<<pos[j].y << " "<<pos[j].z << " ";
    }
  fout << "</DataArray>" << std::endl;
}


template<typename T>
void vtuPointWriter<T>::writeConst( std::ofstream& fout , std::vector<olb::position> pos, std::vector<double> val ){
  fout << "<DataArray type=\"Float32\" Name=\"Constant\" NumberOfComponents=\""<< 1
       << "\">" << std::endl;

    int num = pos.size();
    for ( int j=0; j < num; ++j){
      fout << val[j] << " ";
    }
  fout << "</DataArray>" << std::endl;
}


}

#endif
