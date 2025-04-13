#include "olb3D.h"
namespace olb {

typedef enum {horizontal, vertical, diagonal2d, diagonal3d} SamplingDirection;

template<unsigned int ndim, typename T>
void linePlot( AnalyticalF3D<T,T>& data, size_t ndatapoints, T dist,
  std::string title, std::string ylabel,
  SamplingDirection direction, bool halfDomain = true,
  bool setRange = false, T ymin=0, T ymax=0 )
{
  Gnuplot<T> gplot( title );
  gplot.setLabel( "distance [m]", ylabel );
  int nmin = 0;
  if ( !halfDomain ) nmin = -int( ndatapoints/2 );
  for ( int n = nmin; n <= int(ndatapoints/2); n++ ) {
  T input[ndim] = {0,0,0};
  T distance = 0;
  switch ( direction ) {
    case horizontal:  input[0] = n*dist;                                        distance = n*dist;              break;
    case vertical:    input[1] = n*dist;                                        distance = n*dist;              break;
    case diagonal2d:  input[0] = n*dist; input[1] = n*dist;                     distance = n*dist*std::sqrt(2); break;
    case diagonal3d:  input[0] = n*dist; input[1] = n*dist; input[2] = n*dist;  distance = n*dist*std::sqrt(3); break;
  }
  T output[1];
  data( output, input );
  gplot.setData( distance, output[0] );
  }
  if ( setRange ) gplot.setYrange( ymin, ymax );
    gplot.writePNG( -1, -1, title );
  };
}