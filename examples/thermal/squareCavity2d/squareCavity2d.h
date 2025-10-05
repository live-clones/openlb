#include <olb.h>
using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

bool comparable = true;

int N; // given resolution and current resolution
T Pr, PrTurb, Ra, tau, thermalExpansion, thermalDiffusivity, density, heatCapacity, kinematicViscosity, gravitationalConstant, charL, lx, charU, charUConfig;
T physLambda;
T maxPhysT, epsilon;
T smagoConst;

#ifdef TURBULENT
const int statisticsIntervall = 10; // take the turbulent statistics every 10 time steps after convergence
const int statisticsEnsembles = 20; // take 20 ensembles for the turbulent statistics
#endif

T Tcold, Thot, Tmean;

/// Values from the literature studies from De Vahl Davis (1983):
T compVals[8][6] = {
    { 3.649, 3.696, 1.013, 0.813, 0.178, 1.117          },
    { 16.178, 19.617, 1.212, 0.823, 0.119, 2.238        },
    { 34.730, 68.590, 1.975, 0.855, 0.066, 4.509        },
    { 64.530, 219.36, 3.400, 0.850, 0.036, 8.817        },
    { 164.24, 701.92, 4.831, 0.851, 0.020, 16.790       },
    { 389.88, 2241.37, 5.749, 0.937, 0.011, 30.506      },
    { 503.24, 6820.07, 13.552, 0.966, 0.0064, 57.350    },
    { 2323.00, 21463.00, 9.239, 0.940, 0.491, 103.663   }
};

T compParameters[6];

T charUMultiplier()
{
  if (Ra==1e3) {
    return compVals[0][1];
  }
  else if (Ra==1e4) {
    return compVals[1][1];
  }
  else if (Ra==1e5) {
    return compVals[2][1];
  }
  else if (Ra==1e6) {
    return compVals[3][1];
  }
  else if (Ra==1e7) {
    return compVals[4][1];
  }
  else if (Ra==1e8) {
    return compVals[5][1];
  }
  else if (Ra==1e9) {
    return compVals[6][1];
  }
  else if (Ra==1e10) {
    return compVals[7][1];
  }
  else {
    return charUConfig;
  }
}

void chooseComparisonParameters(T Ra)
{
    int row;
    if (Ra == 1e3) {
      row = 0;
    }
    else if (Ra == 1e4) {
      row = 1;
    }
    else if (Ra == 1e5) {
      row = 2;
    }
    else if (Ra == 1e6) {
      row = 3;
    }
    else if (Ra == 1e7) {
      row = 4;
    }
    else if (Ra == 1e8) {
      row = 5;
    }
    else if (Ra == 1e9) {
      row = 6;
    }
    else if (Ra == 1e10) {
      row = 7;
    } else {
        comparable = false;
    }
    if (comparable) {
      for (int iD = 0; iD < 6; iD++)
      {
        compParameters[iD] = compVals[row][iD];
      }
    }
}

void compareToStudy(Vector<T, 18>& output, Vector<T, 5>& params)
{
  OstreamManager clout(std::cout,"compareToStudy");
  T maxVelX       = params[0];
  T maxVelY       = params[1];
  T maxXVelYCoord = params[2];
  T maxYVelXCoord = params[3];
  T nusselt       = params[4];
  output[0] = maxVelX / physLambda * charL;
  output[1] = compParameters[0];
  output[2] = util::fabs((compParameters[0] - maxVelX / physLambda * charL) / compParameters[0]);

  output[3] = maxVelY / physLambda * charL;
  output[4] = compParameters[1];
  output[5] = util::fabs((compParameters[1] - maxVelY / physLambda * charL) / compParameters[1]);

  output[6] = maxVelY / maxVelX;
  output[7] = compParameters[2];
  output[8] = util::fabs((compParameters[2] - maxVelY / maxVelX)  / compParameters[2]);

  output[9] = maxXVelYCoord/lx;
  output[10] = compParameters[3];
  output[11] = util::fabs((compParameters[3] - maxXVelYCoord / lx) / compParameters[3]);

  output[12] = maxYVelXCoord/lx;
  output[13] = compParameters[4];
  output[14] = util::fabs((compParameters[4] - maxYVelXCoord / lx) / compParameters[4]);

  output[15] = nusselt;
  output[16] = compParameters[5];
  output[17] = util::fabs((compParameters[5] - nusselt) / compParameters[5]);
  if (comparable && singleton::mpi().isMainProcessor())
  {
    clout << "Comparison against De Vahl Davis (1983):" << std::endl;
    clout << "xVelocity in yDir=" <<  maxVelX / physLambda * charL << "; error(rel)=" << (T) util::fabs((compParameters[0] - maxVelX / physLambda * charL) / compParameters[0]) << std::endl;
    clout << "yVelocity in xDir=" <<  maxVelY / physLambda * charL << "; error(rel)=" << (T) util::fabs((compParameters[1] - maxVelY / physLambda * charL) / compParameters[1]) << std::endl;
    clout << "yMaxVel / xMaxVel="  <<  maxVelY / maxVelX << "; error(rel)=" << (T) util::fabs((compParameters[2] - maxVelY / maxVelX)  / compParameters[2]) << std::endl;
    clout << "yCoord of xMaxVel=" <<  maxXVelYCoord/lx << "; error(rel)=" << (T) util::fabs((compParameters[3] - maxXVelYCoord / lx) / compParameters[3]) << std::endl;
    clout << "xCoord of yMaxVel=" <<   maxYVelXCoord/lx << "; error(rel)=" << (T) util::fabs((compParameters[4] - maxYVelXCoord / lx) / compParameters[4]) << std::endl;
    clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) util::fabs((compParameters[5] - nusselt) / compParameters[5]) << std::endl;
    std::fstream fs;
    fs.open("output.txt",
      std::fstream::in | std::fstream::out | std::fstream::app);
    fs << "Comparison against De Vahl Davis (1983):" << std::endl;
    fs << "xVelocity in yDir=" <<  maxVelX / physLambda * charL << "; error(rel)=" << (T) util::fabs((compParameters[0] - maxVelX / physLambda * charL) / compParameters[0]) << std::endl;
    fs << "yVelocity in xDir=" <<  maxVelY / physLambda * charL << "; error(rel)=" << (T) util::fabs((compParameters[1] - maxVelY / physLambda * charL) / compParameters[1]) << std::endl;
    fs << "yMaxVel / xMaxVel="  <<  maxVelY / maxVelX << "; error(rel)=" << (T) util::fabs((compParameters[2] - maxVelY / maxVelX)  / compParameters[2]) << std::endl;
    fs << "yCoord of xMaxVel=" <<  maxXVelYCoord/lx << "; error(rel)=" << (T) util::fabs((compParameters[3] - maxXVelYCoord / lx) / compParameters[3]) << std::endl;
    fs << "xCoord of yMaxVel=" <<   maxYVelXCoord/lx << "; error(rel)=" << (T) util::fabs((compParameters[4] - maxYVelXCoord / lx) / compParameters[4]) << std::endl;
    fs << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) util::fabs((compParameters[5] - nusselt) / compParameters[5]) << std::endl;
    fs.close();
  } else {
    clout << "Comparison against De Vahl Davis (1983) not possible. Received Ra = " << Ra << "." << std::endl;
  }
}
