#include <unordered_map>

template <typename T> class ReferenceData{

    private:
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

    public:
    ReferenceData(){};

    T getCharUMultiplier(T Ra, T charU){
        std::unordered_map<T, T> mult_map =
            {
                {1e3,  compVals[0][1]},
                {1e4,  compVals[1][1]},
                {1e5,  compVals[2][1]},
                {1e6,  compVals[3][1]},
                {1e7,  compVals[4][1]},
                {1e8,  compVals[5][1]},
                {1e9,  compVals[6][1]},
                {1e10, compVals[7][1]}
            };
        
        return mult_map.contains(Ra) ? mult_map[Ra] : charU;
    }

};