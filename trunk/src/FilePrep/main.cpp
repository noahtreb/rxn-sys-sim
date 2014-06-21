#include "main.h"
#include <netcdfcpp.h>

using namespace std;

int main(int argc, const char* argv[]) {
    string fileName;
    int maxDataSavePts, boundHandlingMethod;
    double stoppingTol;
    
    if (argc == 8) {
        for (int i = 0; i < argc; i++) {            
            if (argv[i][0] == '-') {
                switch(argv[i][1]) {
                    case 'n':
                        maxDataSavePts = atoi(argv[i+1]);
                        i++;
                        break;
                    case 'e':
                        stoppingTol = atof(argv[i+1]);
                        i++;
                        break;
                    case 'b':
                        boundHandlingMethod = atof(argv[i+1]);
                        i++;
                        break;
                }
            } else {                
                fileName = argv[i];
            }
        }
    } else {
        fprintf(stderr, "\nUsage: ./FilePrep <fileName> -n <maxDataSavePts> -e <stoppingTol> -b <boundHandlingMethod>\n");
        fprintf(stderr, "                  <fileName> Name of the NetCDF file to be prepared for simulation.\n");
        fprintf(stderr, "    -n      <maxDataSavePts> Maximum number of initial states to be held in the file.\n");
        fprintf(stderr, "    -e         <stoppingTol> Maximum allowable L2 distance between successive probability distributions.\n\n");
        fprintf(stderr, "    -b <boundHandlingMethod> 1 for deletion. 2 for resampling.");
        abort();
    }
        
    NcFile file(fileName.c_str(), NcFile::Write);
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", fileName.c_str());
        abort();
    }
        
    NcDim* dim = file.add_dim("maxDataSavePts", maxDataSavePts);
    
    NcVar* stoppingTolVar = file.add_var("stoppingTol", ncDouble);
    NcVar* numDataSavePtsVar = file.add_var("numDataSavePts", ncInt, file.get_dim("numMdls")); 
    NcVar* boundHandlingMethodVar = file.add_var("boundHandlingMethod", ncInt, file.get_dim("numMdls"));
    file.add_var("numBoundedSpeciesStates", ncInt, file.get_dim("numMdls"));
    file.add_var("speciesStateBounded", ncInt, file.get_dim("numMdls"), file.get_dim("maxSpecies"));
    file.add_var("speciesStateLowerBounds", ncDouble, file.get_dim("numMdls"), file.get_dim("maxSpecies"));
    file.add_var("speciesStateUpperBounds", ncDouble, file.get_dim("numMdls"), file.get_dim("maxSpecies"));
    file.add_var("dataSavePts", ncInt, file.get_dim("numMdls"), dim);
    file.add_var("speciesStateChanges", ncInt, file.get_dim("numMdls"), file.get_dim("maxSpecies"));
    
    file.get_var("state")->rename("initFwdData");    
    file.add_var("initRevData", ncDouble, file.get_dim("numMdls"), file.get_dim("maxTrials"), file.get_dim("maxTimePts"), file.get_dim("maxSaveSpecies"));
    file.add_var("fwdData", ncDouble, file.get_dim("numMdls"), dim, file.get_dim("maxTrials"), file.get_dim("maxTimePts"), file.get_dim("maxSaveSpecies"));
    file.add_var("revData", ncDouble, file.get_dim("numMdls"), dim, file.get_dim("maxTrials"), file.get_dim("maxTimePts"), file.get_dim("maxSaveSpecies"));
    
    file.add_var("initFwdAbsCurr", ncInt, file.get_dim("numMdls"), file.get_dim("maxTimePts"), file.get_dim("maxSaveSpecies"));
    file.add_var("initRevAbsCurr", ncInt, file.get_dim("numMdls"), file.get_dim("maxTimePts"), file.get_dim("maxSaveSpecies"));
    file.add_var("fwdAbsCurr", ncInt, file.get_dim("numMdls"), dim, file.get_dim("maxTimePts"), file.get_dim("maxSaveSpecies"));
    file.add_var("revAbsCurr", ncInt, file.get_dim("numMdls"), dim, file.get_dim("maxTimePts"), file.get_dim("maxSaveSpecies"));
    
    stoppingTolVar->put(&stoppingTol, 1);
    
    int numMdls = (int) file.get_dim("numMdls")->size();    
    for (int i = 0; i < numMdls; i++) {
        numDataSavePtsVar->set_cur(i);
        numDataSavePtsVar->put(&maxDataSavePts, 1);
    }
    
    boundHandlingMethodVar->put(&boundHandlingMethod, 1);
}