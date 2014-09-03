#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "nexus/NeXusFile.hpp"
#include <vector>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
using namespace std;

#include "log.h"
#include "config.h"

//const uint32_t MAX_TOF =4000;
//const uint32_t MAX_DET =1080;
const uint32_t MAX_TOF =3000;
const uint32_t MAX_DET =1080;
typedef  map<uint32_t, vector<uint32_t> > TofDetMap; 
typedef  vector<uint32_t> TofVector; 

uint32_t MapIdx(uint32_t tofidx, uint32_t detidx){
	return tofidx*MAX_DET+detidx;
}

uint32_t TofIdx(uint32_t mapidx){
	return mapidx/MAX_DET;
}

uint32_t DetIdx(uint32_t mapidx){
	return mapidx%MAX_DET;
}

void LoadSimulationFile(uint32_t* cmap, std::string samplefilename){

	std::cout << "LoadSimulationFile "<< std::endl;
	std::ifstream samplefile(samplefilename.c_str());
	string samplebuff;
	getline(samplefile, samplebuff);
	uint64_t tot=0;

	for (int tofidx=0;tofidx<MAX_TOF ;tofidx++){
		getline(samplefile, samplebuff);
		vector<string> substring;
		vector<double> counts;
		boost::split( substring, samplebuff, boost::is_any_of( ";" ), boost::token_compress_on );
		for(uint32_t detidx = 0; detidx < MAX_DET  ; detidx++){
			cmap[MapIdx(tofidx, detidx)] =atoi(substring[detidx+1].c_str());
			cmap[MapIdx(tofidx, detidx)] =(int)(atoi(substring[detidx+1].c_str()));
			tot+=(int)(atoi(substring[detidx+1].c_str()));
		}
	}
	std::cout << "total neutron hit count: " << tot << std::endl;
	samplefile.close();
}

void CreateTofDetMap(uint32_t* cmap, TofDetMap& tofdetmap, TofVector& tofvector){
	int tofidx = 0;


	for(uint32_t tofidx = 0; tofidx < MAX_TOF; tofidx++){
		tofdetmap.insert(make_pair(tofidx, std::vector<uint32_t>()));
		for(uint32_t detidx=0; detidx < MAX_DET; detidx++){
			if( cmap[MapIdx(tofidx, detidx)] > 0)tofdetmap[tofidx].push_back(detidx);
		}
	}


	for(uint32_t tofidx = 0; tofidx < MAX_TOF; tofidx++){
                if(tofdetmap[tofidx].size() > 0) tofvector.push_back(tofidx);
	}


}

void ReadOut(uint32_t tofidx, uint32_t detidx){
	uint32_t tof = 8000+tofidx*8;
        //std::cout << "TOF: " << tof << "; DET IDX: " << detidx << std::endl; 
}

void RandomHit(uint32_t* cmap, TofDetMap& tofdetmap, TofVector& tofvector){
	uint32_t tofidx = tofvector.at(uint32_t(rand()%(tofdetmap.size())));
	uint32_t detidx = tofdetmap[tofidx].at((rand()%(tofdetmap[tofidx].size())));
	//std::cout << "TOF IDX: " << tofidx << std::endl;
	//std::cout << "DET IDX: " << detidx << std::endl;
	//std::cout << "COUNT:   " << cmap[MapIdx(tofidx, detidx)] << std::endl;
	ReadOut(tofidx, detidx);
	cmap[MapIdx(tofidx, detidx)]--;
	//std::cout << "COUNT:   " << cmap[MapIdx(tofidx, detidx)] << std::endl;
	if(cmap[MapIdx(tofidx, detidx)]==0){
                //std::cout << "tofdetmap tofidx " << tofidx << " before remove " << tofdetmap[tofidx].size() << std::endl;
		tofdetmap[tofidx].erase(std::find(tofdetmap[tofidx].begin(), tofdetmap[tofidx].end(), detidx));
                //std::cout << "tofdetmap tofidx " << tofidx << " after remove " << tofdetmap[tofidx].size() << std::endl;
                //std::cout << "Remove TOFIDX: " << tofidx << "; DETIDX: " << detidx << std::endl;
                //if(tofidx == 1) std::cout << "tofidx 1 remove " << tofdetmap[tofidx].size() << std::endl;
		if(tofdetmap[tofidx].size()==0){
			tofdetmap.erase(tofidx);
                        tofvector.erase(std::find(tofvector.begin(), tofvector.end(), tofidx));
                        std::cout << "After  Remove TOF VECTOR: " << tofidx << " size " << tofdetmap.size() << std::endl;
		}
	}
}

int main(int argc, char *argv[]) {
	if (argc != 2){
		std::cout << "Usage: " << argv[0] << " option.txt" << std::endl;
		return 1;
	}
	std::string configfile(argv[1]);
	Config* fConfig = new Config(configfile);
	std::string samplefile  ; 
	samplefile  = fConfig->pString("samplefile") ;  
	std::cout << "read sample file: " << samplefile << std::endl;

	uint32_t *cmap = new uint32_t[MAX_TOF*MAX_DET];
	LoadSimulationFile(cmap, samplefile); 

	TofDetMap* tofdetmap = new TofDetMap;
	TofVector* tofvector = new TofVector;
	CreateTofDetMap(cmap, *tofdetmap, *tofvector);

	srand(time(NULL));
        uint64_t counts=0;
	while(1){
		RandomHit(cmap, *tofdetmap, *tofvector);
                counts++;
                if((counts%100000000)==0)std::cout<< "Read Out Counts: " << counts << std::endl;
                if (tofdetmap->size() == 0)break;
	}

	std::cout << "Total Events: " << counts << std::endl;
}
