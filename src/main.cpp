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
const uint32_t MAX_TOF =4000;
const uint32_t MAX_DET =61440;

typedef  map<uint32_t, vector<uint32_t> > TofDetMap; 
typedef  map<uint32_t, vector<uint32_t> > DetMapping; 
typedef  vector<uint32_t> TofVector; 

typedef vector<uint32_t> EvtLocation;
typedef multimap<uint32_t, EvtLocation> EvtInfo;

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
  
	std::cout << "create tofvector " << tofvector.size() << std::endl;


}

void ReadOut(DetMapping& detmapping, uint32_t tofidx, uint32_t detidx, EvtInfo& evtinfo){
	uint32_t tof = 8000+tofidx*8;
        EvtLocation evtlocation = detmapping[detidx];
	//std::cout << "TOF: " << tof << " detidx " << detidx  << " bank " << evtlocation[0] 
	//	<< " moduleid " << evtlocation[1] << " xid " << evtlocation[2] << " yid " << evtlocation[3] << std::endl;
        evtinfo.insert(make_pair(tof, evtlocation));
}

void RandomHit(uint32_t* cmap, DetMapping& detmapping, TofDetMap& tofdetmap, TofVector& tofvector, EvtInfo& evtinfo){
	uint32_t tofidx = tofvector.at(uint32_t(rand()%(tofvector.size())));
	uint32_t detidx = tofdetmap[tofidx].at((rand()%(tofdetmap[tofidx].size())));
	ReadOut(detmapping, tofidx, detidx, evtinfo);
	cmap[MapIdx(tofidx, detidx)]--;
	if(cmap[MapIdx(tofidx, detidx)]==0){
		tofdetmap[tofidx].erase(std::find(tofdetmap[tofidx].begin(), tofdetmap[tofidx].end(), detidx));
		if(tofdetmap[tofidx].size()==0){
			tofdetmap.erase(tofidx);
			tofvector.erase(std::find(tofvector.begin(), tofvector.end(), tofidx));
			std::cout << "After  Remove TOF VECTOR: " << tofidx << " size " << tofvector.size() << std::endl;
		}
	}
}

void LoadMappingFile(DetMapping& detmapping, std::string mappingfilename){
	std::cout << "LoadMappingFile "<< std::endl;
	std::ifstream mappingfile(mappingfilename.c_str());
	string mappingbuff;
	getline(mappingfile, mappingbuff);


	for(uint32_t detidx=0; detidx<MAX_DET; detidx++){
		getline(mappingfile, mappingbuff);
		vector<string> substring;
		vector<uint32_t> det ;
		boost::split( substring, mappingbuff, boost::is_any_of( ";" ), boost::token_compress_on );
		for(uint32_t i=0;i<5;i++){
			uint32_t pixelid, bankid, moduleid, xid, yid;
			pixelid = atoi(substring[0].c_str());
			bankid  = atoi(substring[1].c_str());
			moduleid= atoi(substring[2].c_str());
			xid     = atoi(substring[3].c_str());
			yid     = atoi(substring[4].c_str());

			detmapping.insert(make_pair(pixelid, std::vector<uint32_t>()));
			detmapping[pixelid].push_back(bankid);
			detmapping[pixelid].push_back(moduleid);
			detmapping[pixelid].push_back(xid);
			detmapping[pixelid].push_back(yid);
			//std::cout << "Pixelid " << pixelid << " bank " << bankid 
			//	<< " moduleid " << moduleid << " xid " << xid << " yid " << yid << std::endl;

		}
	}

}

void PrintEvtInfo(EvtInfo& evtinfo){
 multimap<uint32_t, EvtLocation>::iterator it;
for(it=evtinfo.begin(); it!=evtinfo.end();++it){
EvtLocation evtlocation = (*it).second;
cout << " Tof: " << (*it).first << " Bank: " << evtlocation[0] << " Module: " << evtlocation[1]
<< " Xid: " << evtlocation[2] << " Yid: " << evtlocation[3] << endl;
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
	std::string mappingfile  ; 
	samplefile  = fConfig->pString("samplefile") ;  
	mappingfile  = fConfig->pString("mappingfile") ;  
	std::cout << "read sample file:  " << samplefile << std::endl;
	std::cout << "read mapping file: " << mappingfile << std::endl;

	DetMapping* detmapping = new DetMapping;
	LoadMappingFile(*detmapping, mappingfile); 

	uint32_t *cmap = new uint32_t[MAX_TOF*MAX_DET];
	LoadSimulationFile(cmap, samplefile); 

	TofDetMap* tofdetmap = new TofDetMap;
	TofVector* tofvector = new TofVector;
	CreateTofDetMap(cmap, *tofdetmap, *tofvector);

	srand(time(NULL));
	uint64_t hitcounts=0;
	uint64_t pulsecounts=0;
        EvtInfo* evtinfo = new EvtInfo;
	while(1){
		int hitsinpulse = rand()%20 + 10;
                evtinfo->clear();
	        std::cout << "Pulse: " << pulsecounts  << std::endl;
		for(int i=0; i< hitsinpulse; i++){
			RandomHit(cmap, *detmapping, *tofdetmap, *tofvector, *evtinfo);
			hitcounts++;
			if (tofvector->size() == 0)break;
		}
                PrintEvtInfo(*evtinfo);
		pulsecounts++;
		if (tofvector->size() == 0)break;
	}

	std::cout << "Total Events: " << hitcounts << std::endl;
	std::cout << "Total Pulse: "  << pulsecounts << std::endl;
}
