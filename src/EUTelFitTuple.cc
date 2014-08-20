// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Revised by Simon De Ridder
// Date 2014.08.15

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 */

// this processor is built only if USE_AIDA and USE_GEAR are defined
#if defined(USE_GEAR) && ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )

// eutelescope inlcudes
#include "EUTelFitTuple.h"
#include "EUTelVirtualCluster.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelUtility.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITupleFactory.h>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>

//root includes
#include <TTree.h>
#include <TFile.h>

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

//definition of static members

//name of NTuple
std::string EUTelFitTuple::_TreeName  = "EUFit";

EUTelFitTuple::EUTelFitTuple() : Processor("EUTelFitTuple")
{
	//modify processor description
	_description = "Prepare n-tuple with track fit results" ;


	//register steering parameters:
	//	name, description, class-variable, default value

	//input collection name first:
	registerInputCollection(LCIO::TRACK,
							"InputCollectionName",
							"Name of the input Track collection",
							_inputColName,
							std::string("telescopetracks") ) ;

	//other processor parameters:
	registerProcessorParameter("MissingValue",
								"Value used for missing measurements",
 								_missingValue, static_cast<double> (-100.));

	registerProcessorParameter("UseManualDUT",
								"Flag for manual DUT selection",
								_useManualDUT, static_cast<bool> (false));

	registerProcessorParameter("DUTid",
								"Id of sensor layer which should be used as DUT",
								_DUTid, static_cast<int> (0));

 	registerProcessorParameter ("OutputPath",
								"Path/File where root-file should be stored",
								_path2file, std::string("NTuple.root"));
}

void EUTelFitTuple::init()
{
	//usually a good idea to
	printParameters() ;

	//initialise run number, event number and TLU timestamp
	_nRun = 0 ;
	_nEvt = 0 ;
	_tluTimeStamp = 0 ;

	//check if the GEAR manager pointer is not null!
	if (Global::GEAR == 0x0)
	{
		message<ERROR5> ( "The GearMgr is not available, for an unknown reason." );
		exit(-1);
	}

	//Read geometry information from GEAR
	message<MESSAGE5> ( log() << "Reading telescope geometry description from GEAR ") ;
	_siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
	_siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

	//Take all layers defined in GEAR geometry
	_nPlanes = _siPlanesLayerLayout->getNLayers();

	//Check for DUT
	if( _siPlanesParameters->getSiPlanesType()==_siPlanesParameters->TelescopeWithDUT )
	{
		_iDUT = _nPlanes;
		_nPlanes++;
	} else
	{
		_iDUT = -1;
	}

	//Read position in Z (for sorting)
	_planeSort = new int[_nPlanes];
	_planePosition = new double[_nPlanes];
	for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)
	{
		_planePosition[ipl]=_siPlanesLayerLayout->getLayerPositionZ(ipl);
		_planeSort[ipl]=ipl;
	}
	if(_iDUT>0)
	{
		_planePosition[_iDUT]=_siPlanesLayerLayout->getDUTPositionZ();
		_planeSort[_iDUT]=_iDUT;
	}

	//Sorting of the planes
	bool sorted;
	do{
		sorted=false;
		for(int iz=0; iz<_nPlanes-1 ; iz++)
		{
			if(_planePosition[iz]>_planePosition[iz+1])
			{
				double _posZ = _planePosition[iz];
				_planePosition[iz] = _planePosition[iz+1];
				_planePosition[iz+1] = _posZ;

				int _idZ = _planeSort[iz];
				_planeSort[iz] = _planeSort[iz+1];
				_planeSort[iz+1] = _idZ;

				sorted=true;
			}
		}
	}while(sorted);

	//initialise remaining layer geometry collections
	_planeID         = new vector<int>(_nPlanes);
	_isActive        = new bool[_nPlanes];

  	//Fill remaining layer geometry collections
	for(int iz=0; iz < _nPlanes ; iz++)
    {
		int ipl=_planeSort[iz];
		double resolution;
		if(ipl != _iDUT )
        {
			_planeID->at(ipl) = (_siPlanesLayerLayout->getID(ipl));
			resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);
		}
		else
		{
			_planeID->at(ipl) = (_siPlanesLayerLayout->getID(ipl));
			resolution = _siPlanesLayerLayout->getDUTSensitiveResolution();
		}
		_isActive[iz] = (resolution > 0);
	}

	//Get new DUT position (after sorting)
	for(int iz=0;iz< _nPlanes ; iz++)
	{
		if(_planeSort[iz]==_iDUT)
		{
			_iDUT=iz;
			break;
		}
	}

	//DUT position can be changed by processor parameter
	if(_useManualDUT)
    {
		bool _manualOK=false;
		for(int iz=0; iz < _nPlanes ; iz++)
		{
			if(_planeID->at(iz)==_DUTid)
			{
				_iDUT=iz;
				_manualOK=true;
			}
		}
		if(!_manualOK)
        {
			message<ERROR5> (log() << "Manual DUT flag set, layer not found ID = "
									<< _DUTid
									<< "\n Program will terminate! Correct geometry description!");
			exit(-1);
        }
    }

	//Print out geometry information
	message<MESSAGE5> (log() << "Telescope configuration with " << _nPlanes << " planes");
	for(int ipl=0; ipl < _nPlanes; ipl++)
    {
		stringstream ss;
		if(ipl == _iDUT)
		{
		    ss << "D.U.T.  plane";
		} else
		{
			if(_isActive[ipl])
			{
				ss << "Active  plane";
			} else
			{
				ss << "Passive plane";
			}
		}
		ss << "  ID = " << _planeID->at(ipl)
			<< "  at Z [mm] = " << _planePosition[ipl];
		message<MESSAGE5> ( log() << ss.str() );
    }


	//Allocate arrays for track parameters
	_isMeasured     = new bool[_nPlanes];
	_isFitted       = new bool[_nPlanes];

	_measuredX      = new double[_nPlanes];
	_measuredY      = new double[_nPlanes];
	_measuredZ      = new double[_nPlanes];
	_measuredQ      = new double[_nPlanes];
	_fittedX        = new double[_nPlanes];
	_fittedY        = new double[_nPlanes];
	_fittedZ        = new double[_nPlanes];
	_measuredXLocal = new double[_nPlanes];
	_measuredYLocal = new double[_nPlanes];
	_measuredZLocal = new double[_nPlanes];
	_fittedXLocal   = new double[_nPlanes];
	_fittedYLocal   = new double[_nPlanes];
	_fittedZLocal   = new double[_nPlanes];
	_resXLocal		= new double[_nPlanes];
	_resYLocal		= new double[_nPlanes];
	_resZLocal		= new double[_nPlanes];

	_p_col = new std::vector<int>();
	_p_row = new std::vector<int>();
	_p_tot = new std::vector<int>();
	_p_lv1 = new std::vector<int>();

	//Book histograms
  	bookHistos();

	//initialise EUTelGeometryTelescopeGeoDescription (for coordinate conversion)
    std::string name("EUTelFitTupleGeo.root");
	geo::gGeometry().initializeTGeoDescription(name,false);
}

void EUTelFitTuple::processRunHeader(LCRunHeader* runHeader)
{
	auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
	eutelHeader->addProcessor( type() );

	_nRun++ ;

	//Decode and print out Run Header information - just a check
	_runNr = runHeader->getRunNumber();
	message<MESSAGE5> (log() << "Processing run header " << _nRun
						<< ", run nr " << _runNr);

	const std::string detectorName = runHeader->getDetectorName();
	const std::string detectorDescription = runHeader->getDescription();
	const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();
	message<MESSAGE5> ( log() << detectorName << " : " << detectorDescription ) ;

	int nDet = subDets->size();
	if(nDet)
	{
		message<MESSAGE5> (log() << nDet << " subdetectors defined :");
	}
	stringstream ss;
	for(int idet=0;idet<nDet;idet++)
	{
		message<MESSAGE5> (log()  << idet+1 << " : " << subDets->at(idet));
	}
}

void EUTelFitTuple::processEvent(LCEvent * event)
{
	//Decode event and get its properties
	EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
	if ( euEvent->getEventType() == kEORE )
	{
		message<DEBUG5> ( "EORE found: nothing else to do." );
		return;
	}

	_nEvt++;
	_evtNr        = event->getEventNumber();
	_tluTimeStamp = static_cast<long int> (event->getTimeStamp());

	//get the input tracker collection from the event
	LCCollection* col;
	try
	{
		col = event->getCollection( _inputColName ) ;
	} catch (lcio::DataNotAvailableException& e)
	{
    	streamlog_out( DEBUG5 ) << "Unable to get collection " << _inputColName << " from event "
								<< event->getEventNumber() << " in run " << event->getRunNumber() << endl;
    	return;
  	}

	//get plane ID of DUT
	unsigned int dutPlaneID = std::distance(_planeID->begin(),find(_planeID->begin(), _planeID->end(), _DUTid));

	//Loop over tracks in input collections
	_nTracks = col->getNumberOfElements();
	message<DEBUG5> ( log() << "Total of " << _nTracks << " tracks in input collection " );
	for(int itrack=0; itrack< _nTracks ; itrack++)
    {
		Track * fittrack = dynamic_cast<Track*>(col->getElementAt(itrack));

		//Hit list assigned to track (measured and fitted)
		std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

		//get and print general track information
		int nHit =   trackhits.size();
		_ndf = fittrack->getNdf();
		_chi2 = fittrack->getChi2();
		message<DEBUG5> ( log() << "Track " << itrack << " with " << nHit << " hits, Chi2 = " << _chi2 << "/" << _ndf);
		
		//clear data variables
		clearVariables();
		bool dutHitFound = false;

		//setup cellIdDecoder to decode the hit properties
		CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

		//Loop over hits in this track and fill hit tables
		for(int ihit=0; ihit< nHit ; ihit++)
        {
			TrackerHit * measHit = trackhits.at(ihit);

			//get hit position
			const double * pos = measHit->getPosition();
			
			//find plane number of the hit
		    int sensorID = Utility::getSensorIDfromHit(measHit);
			unsigned int hitPlane = std::distance(_planeID->begin(),find(_planeID->begin(), _planeID->end(), sensorID));
			if (hitPlane==_planeID->size())
			{
				message<DEBUG5> ( log() << "hit not matched to plane, sensorID=" << sensorID);
				continue;
			}

			//find local hit position
			double posLocal[3];
			geo::gGeometry().master2Localtwo(sensorID,pos,posLocal);

			if( (hitCellDecoder(measHit)["properties"] & kFittedHit) == 0 )
		    {//Measured hit
				//Copy positions
				_isMeasured[hitPlane]=true;

				_measuredX[hitPlane]=pos[0];
				_measuredY[hitPlane]=pos[1];
				_measuredZ[hitPlane]=pos[2];
				_measuredXLocal[hitPlane]=posLocal[0];
				_measuredYLocal[hitPlane]=posLocal[1];
				_measuredZLocal[hitPlane]=posLocal[2];

				//save extra hit information if DUT
				if (hitPlane==dutPlaneID)
				{
					dutHitFound = true;
					//save dut hit information
					_dutXLocal = posLocal[0];
					_dutYLocal = posLocal[1];
					_dutZLocal = posLocal[2];
				}

				// Get cluster data
				_measuredQ[hitPlane]=0.;

				EVENT::LCObjectVec rawdata =  measHit->getRawHits();
				if(rawdata.size()>0 && rawdata.at(0)!=NULL )
			    {
					EUTelVirtualCluster* cluster = new EUTelSparseClusterImpl<EUTelGenericSparsePixel>(static_cast<TrackerDataImpl*>(rawdata.at(0)));
					_measuredQ[hitPlane]=cluster->getTotalCharge();

					//save extra cluster information if DUT
					if (hitPlane==dutPlaneID)
					{
						_dutQ = _measuredQ[hitPlane];
						cluster->getClusterSize(_dutClusterSizeX, _dutClusterSizeY);
						IMPL::TrackerDataImpl* clusterContent = cluster->trackerData();
						std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> >
							apixData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(clusterContent));
						EUTelGenericSparsePixel apixPixel;
						for( unsigned int iHit = 0; iHit < apixData->size(); iHit++ )
						{
							apixData->getSparsePixelAt( iHit, &apixPixel);
							_nPixHits++;
							_p_col->push_back( apixPixel.getXCoord() );
							_p_row->push_back( apixPixel.getYCoord() );
							_p_tot->push_back( static_cast< int >(apixPixel.getSignal()) );
							_p_lv1->push_back( static_cast< int >(apixPixel.getTime()) );
						}
					}
		        }
				message<DEBUG5> ( log() << "Measured hit in plane " << hitPlane << " at  X = " << posLocal[0]
										<< ", Y = " << posLocal[1] << ", Q = " << _measuredQ[hitPlane] );
			}
			else
		    {// Fitted hit
				//Copy positions
				_isFitted[hitPlane]=true;

				_fittedX[hitPlane]=pos[0];
				_fittedY[hitPlane]=pos[1];
				_fittedZ[hitPlane]=pos[2];
				_fittedXLocal[hitPlane]=posLocal[0];
				_fittedYLocal[hitPlane]=posLocal[1];
				_fittedZLocal[hitPlane]=posLocal[2];

				message<DEBUG5> ( log() << "Fitted  hit  in plane " << hitPlane << " at  X = "
		                                << pos[0] << ", Y = " << pos[1] );
		    }
        }// End of loop over hits in track

		//continue if necessary hits have been found
		if (dutHitFound && _isFitted[dutPlaneID-1])
		{
			//calculate dutTx and dutTy
			double hitBeforeDUT[] = {_fittedX[dutPlaneID-1], _fittedY[dutPlaneID-1], _fittedZ[dutPlaneID-1]};
			double hitBeforeDUTLocal[3];// hit on plane before DUT in DUT local frame
			geo::gGeometry().master2Localtwo(_DUTid,hitBeforeDUT,hitBeforeDUTLocal);
			_dutTx = atan2(hitBeforeDUTLocal[0]-_fittedXLocal[dutPlaneID], _fittedZLocal[dutPlaneID]-hitBeforeDUTLocal[2]);
			_dutTy = atan2(hitBeforeDUTLocal[1]-_fittedYLocal[dutPlaneID], _fittedZLocal[dutPlaneID]-hitBeforeDUTLocal[2]);
			//calculate residuals
			for (int ipl = 0; ipl<_nPlanes; ipl++)
			{
				_resXLocal[ipl] = _measuredXLocal[ipl]-_fittedXLocal[ipl];
				_resYLocal[ipl] = _measuredYLocal[ipl]-_fittedYLocal[ipl];
				_resZLocal[ipl] = _measuredZLocal[ipl]-_fittedZLocal[ipl];
			}

			// Fill TTree
			_euTree->Fill();
		}//end of if(dutHitFound)
    }// End of loop over tracks
	return;
}


void EUTelFitTuple::check( LCEvent * /* evt */ )
{
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelFitTuple::end()
{
	message<MESSAGE5> (log() << "TTree with "
						<< _euTree->GetEntries() << " rows created");
	//write tree
	_outputFile->Write();

	// Clean memory

	delete [] _planeSort;
	delete [] _planePosition;
	delete _planeID;
	delete [] _isActive;

	delete [] _isMeasured;
	delete [] _isFitted;
	delete [] _measuredX;
	delete [] _measuredY;
	delete [] _measuredZ;
	delete [] _measuredQ;
	delete [] _fittedX;
	delete [] _fittedY;
	delete [] _fittedZ;
	delete [] _measuredXLocal;
	delete [] _measuredYLocal;
	delete [] _measuredZLocal;
	delete [] _fittedXLocal;
	delete [] _fittedYLocal;
	delete [] _fittedZLocal;
	delete [] _resXLocal;
	delete [] _resYLocal;
	delete [] _resZLocal;

	delete _p_col;
	delete _p_row;
	delete _p_tot;
	delete _p_lv1;
}



void EUTelFitTuple::bookHistos()
{
	message<MESSAGE5> ( log() << "Booking TTree \n" );

	_outputFile = new TFile(_path2file.c_str(),"RECREATE");

	_euTree = new TTree(_TreeName.c_str(),_TreeName.c_str());
	_euTree->Branch("Event", &_nEvt);
	_euTree->Branch("RunNr", &_runNr);
	_euTree->Branch("EvtNr", &_evtNr);
	_euTree->Branch("TLUtime", &_tluTimeStamp);
	_euTree->Branch("Tracks", &_nTracks);
	_euTree->Branch("Ndf", &_ndf);
	_euTree->Branch("Chi2", &_chi2);
	for(int ipl=0; ipl<_nPlanes;ipl++)
	{
		stringstream sx;
		stringstream sy;
		stringstream sz;
		stringstream sq;
		sx << "resX_" << _planeID->at(ipl);
		sy << "resY_" << _planeID->at(ipl);
		sz << "resZ_" << _planeID->at(ipl);
		sq << "measQ_" << _planeID->at(ipl);
		_euTree->Branch(sx.str().c_str(), &_resXLocal[ipl]);
		_euTree->Branch(sy.str().c_str(), &_resYLocal[ipl]);
		_euTree->Branch(sz.str().c_str(), &_resZLocal[ipl]);
		_euTree->Branch(sq.str().c_str(), &_measuredQ[ipl]);
	}
	_euTree->Branch("dutXLocal", &_dutXLocal);
	_euTree->Branch("dutYLocal", &_dutYLocal);
	_euTree->Branch("dutZLocal", &_dutZLocal);
	_euTree->Branch("dutTx", &_dutTx);
	_euTree->Branch("dutTy", &_dutTy);
	_euTree->Branch("dutQ", &_dutQ);
	_euTree->Branch("dutClusterSizeX", &_dutClusterSizeX);
	_euTree->Branch("dutClusterSizeY", &_dutClusterSizeY);
	_euTree->Branch("dutNPix", &_nPixHits);
	_euTree->Branch("dutPixX", &_p_col);
	_euTree->Branch("dutPixY", &_p_row);
	_euTree->Branch("dutPixTOT", &_p_tot);
	_euTree->Branch("dutPixTime", &_p_lv1);

	message<DEBUG5> ( log() << "Booking completed \n\n");

	return;
}

void EUTelFitTuple::clearVariables()
{
	//Clear hit tables
		for(int ipl=0;ipl<_nPlanes;ipl++)
        {
			_isMeasured[ipl] =	false;
			_isFitted[ipl] =	false;

			_measuredX[ipl] =	_missingValue;
			_measuredY[ipl] =	_missingValue;
			_measuredZ[ipl] =	_missingValue;
			_measuredQ[ipl] =	_missingValue;

			_fittedX[ipl] =	_missingValue;
			_fittedY[ipl] =	_missingValue;
			_fittedZ[ipl] =	_missingValue;
			_measuredXLocal[ipl] =	_missingValue;
			_measuredYLocal[ipl] =	_missingValue;
			_measuredZLocal[ipl] =	_missingValue;
			_fittedXLocal[ipl] =	_missingValue;
			_fittedYLocal[ipl] =	_missingValue;
			_fittedZLocal[ipl] =	_missingValue;
        }
		
		//Clear DUT variables
		_dutTx=_missingValue;
		_dutTy=_missingValue;
		_dutQ=_missingValue;
		_dutClusterSizeX=_missingValue;
		_dutClusterSizeY=_missingValue;

		//clear pixel variables
		_nPixHits = 0;
		_p_col->clear();
		_p_row->clear();
		_p_tot->clear();
		_p_lv1->clear();
}

#endif // GEAR && AIDA
