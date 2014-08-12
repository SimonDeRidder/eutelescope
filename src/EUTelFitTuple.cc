
// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Revised by Simon De Ridder
// Date 2014.08.12

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// this processor is built only if USE_AIDA and USE_GEAR are defined
#if defined(USE_GEAR) && ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )

// eutelescope inlcudes
#include "EUTelFitTuple.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
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

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static 		members mainly used to name histograms
std::string EUTelFitTuple::_FitTupleName  = "EUFit";


EUTelFitTuple::EUTelFitTuple() : Processor("EUTelFitTuple") {

  // modify processor description
  _description = "Prepare n-tuple with track fit results" ;


  // register steering parameters:
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACK,
                           "InputCollectionName" ,
                           "Name of the input Track collection"  ,
                           _inputColName ,
                           std::string("telescopetracks") ) ;

  // other processor parameters:

  registerProcessorParameter ("MissingValue",
                              "Value used for missing measurements",
                              _missingValue,  static_cast < double > (-100.));


  registerProcessorParameter ("UseManualDUT",
                              "Flag for manual DUT selection",
                              _useManualDUT,  static_cast < bool > (false));

  registerProcessorParameter ("DUTid",
                              "Id of sensor layer which should be used as DUT",
                              _DUTid,  static_cast < int > (0));

  registerProcessorParameter ("DistMax",
                              "Maximum allowed distance between fit and matched DUT hit",
                              _distMax,  static_cast < double > (0.1));


  std::vector<float > initAlign;
  initAlign.push_back(0.);
  initAlign.push_back(0.);
  initAlign.push_back(0.);

  registerProcessorParameter ("DUTalignment",
                              "Alignment corrections for DUT: shift in X, Y and rotation around Z",
                              _DUTalign, initAlign);

}


void EUTelFitTuple::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _tluTimeStamp = 0 ;


  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    message<ERROR5> ( "The GearMgr is not available, for an unknown reason." );
    exit(-1);
  }

  // Read geometry information from GEAR

  message<MESSAGE5> ( log() << "Reading telescope geometry description from GEAR ") ;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));


// Take all layers defined in GEAR geometry
  _nPlanes = _siPlanesLayerLayout->getNLayers();

// Check for DUT

  if( _siPlanesParameters->getSiPlanesType()==_siPlanesParameters->TelescopeWithDUT )
    {
      _iDUT = _nPlanes ;
      _nPlanes++;
    }
  else
    _iDUT = -1 ;

// Read position in Z (for sorting)

  _planeSort = new int[_nPlanes];
  _planePosition   = new double[_nPlanes];

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

  // Binary sorting

  bool sorted;
  do{
    sorted=false;
    for(int iz=0; iz<_nPlanes-1 ; iz++)
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

  }while(sorted);

// Book local geometry arrays

  _planeID         = new vector<int>(_nPlanes);
  _isActive        = new bool[_nPlanes];

// Fill remaining layer parameters

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

  // Get new DUT position (after sorting)

  for(int iz=0;iz< _nPlanes ; iz++)
    if(_planeSort[iz]==_iDUT)
      {
        _iDUT=iz;
        break;
      }

  // DUT position can be changed by processor parameter

  if(_useManualDUT)
    {
      bool _manualOK=false;

      for(int iz=0; iz < _nPlanes ; iz++)
        if(_planeID->at(iz)==_DUTid)
          {
            _iDUT=iz;
            _manualOK=true;
          }

      if(!_manualOK)
        {
          message<ERROR5> ( log() << "Manual DUT flag set, layer not found ID = "
                           << _DUTid
                           << "\n Program will terminate! Correct geometry description!");
          exit(-1);
        }
    }

  // Print out geometry information

  message<MESSAGE5> ( log() << "Telescope configuration with " << _nPlanes << " planes" );


  for(int ipl=0; ipl < _nPlanes; ipl++)
    {
      stringstream ss ;
      if(ipl == _iDUT)
        ss << "D.U.T.  plane" ;
      else
        if(_isActive[ipl])
          ss << "Active  plane" ;
        else
          ss << "Passive plane" ;

      ss << "  ID = " << _planeID->at(ipl)
         << "  at Z [mm] = " << _planePosition[ipl];

      message<MESSAGE5> ( log() << ss.str() );
    }


  // Allocate arrays for track fitting

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


	// Book histograms
  	bookHistos();

	// initialise EUTelGeometryTelescopeGeoDescription
    std::string name("EUTelFitTupleGeo.root");
	geo::gGeometry().initializeTGeoDescription(name,false);
}

void EUTelFitTuple::processRunHeader( LCRunHeader* runHeader) {

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  _runNr = runHeader->getRunNumber();

  message<MESSAGE5> ( log() << "Processing run header " << _nRun
                     << ", run nr " << _runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  message<MESSAGE5> ( log() << detectorName << " : " << detectorDescription ) ;

  int nDet = subDets->size();

  if(nDet)message<MESSAGE5> ( log() << nDet << " subdetectors defined :" );
  stringstream ss;
  for(int idet=0;idet<nDet;idet++)  message<MESSAGE5> (log()  << idet+1 << " : " << subDets->at(idet) );


}

void EUTelFitTuple::processEvent( LCEvent * event )
{
	EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
	if ( euEvent->getEventType() == kEORE )
	{
		message<DEBUG5> ( "EORE found: nothing else to do." );
		return;
	}

	_nEvt ++ ;
	_evtNr        = event->getEventNumber();
	_tluTimeStamp = static_cast<long int> (event->getTimeStamp());

	LCCollection* col;
	try
	{
		col = event->getCollection( _inputColName ) ;
	} catch (lcio::DataNotAvailableException& e)
	{
    	streamlog_out( DEBUG5 ) << "Not able to get collection " << _inputColName << "from event "
								<< event->getEventNumber() << " in run " << event->getRunNumber() << endl;
    	return;
  	}

	//get hitplane number of DUT
	unsigned int dutHitPlane = std::distance(_planeID->begin(),find(_planeID->begin(), _planeID->end(), _DUTid));

	// Loop over tracks in input collections

	int nTrack = col->getNumberOfElements()  ;

	message<DEBUG5> ( log() << "Total of " << nTrack << " tracks in input collection " );

	for(int itrack=0; itrack< nTrack ; itrack++)
    {
		Track * fittrack = dynamic_cast<Track*>( col->getElementAt(itrack) ) ;

		// Hit list assigned to track

		std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

		// Copy hits assign to the track to local table
		// Assign hits to sensor planes


		int nHit =   trackhits.size();

		message<DEBUG5> ( log() << "Track " << itrack << " with " << nHit << " hits, Chi2 = "
                                << fittrack->getChi2() << "/" << fittrack->getNdf());


		// Clear plane tables

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
		
		// Clear DUT variables
		double dutXLocal=_missingValue;
		double dutYLocal=_missingValue;
		double dutZLocal=_missingValue;
		double dutTx=_missingValue;
		double dutTy=_missingValue;
		double dutQ=_missingValue;
		int dutClusterSizeX=_missingValue;
		int dutClusterSizeY=_missingValue;
		bool dutHitFound = false;
		

		// setup cellIdDecoder to decode the hit properties
		CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

		// Loop over hits and fill hit tables
		for(int ihit=0; ihit< nHit ; ihit++)
        {
			TrackerHit * measHit = trackhits.at(ihit);

			// Hit position
			const double * pos = measHit->getPosition();
			
			// find plane number of the hit
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
		    {
				// Measured hits
				_isMeasured[hitPlane]=true;

				_measuredX[hitPlane]=pos[0];
				_measuredY[hitPlane]=pos[1];
				_measuredZ[hitPlane]=pos[2];
				_measuredXLocal[hitPlane]=posLocal[0];
				_measuredYLocal[hitPlane]=posLocal[1];
				_measuredZLocal[hitPlane]=posLocal[2];

				// Get cluster charge
				_measuredQ[hitPlane]=0.;

				EVENT::LCObjectVec rawdata =  measHit->getRawHits();

				if(rawdata.size()>0 && rawdata.at(0)!=NULL )
		        {
					EUTelVirtualCluster * cluster = new EUTelFFClusterImpl ( static_cast<TrackerDataImpl*> (rawdata.at(0))) ;
					_measuredQ[hitPlane]=cluster->getTotalCharge();
					if (hitPlane==dutHitPlane)
					{//save DUT cluster information
						dutHitFound = true;
						dutQ = _measuredQ[hitPlane];
						cluster->getClusterSize(dutClusterSizeX,dutClusterSizeY);
					}
		        }
				message<DEBUG5> ( log() << "Measured hit in plane " << hitPlane << " at  X = "
		                                << pos[0] << ", Y = " << pos[1] << ", Q = " << _measuredQ[hitPlane] );
			}
			else
		    {
				// Fitted hits

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

		//check values of hits for missing fitted hits (probably too strict)
		bool telPlanesOk = true;
		for (int ipl = 0; ipl<_nPlanes; ipl++)
		{
			telPlanesOk = telPlanesOk && _isMeasured[ipl] && _isFitted[ipl];
		}

		//continue if all hits have been found
		if (dutHitFound && telPlanesOk)
		{
			//copy local hits to dut variables(technically redundant)
			dutXLocal = _measuredXLocal[dutHitPlane];
			dutYLocal = _measuredYLocal[dutHitPlane];
			dutZLocal = _measuredZLocal[dutHitPlane];
			//calculate dutTx and dutTy
			double hitBeforeDUT[] = {_fittedX[dutHitPlane-1], _fittedY[dutHitPlane-1], _fittedZ[dutHitPlane-1]};
			double hitBeforeDUTLocal[3];// hit on plane before DUT in DUT local frame
			geo::gGeometry().master2Localtwo(_DUTid,hitBeforeDUT,hitBeforeDUTLocal);
			dutTx = atan2(hitBeforeDUTLocal[0]-_fittedXLocal[dutHitPlane], _fittedZLocal[dutHitPlane]-hitBeforeDUTLocal[2]);
			dutTy = atan2(hitBeforeDUTLocal[1]-_fittedYLocal[dutHitPlane], _fittedZLocal[dutHitPlane]-hitBeforeDUTLocal[2]);
			// Fill n-tuple
			int icol=0;
			_FitTuple->fill(icol++,_nEvt);
			_FitTuple->fill(icol++,_runNr);
			_FitTuple->fill(icol++,_evtNr);
			_FitTuple->fill(icol++,_tluTimeStamp);
			_FitTuple->fill(icol++,nTrack);
			_FitTuple->fill(icol++,fittrack->getNdf());
			_FitTuple->fill(icol++,fittrack->getChi2());

			for(int ipl=0; ipl<_nPlanes;ipl++)
			{
				_FitTuple->fill(icol++,_measuredXLocal[ipl]-_fittedXLocal[ipl]);
				_FitTuple->fill(icol++,_measuredYLocal[ipl]-_fittedYLocal[ipl]);
				_FitTuple->fill(icol++,_measuredZLocal[ipl]-_fittedZLocal[ipl]);
				_FitTuple->fill(icol++,_measuredQ[ipl]);
			}
			_FitTuple->fill(icol++,dutXLocal);
			_FitTuple->fill(icol++,dutYLocal);
			_FitTuple->fill(icol++,dutZLocal);
			_FitTuple->fill(icol++,dutTx);
			_FitTuple->fill(icol++,dutTy);
			_FitTuple->fill(icol++,dutQ);
			_FitTuple->fill(icol++,dutClusterSizeX);
			_FitTuple->fill(icol++,dutClusterSizeY);

			_FitTuple->addRow();
		}//end of if(dutHitFound)
    }// End of loop over tracks
	return;
}



void EUTelFitTuple::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelFitTuple::end(){

  //   std::cout << "EUTelFitTuple::end()  " << name()
  //        << " processed " << _nEvt << " events in " << _nRun << " runs "
  //        << std::endl ;


  message<MESSAGE5> ( log() << "N-tuple with "
                     << _FitTuple->rows() << " rows created" );


  // Clean memory

  delete [] _planeSort ;
  delete [] _planePosition ;
  delete _planeID ;
  delete [] _isActive ;

  delete [] _isMeasured ;
  delete [] _isFitted ;
  delete [] _measuredX  ;
  delete [] _measuredY  ;
  delete [] _measuredZ  ;
  delete [] _measuredQ  ;
  delete [] _fittedX ;
  delete [] _fittedY ;
  delete [] _fittedZ ;
  delete [] _measuredXLocal  ;
  delete [] _measuredYLocal  ;
  delete [] _measuredZLocal  ;
  delete [] _fittedXLocal ;
  delete [] _fittedYLocal ;
  delete [] _fittedZLocal ;


}



void EUTelFitTuple::bookHistos()
{


  message<MESSAGE5> ( log() << "Booking fit n-tuple \n" );

  std::vector<std::string> _columnNames;
  std::vector<std::string> _columnType;

  _columnNames.push_back("Event");
  _columnType.push_back("int");

  _columnNames.push_back("RunNr");
  _columnType.push_back("int");

  _columnNames.push_back("EvtNr");
  _columnType.push_back("int");

  _columnNames.push_back("TLUtime");
  _columnType.push_back("long int");

  _columnNames.push_back("Track");
  _columnType.push_back("int");

  _columnNames.push_back("Ndf");
  _columnType.push_back("int");

  _columnNames.push_back("Chi2");
  _columnType.push_back("float");

  const char * _varName[] = { "resX", "resY", "resZ", "measQ" };

  for(int ipl=0; ipl<_nPlanes;ipl++)
  {
    for(int ivar=0; ivar<4;ivar++)
      {
        stringstream ss;
        ss << _varName[ivar] << "_" << _planeID->at(ipl);
        _columnNames.push_back(ss.str());
        _columnType.push_back("double");
      }
  }

  // DUT variables

  _columnNames.push_back("dutXLocal");
  _columnType.push_back("double");

  _columnNames.push_back("dutYLocal");
  _columnType.push_back("double");

  _columnNames.push_back("dutZLocal");
  _columnType.push_back("double");

  _columnNames.push_back("dutTx");
  _columnType.push_back("double");

  _columnNames.push_back("dutTy");
  _columnType.push_back("double");

  _columnNames.push_back("dutQ");
  _columnType.push_back("double");

  _columnNames.push_back("dutClusterSizeX");
  _columnType.push_back("int");

  _columnNames.push_back("dutClusterSizeY");
  _columnType.push_back("int");

  _FitTuple=AIDAProcessor::tupleFactory(this)->create(_FitTupleName, _FitTupleName, _columnNames, _columnType, "");


  message<DEBUG5> ( log() << "Booking completed \n\n");

  return;
}

#endif // GEAR && AIDA
