// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sys/stat.h>
// ROOT


// HESS
#include <sash/HESSArray.hh>
#include <sash/Telescope.hh>
#include <sash/TelescopeConfig.hh>
#include <sash/Pixel.hh>
#include <sash/ListIterator.hh>
#include <sash/List.hh>
#include <sash/PixelConfig.hh>
#include <sash/RunHeader.hh>
#include <sash/DataSet.hh>
#include <sash/Pointer.hh>
#include <sash/PointerSet.hh>
#include <sash/PointerSetIterator.hh>

#include <calibration/TelescopeNSB.hh>
#include <calibration/PixelNSB.hh>

#include <utilities/SimpleStats.hh>

#include <sashfile/HandlerC.hh>
#include <stash/Coordinate.hh>
void highpixel_eval_UTC_single(int nrun);

void highpixel_eval_UTC_single(int nrun) {
  stringstream ss_high_pix;
  string home_path = "/lfs/l7/hess/users/sghosh/";

  int firstrun_in_folder = nrun -nrun % 200;
  std::cout << firstrun_in_folder << std::endl;

  stringstream make_eval_folder;
  string eval_folder = "eval_high_data_utc";
  make_eval_folder << "mkdir -p /lfs/l7/hess/users/sghosh/" << eval_folder << "/";
  std::cout <<  make_eval_folder.str() << std::endl;
  const int dir_make_eval_folder = system( make_eval_folder.str().c_str());
  if (dir_make_eval_folder == -1) {
    printf("Error in creating evaluation directory!");
    exit(1);
  }

  int chdir_eval_folder = chdir(eval_folder.c_str());


  SashFile::HandlerC han("");


  if (!han.ConnectFile("NSB",nrun)) {
    std::cout << "Can't find the NSB file for run " << nrun << std::endl;
    return;
  }

  Sash::DataSet *ds_run = han.GetDataSet("run");
  if (!ds_run) {
    std::cout << "Can't find the DataSet [run]" << std::endl;
    return;
  }


  ds_run->GetEntry(0); // Load the first entry in memory. This contain the information about the run and the list of telescope that where involved.

  Sash::HESSArray *fhess = &Sash::HESSArray::GetHESSArray(); // Retrieve the singleton class that is the access point to every other storage class

  const Sash::RunHeader *runh = fhess->Get<Sash::RunHeader>();

  //Get zenith of run:
  const Stash::Coordinate &coord =  runh->GetObservationPosition().GetCoordinate(*fhess->GetGroundSystem());
  Double_t alt = coord.GetBeta().GetDegrees();
  Double_t az  = coord.GetLambda().GetDegrees();




  //Create directories


  ofstream high_pix;
  stringstream nrun_high_pix;
  stringstream run_numbers;
  run_numbers << "run"<<firstrun_in_folder <<"-"<<firstrun_in_folder + 199 <<"/";
  nrun_high_pix << "mkdir -p high_pixel/"<< run_numbers.str() << nrun;
  std::cout <<  nrun_high_pix.str() << std::endl;

  const int dir_high_pix = system( nrun_high_pix.str().c_str());

  if (dir_high_pix == -1) {
    printf("Error in creating high pixel directory!");
    exit(1);
  }
  // Now loop on the telescope in run
  const Sash::PointerSet<Sash::Telescope>& telsinrun = runh->GetTelsInRun();
  for( Sash::PointerSet<Sash::Telescope>::iterator tel = telsinrun.begin(); tel != telsinrun.end(); ++tel) {
    int telId = (*tel)->GetId();
    if (telId != 5){  //REMOVE IF YOU LOOK IN CT1-4 AS WELL
      continue;
    }
    // At the moment, only the NSB derived from Pedestal width is valid for CT1-4 and HESS2_2048 camera
    // For FlashCam, this is the one derived on event that is valid
    std::string ds_suffix;
    std::string class_name;
    if ( (*tel)->GetConfig()->GetSetupID() == Sash::TelescopeConfig::FlashCam_1764 ) {
      ds_suffix = "EvtNSB";
      class_name = "NSBEvent";
    }
    else {
      ds_suffix = "PedNSB";
      class_name = "NSBPed";
    }

    // Retrieve the DataSet

    std::ostringstream oss_ds_name;
    oss_ds_name << "CT" << telId << "_" << ds_suffix;

    std::cout << "Running for CT" << telId << std::endl;

    Sash::DataSet *ds_nsb = han.GetDataSet(oss_ds_name.str());
    if (!ds_nsb) {
      std::cout << "Can't find the DataSet [" << oss_ds_name.str() << "]" << std::endl;
      continue;
    }

    //Open and name the files
    ss_high_pix << home_path << eval_folder  << "/high_pixel/" << run_numbers.str().c_str() << nrun << "/high_pix_"<< nrun <<"_CT_"<<telId<<".txt";
    std::cout << ss_high_pix.str()<<std::endl;
    high_pix.open(ss_high_pix.str().c_str());

    const int nevents = (int)ds_nsb->GetEntries();
    int pix_count = 0;
    for (Sash::Pointer<Sash::Pixel> pix = (*tel)->beginPixel(), end_pix = (*tel)->endPixel(); pix != end_pix; ++pix) {
      ++pix_count;
    }
    const int size = pix_count;

    float pix_cur[size]; //last saved value of pixels
    float threshold = 0; //MHz, this is the threshold for high pixel; SG: Chnaged from 2000 to 0 to get all (unfiltered) data

    if (nrun < 153001) {
        threshold = 0;
    }

      // Now loop on all the dataset entry (one per pedestal evaluation)
    Utilities::SimpleStats uss_nsb;
    for ( int ie = 0 ; ie < ds_nsb->GetEntries(); ++ie ) { // For pedestal, the first entry is not necessarily a great estimate because the pedestal is not estimated on lot of events.... skip it ?
      ds_nsb->GetEntry( ie ); // Load the entry in memory


      Utilities::SimpleStats uss_nsb_evt; // new each event / time stamp

      // Get the Telescope NSB information
      const Calibration::TelescopeNSB *TelNSB = (*tel)->Get<Calibration::TelescopeNSB>(class_name.c_str()); // Get retrieve a const pointer (Handle will create the object if it does not exist in memory)
      if (!TelNSB) {
	std::cout << "Can't find the Calibration::TelescopeNSB class named " << class_name << " for CT" << telId << std::endl;
	continue;
      }
      TelNSB->LoadAllMembers(); // Not needed because TelescopeNSB contain simple data but nice to keep in mind that sometimes it is needed to force complex stored information like pointers) to be loaded in memory
      Sash::Time teltime_nsb = TelNSB->GetTimeStamp();
      pix_count = 0;
    for (Sash::Pointer<Sash::Pixel> pix = (*tel)->beginPixel(), end_pix = (*tel)->endPixel(); pix != end_pix; ++pix) {
	const Calibration::PixelNSB *PixNSB = TelNSB->GetMonitor(pix);

	if (!PixNSB) {
	  std::cout << "Can't find the PixelNSB class for CT" << telId << " and Pixel: "  << pix->GetConfig()->GetPixelID() << std::endl;
	  continue;
	}
	if ( teltime_nsb.GetTime() == 0){
	  continue;
	}
	if ( PixNSB->GetFlag() != 0 ) {
	  pix_cur[pix_count] = PixNSB->GetNSBValue();
	  
	  uss_nsb.Fill( PixNSB->GetNSBValue() );
	  uss_nsb_evt.Fill( PixNSB->GetNSBValue() );
	  
	  if (pix_cur[pix_count]>threshold) {
	    high_pix << int(pix_count) << ";" << pix_cur[pix_count]<< ";"<< teltime_nsb.GetUTC() //<< ";"
		     //<<(*pix).GetConfig()->GetPos().GetX() << ";" << (*pix).GetConfig()->GetPos().GetY()
		     << std::endl;
	  } 
	}
	++pix_count;
      }
     }
    pix_count = 0;

high_pix.close();
    ss_high_pix.str("");

    std::cout << "CT" << telId << " Average Camera NSB is " << uss_nsb.GetMean() << " MHz (RMS: " << uss_nsb.GetRMS() << ")" << std::endl;

  }
  
  
}
