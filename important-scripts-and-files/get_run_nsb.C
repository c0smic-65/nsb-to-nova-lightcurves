// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sys/stat.h>

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

void get_run_nsb(int ntel, int nrun, float threshold);

void get_run_nsb(int ntel, int nrun, float threshold) {
  SashFile::HandlerC han("");

  if (!han.ConnectFile("NSB",nrun)) {
    std::cerr << "Can't find the NSB file for run " << nrun << std::endl;
    return;
  }

  Sash::DataSet *ds_run = han.GetDataSet("run");
  if (!ds_run) {
    std::cerr << "Can't find the DataSet [run]" << std::endl;
    return;
  }

  ds_run->GetEntry(0);
  Sash::HESSArray *fhess = &Sash::HESSArray::GetHESSArray(); 
    
  const Sash::RunHeader *runh = fhess->Get<Sash::RunHeader>();

  // Loop on telescopes in run, get only the ones that are available
  const Sash::PointerSet<Sash::Telescope>& telsinrun = runh->GetTelsInRun();
  for( Sash::PointerSet<Sash::Telescope>::iterator tel = telsinrun.begin(); tel != telsinrun.end(); ++tel) {
    int telId = (*tel)->GetId();
    if (telId != ntel){  //REMOVE IF YOU LOOK IN CT1-4 AS WELL
      continue;
    }
    // Distinguish between Pedestal NSB for CT1-4 and FlashCam for CT5
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

    Sash::DataSet *ds_nsb = han.GetDataSet(oss_ds_name.str());
    if (!ds_nsb) {
      std::cerr << "Can't find the DataSet [" << oss_ds_name.str() << "]" << std::endl;
      continue;
    }

    // Get number of pixels with nsb entries and create empty array
    const int nevents = (int)ds_nsb->GetEntries();
    int pix_count = 0;
    for (Sash::Pointer<Sash::Pixel> pix = (*tel)->beginPixel(), end_pix = (*tel)->endPixel(); pix != end_pix; ++pix) {
      ++pix_count;
    }
    const int size = pix_count;

    float pix_cur[size]; //last saved value of pixels

    // Iterate over entries
    for ( int ie = 0 ; ie < ds_nsb->GetEntries(); ++ie ) {
      ds_nsb->GetEntry( ie ); // Load the entry in memory

      // Get the Telescope NSB information
      const Calibration::TelescopeNSB *TelNSB = (*tel)->Get<Calibration::TelescopeNSB>(class_name.c_str());
      if (!TelNSB) {
        std::cerr << "Can't find the Calibration::TelescopeNSB class named " << class_name << " for CT" << telId << std::endl;
        continue;
      }
      TelNSB->LoadAllMembers(); 
      Sash::Time teltime_nsb = TelNSB->GetTimeStamp();
      pix_count = 0;
      for (Sash::Pointer<Sash::Pixel> pix = (*tel)->beginPixel(), end_pix = (*tel)->endPixel(); pix != end_pix; ++pix) {
       const Calibration::PixelNSB *PixNSB = TelNSB->GetMonitor(pix);
        if (!PixNSB) {
          std::cerr << "Can't find the PixelNSB class for CT" << telId << " and Pixel: "  << pix->GetConfig()->GetPixelID() << std::endl;
          continue;
        }
        if ( teltime_nsb.GetTime() == 0){
          continue;
        }
        if ( PixNSB->GetFlag() != 0 ) {
          pix_cur[pix_count] = PixNSB->GetNSBValue();
          if (pix_cur[pix_count]>threshold) {
            std::cout << int(pix_count) << ";" << pix_cur[pix_count]<< ";"<< teltime_nsb.GetUTC()<< std::endl;
          } 
        }
    ++pix_count;
      }
    }
    pix_count = 0;
  }
}