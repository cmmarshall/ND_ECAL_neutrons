#include "TFile.h"
#include "TTree.h"
#include <stdio.h>
#include <iostream>

int main( int argc, char const *argv[] ) 
{

  std::string topdir = "/pnfs/dune/persistent/users/marshalc/neutronSim";
  int first = 0;
  int last = 0;
  std::string horn = "FHC";
  std::string neutrino = "neutrino";
  std::string geometry = "DetEnclosure";
  bool grid = false;

  int i = 0;
  while( i < argc ) {
    if( argv[i] == std::string("--topdir") ) {
      topdir = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--geom") ) {
      geometry = argv[i+1];
      i += 2;
    } else if( argv[i] == std::string("--first") ) {
      first = atoi(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--last") ) {
      last = atoi(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--rhc") ) {
      horn = "RHC";
      neutrino = "antineutrino";
      i += 1;
    } else if( argv[i] == std::string("--grid") ) {
      grid = true;
      i += 1;
    } else i++;
  }

  double pot = 0.;
  for( int run = first; run <= last; ++run ) {
    TFile * tf;
    if( grid ) tf =  new TFile( Form("%s.%d.ghep.root",neutrino.c_str(),run) );
    else       tf  = new TFile( Form("%s/GENIE/%s/%s/%s.%d.ghep.root",topdir.c_str(),horn.c_str(),geometry.c_str(),neutrino.c_str(),run) );
    if( tf ) {
      TTree * gtree = (TTree*) tf->Get( "gtree" );
      if( gtree ) pot += gtree->GetWeight();
      tf->Close();
    }
    delete tf;
  }

  std::cout << pot << std::endl;
}
