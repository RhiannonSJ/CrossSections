/*
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : February 2017
 *
 *
*/

#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include "TNtuple.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <Math/GSLRndmEngines.h>
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TObjArray.h"

// -------------------------------------------------------------------------
// Smearing function:
// For the 2D histogram, smear the cos theta bins by 5 deg and the Tmu bins
// by 10%
// -------------------------------------------------------------------------

void Smear ( TTree   *tree,
             TNtuple *info_s,
             TNtuple *info_bg,
             TH2D    *h_unsmeared, 
             TH2D    *h_cut, 
             TH2D    *h_smeared ); 

// -------------------------------------------------------------------------
// Characterisation function:
// For the 2D histogram, draw 1D histograms slice-by-slice ( bin-by-bin )
// and characterise the signal and background
// -------------------------------------------------------------------------

void Characterisation ( TH2D    *h_smeared,
                        TNtuple *info_s,
                        TNtuple *info_bg ); 

// -------------------------------------------------------------------------
// Efficieny function:
// Loop over x and y bins and calculate the efficiencies and purities on a 
// bin-by-bin basis and plot a 2D efficiency graph and a 2D purity graph
// -------------------------------------------------------------------------

void Efficiency ( TH2D    *h_smeared, 
                  TH2D    *h_unsmeared, 
                  TH2D    *h_cut,
                  TNtuple *info_s ); 

