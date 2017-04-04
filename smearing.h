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
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <Math/GSLRndmEngines.h>
#include "TLatex.h"
#include "TStyle.h"
#include "TPad.h"
#include "TColor.h"
#include "TObjArray.h"

// -------------------------------------------------------------------------
// Smearing function:
// For the 2D histogram, smear the cos theta bins by 5 deg and the Tmu bins
// by 10%
// -------------------------------------------------------------------------

void Smear ( TTree *tree, 
             TH2D *h_unsmeared, 
             TH2D *h_smeared ); 

// -------------------------------------------------------------------------
// Slices function:
// For the 2D histogram, draw 1D histograms slice-by-slice ( bin-by-bin )
// -------------------------------------------------------------------------

void Slices ( TH2D *h_unsmeared, 
              TH2D *h_smeared );

// -------------------------------------------------------------------------
// Characterisation function:
// For the 2D histogram, draw 1D histograms slice-by-slice ( bin-by-bin )
// and characterise the signal and background
// -------------------------------------------------------------------------

void Characterisation ( TH2D *h_smeared,
                        TTree * tree ); 

