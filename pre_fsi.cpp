/*
 *--------------------------------------------------------------
 * An exercise to smear and attempt to unfold a simple histrogram 
 *
 *--------------------------------------------------------------
 * Author: Rhiannon Jones
 * Date  : March 2017
 *
 * For now, this is a simple macro for neutrino-mode only
 * I will amend this entire system to be a class structure with 
 * constructors and destructors for neutrino and antineutrino-modes
 *
*/

#include "pre_fsi.h"

using namespace std; 

int pre_fsi() {

    // -------------------------------------------------------------------------
    //                              Open event files
    // -------------------------------------------------------------------------
    
    // Open Default
    TFile f1("/hepstore/rjones/Exercises/Flavours/Default+MEC/sbnd/1M/gntp.10000.gst.root");
    if (f1.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Default event file is open " << endl;
    }

    // -------------------------------------------------------------------------
    //                            Parameters to consider
    // -------------------------------------------------------------------------
    
    // Number of MiniBooNE bins (data) 
    // costhetamu : 20
    // Tmu        : 18
    // 
    // Amount to smear 
    // costhetamu : 5 degrees
    // Tmu        : 10%
    //
    // Kinetic energy cuts
    // mu         : > 50MeV
    // charged pi : > 50MeV
   
    // -------------------------------------------------------------------------
    //                  Get the tree, normalisation and the branches
    // -------------------------------------------------------------------------
    TTree *def_tree = (TTree*) f1.Get( "gst" );
 
    // Define the log and normal histograms for the unsmeared distributions
    // El - m_mu : kinetic energy of the final state primary lepton (fspl : muon == 13 )
    // cthl      : costheta of the final state primary lepton
    // h_un      : unsmeared histogram
    // h_sm      : smeared histogram
    //
    TCanvas *c1 = new TCanvas("c1","",800,600);
    //TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} distribution before smearing ",20,-1,1,18,0,2);
    TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} event rates, truth level ",20,-1,1,18,0,2);
    def_tree->Draw("( El - 0.10566 ):cthl>>h_un","fspl == 13 && cc && (nfpip + nfpim + nfpi0 == 0)","colz"); 
   
    c1->SetRightMargin(0.13);

    h_un->SetStats(kFALSE);
    h_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un->GetYaxis()->SetTitle("T_{#mu}");
    c1->SetLogz();
    //c1->SaveAs("distribution_before.png");
    c1->SaveAs("distribution_before_log.png");
    
    // The same histogram definitions for the smeared distributions
    //TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, impurity ",20,-1,1,18,0,2);
    TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} event rates, after smearing ",20,-1,1,18,0,2);

    // Initiate the TNtuple to hold the interesting features
    TNtuple *nt = new TNtuple("nt","Important information about signal","ev:cth:T:neu:cc:nc:npip:npim:npi0:qel:res:mec");
    nt->SetDirectory(0);

    // Take h_un and smear it
    Smear(def_tree, nt, h_un, h_sm );

    // The smeared histogram and an ntuple containing the important information
    // to help characterise the signal
    Characterisation( h_sm, nt );

    /*
    // The efficiency function 
    Efficiency( h_sm, h_un, h_c, nt_s );
    */

    delete h_sm;
    delete h_un;
    delete c1;
    delete nt;

    return 0;
}

// -------------------------------------------------------------------------
//                            Function definitions
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
//                             Smearing function
// -------------------------------------------------------------------------

void Smear ( TTree   *tree, 
             TNtuple *info,
             TH2D  *h_unsmeared,
             TH2D  *h_smeared ){

    // Get the branches
    TBranch *b_nf    = tree->GetBranch( "nf" );
    TBranch *b_neu   = tree->GetBranch( "neu" );
    TBranch *b_cthl  = tree->GetBranch( "cthl" );
    TBranch *b_El    = tree->GetBranch( "El" );
    TBranch *b_fspl  = tree->GetBranch( "fspl" );
    TBranch *b_cthf  = tree->GetBranch( "cthf" );
    TBranch *b_Ef    = tree->GetBranch( "Ef" );
    TBranch *b_pdgf  = tree->GetBranch( "pdgf" );
    TBranch *b_cc    = tree->GetBranch( "cc" );
    TBranch *b_nc    = tree->GetBranch( "nc" );
    TBranch *b_qel   = tree->GetBranch( "qel" );
    TBranch *b_res   = tree->GetBranch( "res" );
    TBranch *b_mec   = tree->GetBranch( "mec" );
    TBranch *b_nipip = tree->GetBranch( "nipip" );
    TBranch *b_nipim = tree->GetBranch( "nipim" );
    TBranch *b_nipi0 = tree->GetBranch( "nipi0" );
    TBranch *b_nfpip = tree->GetBranch( "nfpip" );
    TBranch *b_nfpim = tree->GetBranch( "nfpim" );
    TBranch *b_nfpi0 = tree->GetBranch( "nfpi0" );

    // Number of events in the TTree
    int n_values = tree->GetEntries();
    //int n_values = 1000; // Using 100 events for debugging purposes
  
    double m_mu = 0.10566; // Muon mass, GeV
    double m_pi = 0.13957; // Charged pion mass, GeV
    
    // Vectors to fill for the impurity implementation
    vector< double > T_mu_vect;
    vector< double > T_pi_vect;
    vector< double > cos_mu_vect;
    vector< double > cos_pi_vect;
  
    vector< int > Impurity;
   
    // Loop over the entries of the tree and calculate the kinetic energies 
    // of the muons and pions and define the impurity
    for ( int i = 0; i < n_values; ++i ){
        
        tree->GetEntry(i);

        int nf    = b_nf->GetLeaf("nf")->GetValue();
        int fspl  = b_fspl->GetLeaf("fspl")->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf("nfpi0")->GetValue();
        int nfpip = b_nfpip->GetLeaf("nfpip")->GetValue();
        int nfpim = b_nfpim->GetLeaf("nfpim")->GetValue();

        double e_mu   = b_El->GetLeaf("El")->GetValue();
        double cos_mu = b_cthl->GetLeaf("cthl")->GetValue();
        
        double T_mu;

        // Calculate the kinetic energy for muons
        if ( fspl == 13 ){  
         
            // Energy of the final state primary lepton
            T_mu = e_mu - m_mu;

            T_mu_vect.push_back(T_mu);
            cos_mu_vect.push_back(cos_mu);

        }
        // If the final state primary is a lepton, push back a number that will
        // be removed in the cuts later
        else if ( fspl == 11 ){
            T_mu_vect.push_back(-99999);
            cos_mu_vect.push_back(-99999);
        }
        else{
            T_mu_vect.push_back(-99999);
            cos_mu_vect.push_back(-99999);
        }

        
        if ( nfpip + nfpim == 1 && nfpi0 == 0 ){
    
            // For all the final state hadronic particles, get their pdg code
            for ( int j = 0; j < nf; ++j ) {
         
                b_pdgf->GetEntry(i);
                b_cthf->GetEntry(i);
                b_Ef->GetEntry(i);
                
                int pdgf      = b_pdgf->GetLeaf("pdgf")->GetValue(j);
                double e_pi   = b_Ef->GetLeaf("Ef")->GetValue(j);
                double cos_pi = b_cthf->GetLeaf("cthf")->GetValue(j); 
            
                
                // Calculate the kinetic energy of the pions
                if (pdgf == 211 || pdgf == -211 ){

                    double T_pi = e_pi - m_pi;

                    T_pi_vect.push_back(T_pi);
                    cos_pi_vect.push_back(cos_pi);
                }
            }
        }
        else{
            T_pi_vect.push_back(-99999); 
            cos_pi_vect.push_back(-99999); 
        }


        // Apply the impurity cut
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 && nfpip + nfpim == 1 && nfpi0 == 0 ){
            int random;
            random = rand() % 5 + 1;

            // If the random number = 5 then push back a 1 onto the vector
            // If the random number < 5 push back a 0
            if( random == 5 ){
                Impurity.push_back(1);       
            }
            else{
                Impurity.push_back(0);
            }
        }
        else{
            Impurity.push_back(0);
        }
    }

    // Vectors to hold the smeared values of T and cos theta
    vector< double > T_mu_prime;
    vector< double > cos_mu_prime;

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();    
    _random_gen->SetSeed( time( NULL ) ); 
    
    // Event by event, generate Tmu_prime and Tpi_prime: lognormal
    // Then find thetamu_prime and thetapi_prime: gaussian
    for ( int i = 0; i < n_values; ++i ){

        tree->GetEntry(i);

        double El    = b_El->GetLeaf( "El" )->GetValue();
        double cthl  = b_cthl->GetLeaf( "cthl" )->GetValue();

        // -------------------------------------------------------
        //                      Kinetic energy
        // -------------------------------------------------------
        // Calculate the mean and sigma for the LogNormal function
        //      zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
        //      sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
        //      m     = expectation value = El - m_mu
        //      var   = variance = s.d.^2 = ( El - m_mu * 0.1 ) ^ 2
    
        double var_mu     = TMath::Power( ( El - m_mu ) * 0.1, 2 ); 
        double sigma_mu   = TMath::Sqrt( TMath::Log( 1 + ( var_mu / TMath::Power( ( El - m_mu ), 2 ) ) ) );
        double zeta_mu    = TMath::Log( ( El - m_mu ) * ( 1 / TMath::Sqrt( 1 + ( var_mu / TMath::Power( ( El - m_mu ), 2 ) ) ) ) );
        double lognorm_mu = _random_gen->LogNormal( zeta_mu, sigma_mu );

        T_mu_prime.push_back(lognorm_mu);
            
        // -------------------------------------------------------
        //                      Cos theta
        // -------------------------------------------------------
        
        // Calculate the mean and sigma for the LogNormal function
        //      theta = acos(costtheta)
        //      var   = 5 degrees
    
        double sd_thetamu = TMath::Pi() / 36; // 5 degrees
        double gaus_theta = TMath::ACos( cthl ) + _random_gen->Gaussian( sd_thetamu );
        double gaus_costheta = TMath::Cos( gaus_theta ); 

        cos_mu_prime.push_back(gaus_costheta);

    } 

    // Implement the smearing factors and plot the histograms both before and after to observe the effect
    // Let's start with Gaussian smearing for now and move onto log normal smearing later
    // Now try the smearing
    TCanvas *c2 = new TCanvas("c2","",800,600);

    // Fill h_smeared normally and with the impurities
    for ( int i = 0; i < n_values; ++i ){        
        
        tree->GetEntry(i); 

        int neu   = b_neu->GetLeaf( "neu" )->GetValue();
        int cc    = b_cc->GetLeaf( "cc" )->GetValue(); 
        int nc    = b_nc->GetLeaf( "nc" )->GetValue(); 
        int qel   = b_qel->GetLeaf( "qel" )->GetValue();
        int res   = b_res->GetLeaf( "res" )->GetValue();
        int mec   = b_mec->GetLeaf( "mec" )->GetValue();
        int fspl  = b_fspl->GetLeaf( "fspl" )->GetValue(); 
        int nipip = b_nipip->GetLeaf( "nipip" )->GetValue();
        int nipim = b_nipim->GetLeaf( "nipim" )->GetValue();
        int nipi0 = b_nipi0->GetLeaf( "nipi0" )->GetValue();
        int nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        int nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();

        // Set the energy and impurity cuts
        // Filling the signal ntuple
        if ( fspl == 13 
             && cc
             && ( nfpip + nfpim + nfpi0 == 0 )
             && T_mu_vect[i] > 0.05 ){
                
            // Filling the cuts histogram
            h_smeared->Fill(cos_mu_prime[i], T_mu_prime[i]);

            // Fill the nTuple for categorisation
            info->Fill( i, cos_mu_prime[i], T_mu_prime[i], neu, cc, nc, nipip, nipim, nipi0, qel, res, mec);
        }
        else if ( Impurity[i] == 1 
            && T_pi_vect[i] > 0.05 ){

            h_smeared->Fill(cos_pi_vect[i], T_pi_vect[i]);
            info->Fill( i, cos_pi_vect[i], T_pi_vect[i], neu, cc, nc, nipip, nipim, nipi0, qel, res, mec);
        }


    }

    c2->SetRightMargin(0.13);
   
    h_smeared->Draw("colz");
    h_smeared->SetStats(kFALSE);
    h_smeared->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_smeared->GetYaxis()->SetTitle("T_{#mu}");
    c2->SetLogz();
    //c2->SaveAs("distribution_after_impure.png");
    c2->SaveAs("smeared_cc0pi_2D.png");

    delete c2;
}

// ---------------------------------------------------------------------------
//                          Characterisation function
// ---------------------------------------------------------------------------

void Characterisation ( TH2D *h_smeared, TNtuple *info ){

    int n_entries  = info->GetEntries();

    // Get the  signal leaves
    TLeaf *l_ev   = info->GetLeaf( "ev" );
    TLeaf *l_cth  = info->GetLeaf( "cth" );
    TLeaf *l_T    = info->GetLeaf( "T" );
    TLeaf *l_neu  = info->GetLeaf( "neu" );
    TLeaf *l_cc   = info->GetLeaf( "cc" );
    TLeaf *l_nc   = info->GetLeaf( "nc" );
    TLeaf *l_qel  = info->GetLeaf( "qel" );
    TLeaf *l_res  = info->GetLeaf( "res" );
    TLeaf *l_mec  = info->GetLeaf( "mec" );
    TLeaf *l_npip = info->GetLeaf( "npip" );
    TLeaf *l_npim = info->GetLeaf( "npim" );
    TLeaf *l_npi0 = info->GetLeaf( "npi0" );

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_smeared->GetNbinsX(); // Cos theta
    int y_bins = h_smeared->GetNbinsY(); // Tmu

    vector< TH1D* > signal_h;
    
    TCanvas *c_Tmu         = new TCanvas ( "c_Tmu", "", 800, 600 );
    TLegend *leg_T         = new TLegend( 0.12, 0.70, 0.28, 0.88 );

    TH1D *h_Tmu_muccqe     = new TH1D ( "h_Tmu_muccqe",     "", x_bins, -1, 1 );
    TH1D *h_Tmu_muccmec    = new TH1D ( "h_Tmu_muccmec",    "", x_bins, -1, 1 );
    TH1D *h_Tmu_muccres1pi = new TH1D ( "h_Tmu_muccres1pi", "", x_bins, -1, 1 );
    TH1D *h_Tmu_mucc1pi    = new TH1D ( "h_Tmu_mucc1pi",    "", x_bins, -1, 1 );
    TH1D *h_Tmu_muccother  = new TH1D ( "h_Tmu_muccother",  "", x_bins, -1, 1 );
    TH1D *h_Tmu_mubarcc    = new TH1D ( "h_Tmu_mubarcc",    "", x_bins, -1, 1 );
    TH1D *h_Tmu_mubarnc    = new TH1D ( "h_Tmu_mubarnc",    "", x_bins, -1, 1 );
    
    signal_h.push_back(h_Tmu_muccqe);
    signal_h.push_back(h_Tmu_muccmec);
    signal_h.push_back(h_Tmu_muccres1pi);
    signal_h.push_back(h_Tmu_mucc1pi);
    signal_h.push_back(h_Tmu_muccother);
    signal_h.push_back(h_Tmu_mubarcc);
    signal_h.push_back(h_Tmu_mubarnc);

    // leg_T->AddEntry( h_Tmu_sm, " Total signal ", "f" );
    leg_T->AddEntry( h_Tmu_muccqe,     " #nu_{#mu} CCQE ", "f" );
    leg_T->AddEntry( h_Tmu_muccmec,    " #nu_{#mu} CCMEC ", "f" );
    leg_T->AddEntry( h_Tmu_muccres1pi, " #nu_{#mu} CCRES ", "f" );
    leg_T->AddEntry( h_Tmu_mucc1pi,    " #nu_{#mu} CC non-RES 1#pi ", "f" );
    leg_T->AddEntry( h_Tmu_muccother,  " #nu_{#mu} CCOther ", "f" );
    leg_T->AddEntry( h_Tmu_mubarcc,    " #bar{ #nu_{#mu} } CC ", "f" );
    leg_T->AddEntry( h_Tmu_mubarnc,    " #bar{ #nu_{#mu} } NC ", "f" );
   
    for ( int i = 1; i <= y_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_T;
        double up_edge_T;

        low_edge_T = h_smeared->GetYaxis()->GetBinLowEdge(i);
        up_edge_T  = h_smeared->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        // Filenames for signal histograms
        stringstream conv;
        conv.clear();

        string title;
        title.clear();

        char file_name[1024];

        conv << setprecision(4) << "Pre_FSI_Tmu_slice_" << low_edge_T << "_" << up_edge_T << ".png";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        stringstream hist;
        hist.clear();

        char hist_name[1024];
        string temp;

        hist << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "_" << up_edge_T << ", pre-FSI ";
        temp = hist.str();

        strcpy( hist_name, temp.c_str() );

        for( unsigned int i = 0; i < signal_h.size(); ++i ){
            signal_h[i]->Reset("M");
        }  

        THStack *hsT         = new THStack("hsT","pre-FSI ");
        
        for ( int k = 0; k < n_entries; ++k ){
       
            info->GetEntry(k);

            double cth = l_cth->GetValue();
            double T   = l_T->GetValue();
            int neu    = l_neu->GetValue();
            int ev     = l_ev->GetValue();
            int cc     = l_cc->GetValue();
            int nc     = l_nc->GetValue();
            int qel    = l_qel->GetValue();
            int res    = l_res->GetValue();
            int mec    = l_mec->GetValue();
            int npip   = l_npip->GetValue();
            int npim   = l_npim->GetValue();
            int npi0   = l_npi0->GetValue();

            // Muon neutrino in this bin
            if ( T >= low_edge_T 
              && T < up_edge_T
              && neu == 14
              && cc ){
                
                // CCQE
                if ( qel ) {
                    h_Tmu_muccqe->Fill( cth,1 );
                }
                
                // CCMEC
                else if ( mec ){
                    h_Tmu_muccmec->Fill( cth, 1 );
                }
                
                // CCRES, 1Pi
                else if ( res 
                     && ( npip + npim == 1 ) ){
                    h_Tmu_muccres1pi->Fill( cth, 1 );
                }
                // CC non-RES, 1Pi
                else if ( !res
                     && ( npip + npim == 1 ) ){
                    h_Tmu_mucc1pi->Fill( cth, 1 );
                }
                
                // CCOther
                else{
                    h_Tmu_muccother->Fill( cth, 1 );
                }
            }    
            else if ( T >= low_edge_T
                   && T < up_edge_T
                   && neu == -14 ){
                
                // Numubar CC
                if ( cc ){
                    h_Tmu_mubarcc->Fill( cth, 1 );
                }

                // Numubar NC
                else if ( nc ){
                    h_Tmu_mubarnc->Fill( cth, 1 );
                }
            } 
        }

        int pal[12];
        pal[0]  = kRed;
        pal[1]  = kGreen + 2;
        pal[2]  = kOrange + 7;
        pal[3]  = kBlue;
        pal[4]  = kMagenta + 1;
        pal[5]  = kCyan + 2;
        pal[6]  = kYellow + 1;
        pal[7]  = kAzure + 2;
        pal[8]  = kViolet + 2;
        pal[9]  = kTeal - 1;
        pal[10] = kPink + 2;
        pal[11] = kSpring + 2;
        
        gStyle->SetPalette( 12, pal );
        gStyle->SetHatchesLineWidth( 1 );
        gStyle->SetHatchesSpacing( 0.5 );
        gStyle->SetTitleOffset(1.5, "Y");
        gStyle->SetOptStat( 0 );

        double norm_Tmu = 0;
       
        for ( unsigned int i = 0; i < signal_h.size(); ++i ){
            norm_Tmu += signal_h[i]->Integral();
        }

        for ( unsigned int i = 0; i < signal_h.size(); ++i ){

            signal_h[i]->SetFillColor(pal[i]);
            signal_h[i]->SetLineColor(pal[i]);

            signal_h[i]->SetLineWidth(1.5);
            signal_h[i]->Scale(1/norm_Tmu);

            hsT->Add(signal_h[i]);
            
        }
        
        // Fill the stacks 
        // Signal
        hsT->Draw();
        
        hsT->SetTitle(hist_name);
        hsT->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        hsT->GetYaxis()->SetTitle("Number of SBND events");   
        
        leg_T->Draw();
        c_Tmu->SaveAs(file_name);

        delete hsT;
    } 

    for ( unsigned int i = 0; i < signal_h.size(); ++i ){
        delete signal_h[i];
    }
    delete c_Tmu;
    
    delete leg_T;
   
    // ------------------------------------------------------------------------------
    //                              Cos slices
    // ------------------------------------------------------------------------------
    
    vector< TH1D* > signal_c_h;
    
    TCanvas *c_cosmu = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_cos = new TLegend( 0.72, 0.70, 0.88, 0.88 );
 
    TH1D *h_cosmu_muccqe     = new TH1D ( "h_cosmu_muccqe",     "", y_bins, 0, 2 );
    TH1D *h_cosmu_muccmec    = new TH1D ( "h_cosmu_muccmec",    "", y_bins, 0, 2 );
    TH1D *h_cosmu_muccres1pi = new TH1D ( "h_cosmu_muccres1pi", "", y_bins, 0, 2 );
    TH1D *h_cosmu_mucc1pi    = new TH1D ( "h_cosmu_mucc1pi",    "", y_bins, 0, 2 );
    TH1D *h_cosmu_muccother  = new TH1D ( "h_cosmu_muccother",  "", y_bins, 0, 2 );
    TH1D *h_cosmu_mubarcc    = new TH1D ( "h_cosmu_mubarcc",    "", y_bins, 0, 2 );
    TH1D *h_cosmu_mubarnc    = new TH1D ( "h_cosmu_mubarnc",    "", y_bins, 0, 2 );
    
    signal_c_h.push_back(h_cosmu_muccqe);
    signal_c_h.push_back(h_cosmu_muccmec);
    signal_c_h.push_back(h_cosmu_muccres1pi);
    signal_c_h.push_back(h_cosmu_mucc1pi);
    signal_c_h.push_back(h_cosmu_muccother);
    signal_c_h.push_back(h_cosmu_mubarcc);
    signal_c_h.push_back(h_cosmu_mubarnc);

    // leg_cos->AddEntry( h_cosmu_sm, " Total signal ", "f" );
    leg_cos->AddEntry( h_cosmu_muccqe,     " #nu_{#mu} CCQE ", "f" );
    leg_cos->AddEntry( h_cosmu_muccmec,    " #nu_{#mu} CCMEC ", "f" );
    leg_cos->AddEntry( h_cosmu_muccres1pi, " #nu_{#mu} CCRES ", "f" );
    leg_cos->AddEntry( h_cosmu_mucc1pi,    " #nu_{#mu} CC non-RES 1#pi ", "f" );
    leg_cos->AddEntry( h_cosmu_muccother,  " #nu_{#mu} CCOther ", "f" );
    leg_cos->AddEntry( h_cosmu_mubarcc,    " #bar{ #nu_{#mu} } CC ", "f" );
    leg_cos->AddEntry( h_cosmu_mubarnc,    " #bar{ #nu_{#mu} } NC ", "f" );
   
    for ( int i = 1; i <= x_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_cos;
        double up_edge_cos;

        low_edge_cos = h_smeared->GetXaxis()->GetBinLowEdge(i);
        up_edge_cos  = h_smeared->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        // Filenames for signal histograms
        stringstream conv_c;
        conv_c.clear();

        string title_c;
        title_c.clear();

        char file_name_c[1024];

        conv_c << setprecision(4) << "pre_FSI_cosmu_slice_" << low_edge_cos << "_" << up_edge_cos << ".png";
        title_c = conv_c.str();
        
        strcpy( file_name_c, title_c.c_str() );

        // For the title of the histogram
        stringstream hist_c;
        hist_c.clear();

        char hist_name_c[1024];
        string temp_c;

        hist_c << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "_" << up_edge_cos << ", pre-FSI ";
        temp_c = hist_c.str();

        strcpy( hist_name_c, temp_c.c_str() );


        for( unsigned int i = 0; i < signal_c_h.size(); ++i ){
            signal_c_h[i]->Reset("M");
        }  

        THStack *hscos = new THStack("hscos","pre_FSI");
            
        for ( int k = 0; k < n_entries; ++k ){
       
            info->GetEntry(k);

            int ev     = l_ev->GetValue();
            double cth = l_cth->GetValue();
            double T   = l_T->GetValue();
            int neu    = l_neu->GetValue();
            int cc     = l_cc->GetValue();
            int nc     = l_nc->GetValue();
            int qel    = l_qel->GetValue();
            int res    = l_res->GetValue();
            int mec    = l_mec->GetValue();
            int npip   = l_npip->GetValue();
            int npim   = l_npim->GetValue();
            int npi0   = l_npi0->GetValue();

            // Muon neutrino in this bin
            if ( cth >= low_edge_cos 
              && cth < up_edge_cos
              && neu == 14
              && cc ){
                
                // CCQE
                if ( qel ) {
                    h_cosmu_muccqe->Fill( T, 1 );
                }
                
                // CCMEC
                else if ( mec ){
                    h_cosmu_muccmec->Fill( T, 1 );
                }
                
                // CCRES, 1Pi
                else if ( res 
                     && ( npip + npim == 1 ) ){
                    h_cosmu_muccres1pi->Fill( T, 1 );
                }
                // CC non-RES, 1Pi
                else if ( !res
                     && ( npip + npim == 1 ) ){
                    h_cosmu_mucc1pi->Fill( T, 1 );
                }
                
                // CCOther
                else{
                    h_cosmu_muccother->Fill( T, 1 );
                }
            }    
            else if ( cth >= low_edge_cos
                   && cth < up_edge_cos
                   && neu == -14 ){
                
                // Numubar CC
                if ( cc ){
                    h_cosmu_mubarcc->Fill( T, 1 );
                }

                // Numubar NC
                else if ( nc ){
                    h_cosmu_mubarnc->Fill( T, 1 );
                }
            } 
        }
        
        int pal[12];
        pal[0]  = kRed;
        pal[1]  = kGreen + 2;
        pal[2]  = kOrange + 7;
        pal[3]  = kBlue;
        pal[4]  = kMagenta + 1;
        pal[5]  = kCyan + 2;
        pal[6]  = kYellow + 1;
        pal[7]  = kAzure + 2;
        pal[8]  = kViolet + 2;
        pal[9]  = kTeal - 1;
        pal[10] = kPink + 2;
        pal[11] = kSpring + 2;
        
        gStyle->SetPalette( 12, pal );
        gStyle->SetHatchesLineWidth( 1 );
        gStyle->SetHatchesSpacing( 0.5 );
        gStyle->SetTitleOffset(1.5, "Y");
        gStyle->SetOptStat(0);

        double norm_cosmu = 0;
       
        for ( unsigned int i = 0; i < signal_c_h.size(); ++i ){
            norm_cosmu += signal_c_h[i]->Integral();
        }
        
        for ( unsigned int i = 0; i < signal_c_h.size(); ++i ){

            signal_c_h[i]->SetFillColor(pal[i]);
            signal_c_h[i]->SetLineColor(pal[i]);

            signal_c_h[i]->SetLineWidth(1.5);
            signal_c_h[i]->Scale(1/norm_cosmu);

            hscos->Add(signal_c_h[i]);
        }
        
        // Fill the stacks 
        // Signal
        hscos->Draw();
        
        hscos->SetTitle(hist_name_c);
        hscos->GetXaxis()->SetTitle("T_{#mu}");   
        hscos->GetYaxis()->SetTitle("Number of SBND events");   
    
        leg_cos->Draw();
        c_cosmu->SaveAs(file_name_c);

        delete hscos;
    } 

    for ( unsigned int i = 0; i < signal_c_h.size(); ++i ){
        delete signal_c_h[i];
    }
    
    delete c_cosmu;
    
    delete leg_cos;

}


