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

#include "post_fsi.h"

using namespace std; 

int less_impure_post_fsi() {

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
        cout << " Default+MEC event file is open " << endl;
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
   

    gStyle->SetPalette(55);

    h_un->SetStats(kFALSE);
    h_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un->GetYaxis()->SetTitle("T_{#mu}");
    c1->SetRightMargin(0.13);
    c1->SetLogz();
    //c1->SaveAs("distribution_before.png");
    c1->SaveAs("distribution_before_log.png");
    
    // The same histogram definitions for the smeared distributions
    //TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, impurity ",20,-1,1,18,0,2);
    TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} event rates, after smearing ",20,-1,1,18,0,2);

    // Histogram to hold only the cut values
    TH2D *h_c = new TH2D("h_c"," Cut applied ",20,-1,1,18,0,2);
    def_tree->Draw("( El - 0.10566 ):cthl>>h_c","( El - 0.10566 ) > 0.05 && fspl == 13 && cc && (nfpip + nfpim + nfpi0 == 0)","colz"); 
    
    // Initiate the TNtuple to hold the interesting features
    TNtuple *nt_s = new TNtuple("nt_s","Important information about signal","ev:cth:T:cc:nc:np:npip:npim:npi0:nkp:nkm:nk0:qel:res:mec");
    nt_s->SetDirectory(0);

    TNtuple *nt_bg = new TNtuple("nt_bg","Important information about background","ev:cth:T:cc:nc:np:npip:npim:npi0:nkp:nkm:nk0:qel:res:mec");
    nt_bg->SetDirectory(0);
    
    // Take h_un and smear it
    Smear(def_tree, nt_s, nt_bg, h_un, h_sm, h_c);

    // The smeared histogram and an ntuple containing the important information
    // to help characterise the signal
    Characterisation( h_sm, nt_s, nt_bg );
 
    // The efficiency function 
    Efficiency( h_sm, h_un, h_c, nt_s );

    delete h_sm;
    delete h_un;
    delete c1;
    delete nt_s;
    delete nt_bg;

    return 0;
}

// -------------------------------------------------------------------------
//                            Function definitions
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
//                             Smearing function
// -------------------------------------------------------------------------

void Smear ( TTree   *tree, 
             TNtuple *info_s,
             TNtuple *info_bg,
             TH2D  *h_unsmeared,
             TH2D  *h_smeared,
             TH2D  *h_cut ){

    ofstream file1;
    file1.open("unsmeared_results.tex"); 
    file1 << " \\begin{tabular}{ | * {20}{c} | } " << endl;
    
    ofstream file2;
    file2.open("impure_smeared_results.tex"); 
    file2 << " \\begin{tabular}{ | * {20}{c} | } " << endl;
     
    ofstream file3;
    file3.open("impure_results_difference.tex"); 
    file3 << " \\begin{tabular}{ | * {20}{c} | } " << endl;

    // Get the branches
    TBranch *b_nf    = tree->GetBranch( "nf" );
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
    TBranch *b_nfp   = tree->GetBranch( "nfp" );
    TBranch *b_nfpip = tree->GetBranch( "nfpip" );
    TBranch *b_nfpim = tree->GetBranch( "nfpim" );
    TBranch *b_nfpi0 = tree->GetBranch( "nfpi0" );
    TBranch *b_nfkp  = tree->GetBranch( "nfkp" );
    TBranch *b_nfkm  = tree->GetBranch( "nfkm" );
    TBranch *b_nfk0  = tree->GetBranch( "nfk0" );

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
  
    vector< int > Impurity; // Smearing: 2 % 

    vector< int > CC1cpi_Imp;  // BG : 2 %
    vector< int > CC1pi0_Imp;  // BG : 1 %
    vector< int > CCOther_Imp; // BG : 1 %

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
        // Signal impurity
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
       
        // Apply the impurity cut
        // CC1cpi 
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 && nfpip + nfpim == 1 && nfpi0 == 0 ){
            int random1;
            random1 = rand() % 10 + 1;

            // If the random number = 5 then push back a 1 onto the vector
            // If the random number < 5 push back a 0
            if( random1 == 10 ){
                CC1cpi_Imp.push_back(1);       
            }
            else{
                CC1cpi_Imp.push_back(0);
            }
        }
        else{
            CC1cpi_Imp.push_back(0);
        }
        
        // Apply the impurity cut
        // CC1pi0 
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 && nfpip + nfpim == 0 && nfpi0 == 1 ){
            int random2;
            random2 = rand() % 10 + 1;

            // If the random number = 5 then push back a 1 onto the vector
            // If the random number < 5 push back a 0
            if( random2 == 10 ){
                CC1pi0_Imp.push_back(1);       
            }
            else{
                CC1pi0_Imp.push_back(0);
            }
        }
        else{
            CC1pi0_Imp.push_back(0);
        }
        
        // Apply the impurity cut
        // CCOther 
        if ( b_cc->GetLeaf("cc")->GetValue() != 0 && nfpip + nfpim + nfpi0 > 1 ){
            int random3;
            random3 = rand() % 10 + 1;

            // If the random number = 5 then push back a 1 onto the vector
            // If the random number < 5 push back a 0
            if( random3 == 10 ){
                CCOther_Imp.push_back(1);       
            }
            else{
                CCOther_Imp.push_back(0);
            }
        }
        else{
            CCOther_Imp.push_back(0);
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

    // Fill h_smeared normally and with the impurities and THEN clone it
    for ( int i = 0; i < n_values; ++i ){        
        
        tree->GetEntry(i); 
       
        int cc    = b_cc->GetLeaf( "cc" )->GetValue(); 
        int nc    = b_nc->GetLeaf( "nc" )->GetValue(); 
        int qel   = b_qel->GetLeaf( "qel" )->GetValue();
        int res   = b_res->GetLeaf( "res" )->GetValue();
        int mec   = b_mec->GetLeaf( "mec" )->GetValue();
        int fspl  = b_fspl->GetLeaf( "fspl" )->GetValue(); 
        int nfp   = b_nfp->GetLeaf( "nfp" )->GetValue();
        int nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        int nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        int nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();
        int nfkp  = b_nfkp->GetLeaf( "nfkp" )->GetValue();
        int nfkm  = b_nfkm->GetLeaf( "nfkm" )->GetValue();
        int nfk0  = b_nfk0->GetLeaf( "nfk0" )->GetValue();

        // Set the energy and impurity cuts
        // Filling the signal ntuple
        if ( fspl == 13 
          && cc
          && ( nfpip + nfpim + nfpi0 == 0 )
          && T_mu_vect[i] > 0.05 ){
                
            // Filling the cuts histogram
            h_smeared->Fill(cos_mu_prime[i], T_mu_prime[i]);
            info_s->Fill( i, cos_mu_prime[i], T_mu_prime[i], cc, nc, nfp, nfpip, nfpim, nfpi0, nfkp, nfkm, nfk0, qel, res, mec );
        }
    
        // Filling the background ntuple
        else if ( Impurity[i] == 1 
               && T_pi_vect[i] > 0.05 ){

            h_smeared->Fill(cos_pi_vect[i], T_pi_vect[i]);
            info_bg->Fill( i, cos_pi_vect[i], T_pi_vect[i], cc, nc, nfp, nfpip, nfpim, nfpi0, nfkp, nfkm, nfk0, qel, res, mec );
        }
        else if ( fspl == 13
            && T_mu_prime[i] > 0.05
            && T_pi_vect[i] > 0.05
            && CC1cpi_Imp[i] == 1 ){
         
            info_bg->Fill( i, cos_mu_prime[i], T_mu_prime[i], cc, nc, nfp, nfpip, nfpim, nfpi0, nfkp, nfkm, nfk0, qel, res, mec );
        }
        else if ( fspl == 13
            && T_mu_prime[i] > 0.05
            && T_pi_vect[i] > 0.05
            && CC1pi0_Imp[i] == 1 ){
         
            info_bg->Fill( i, cos_mu_prime[i], T_mu_prime[i], cc, nc, nfp, nfpip, nfpim, nfpi0, nfkp, nfkm, nfk0, qel, res, mec );
        }
        else if ( fspl == 13
            && T_mu_vect[i] > 0.05
            && CCOther_Imp[i] == 1 ){
         
            info_bg->Fill( i, cos_mu_prime[i], T_mu_prime[i], cc, nc, nfp, nfpip, nfpim, nfpi0, nfkp, nfkm, nfk0, qel, res, mec );
        }
    }

    int nX = h_smeared->GetNbinsX();
    int nY = h_smeared->GetNbinsY();
 
    // Write the bin contents to .tex files
    for ( int i = nY; i >= 1; --i){
        for ( int j = 1; j <= nX - 1; ++j ){
            
            file1 << h_unsmeared->GetBinContent(j,i) << " & ";
            file2 << h_smeared->GetBinContent(j,i) << " & "; 
            file3 << ( h_smeared->GetBinContent(j,i) - h_unsmeared->GetBinContent(j,i) ) << " & ";
        }

        file1 << h_unsmeared->GetBinContent(nX, i) << " \\\\ " << endl;
        file2 << h_smeared->GetBinContent(nX, i) << " \\\\ " << endl;
        file3 << ( h_smeared->GetBinContent(nX, i) - h_unsmeared->GetBinContent(nX, i) ) << " \\\\ " << endl;
        
    }
    file1 << " \\end{tabular} " << endl;
    file2 << " \\end{tabular} " << endl;
    file3 << " \\end{tabular} " << endl;

    c2->SetRightMargin(0.13);
   
    h_smeared->Draw("colz");
    h_smeared->SetStats(kFALSE);
    h_smeared->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_smeared->GetYaxis()->SetTitle("T_{#mu}");
    c2->SetLogz();
    //c2->SaveAs("distribution_after_impure.png");
    c2->SaveAs("distribution_after_log_impure.png");

    delete c2;
}

// ---------------------------------------------------------------------------
//                          Characterisation function
// ---------------------------------------------------------------------------

void Characterisation ( TH2D *h_smeared, TNtuple *info_s, TNtuple *info_bg ){

    int n_entries  = info_s->GetEntries();
    int bg_entries = info_bg->GetEntries();

    // Get the  signal leaves
    TLeaf *l_s_ev   = info_s->GetLeaf( "ev" );
    TLeaf *l_s_cth  = info_s->GetLeaf( "cth" );
    TLeaf *l_s_T    = info_s->GetLeaf( "T" );
    TLeaf *l_s_cc   = info_s->GetLeaf( "cc" );
    TLeaf *l_s_nc   = info_s->GetLeaf( "nc" );
    TLeaf *l_s_qel  = info_s->GetLeaf( "qel" );
    TLeaf *l_s_res  = info_s->GetLeaf( "res" );
    TLeaf *l_s_mec  = info_s->GetLeaf( "mec" );
    TLeaf *l_s_np   = info_s->GetLeaf( "np" );
    TLeaf *l_s_npip = info_s->GetLeaf( "npip" );
    TLeaf *l_s_npim = info_s->GetLeaf( "npim" );
    TLeaf *l_s_npi0 = info_s->GetLeaf( "npi0" );
    TLeaf *l_s_nkp  = info_s->GetLeaf( "nkp" );
    TLeaf *l_s_nkm  = info_s->GetLeaf( "nkm" );
    TLeaf *l_s_nk0  = info_s->GetLeaf( "nk0" );

    // Get the  background leaves
    TLeaf *l_bg_ev   = info_bg->GetLeaf( "ev" );
    TLeaf *l_bg_cth  = info_bg->GetLeaf( "cth" );
    TLeaf *l_bg_T    = info_bg->GetLeaf( "T" );
    TLeaf *l_bg_cc   = info_bg->GetLeaf( "cc" );
    TLeaf *l_bg_nc   = info_bg->GetLeaf( "nc" );
    TLeaf *l_bg_qel  = info_bg->GetLeaf( "qel" );
    TLeaf *l_bg_res  = info_bg->GetLeaf( "res" );
    TLeaf *l_bg_mec  = info_bg->GetLeaf( "mec" );
    TLeaf *l_bg_np   = info_bg->GetLeaf( "np" );
    TLeaf *l_bg_npip = info_bg->GetLeaf( "npip" );
    TLeaf *l_bg_npim = info_bg->GetLeaf( "npim" );
    TLeaf *l_bg_npi0 = info_bg->GetLeaf( "npi0" );
    TLeaf *l_bg_nkp  = info_bg->GetLeaf( "nkp" );
    TLeaf *l_bg_nkm  = info_bg->GetLeaf( "nkm" );
    TLeaf *l_bg_nk0  = info_bg->GetLeaf( "nk0" );

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_smeared->GetNbinsX(); // Cos theta
    int y_bins = h_smeared->GetNbinsY(); // Tmu

    TCanvas *c_Tmu           = new TCanvas ( "c_Tmu", "", 800, 600 );
    
    TLegend *leg_T_s         = new TLegend( 0.12, 0.70, 0.28, 0.88);
    TLegend *leg_T_bg        = new TLegend( 0.12, 0.74, 0.28, 0.88 );
    TLegend *leg_T_sbg       = new TLegend( 0.12, 0.74, 0.28, 0.88 );
    TLegend *leg_T_sbg_split = new TLegend( 0.12, 0.60, 0.28, 0.88 );
    
    TH1D *h_Tmu_sm        = new TH1D ( "h_Tmu_sm", "", x_bins, -1, 1 );
    TH1D *h_Tmu_bg        = new TH1D ( "h_Tmu_bg", "", x_bins, -1, 1 );
 
    vector< TH1D* > signal_h;
    vector< TH1D* > background_h;
    
    TH1D *h_Tmu_0p_s     = new TH1D ( "h_Tmu_0p_s", "", x_bins, -1, 1 );
    TH1D *h_Tmu_1p_s     = new TH1D ( "h_Tmu_1p_s", "", x_bins, -1, 1 );
    TH1D *h_Tmu_2p_s     = new TH1D ( "h_Tmu_2p_s", "", x_bins, -1, 1 );
    TH1D *h_Tmu_3p_s     = new TH1D ( "h_Tmu_3p_s", "", x_bins, -1, 1 );
    TH1D *h_Tmu_mp_s     = new TH1D ( "h_Tmu_mp_s", "", x_bins, -1, 1 );
    
    TH1D *h_Tmu_1pip_bg   = new TH1D ( "h_Tmu_1pip_bg", "", x_bins, -1, 1 );
    TH1D *h_Tmu_1pim_bg   = new TH1D ( "h_Tmu_1pim_bg", "", x_bins, -1, 1 );
    TH1D *h_Tmu_1pi0_bg   = new TH1D ( "h_Tmu_1pi0_bg", "", x_bins, -1, 1 );
    TH1D *h_Tmu_2cpi_bg   = new TH1D ( "h_Tmu_2cpi_bg", "", x_bins, -1, 1 );
    TH1D *h_Tmu_n1cpi_bg  = new TH1D ( "h_Tmu_n1cpi_bg", "", x_bins, -1, 1 );
    TH1D *h_Tmu_other_bg  = new TH1D ( "h_Tmu_other_bg", "", x_bins, -1, 1 );
    TH1D *h_Tmu_nother_bg = new TH1D ( "h_Tmu_nother_bg", "", x_bins, -1, 1 );
    
    signal_h.push_back(h_Tmu_mp_s);
    signal_h.push_back(h_Tmu_3p_s);
    signal_h.push_back(h_Tmu_2p_s);
    signal_h.push_back(h_Tmu_1p_s);
    signal_h.push_back(h_Tmu_0p_s);

    background_h.push_back(h_Tmu_nother_bg);
    background_h.push_back(h_Tmu_n1cpi_bg);
    background_h.push_back(h_Tmu_other_bg);
    background_h.push_back(h_Tmu_2cpi_bg);
    background_h.push_back(h_Tmu_1pi0_bg);
    background_h.push_back(h_Tmu_1pim_bg);
    background_h.push_back(h_Tmu_1pip_bg);
    
    leg_T_s->AddEntry( h_Tmu_0p_s, " #nu_{#mu} CC0#pi, 0 p ", "f" );
    leg_T_s->AddEntry( h_Tmu_1p_s, " #nu_{#mu} CC0#pi, 1 p ", "f" );
    leg_T_s->AddEntry( h_Tmu_2p_s, " #nu_{#mu} CC0#pi, 2 p ", "f" );
    leg_T_s->AddEntry( h_Tmu_3p_s, " #nu_{#mu} CC0#pi, 3 p ", "f" );
    leg_T_s->AddEntry( h_Tmu_mp_s, " #nu_{#mu} CC0#pi, >3 p ", "f" );
   
    leg_T_bg->AddEntry( h_Tmu_1pip_bg, " #nu_{#mu} CC1#pi^{ + } ", "f" );
    leg_T_bg->AddEntry( h_Tmu_1pim_bg, " #nu_{#mu} CC1#pi^{ - } ", "f" );
    leg_T_bg->AddEntry( h_Tmu_1pi0_bg, " #nu_{#mu} CC1#pi^{ 0 } ", "f" );
    leg_T_bg->AddEntry( h_Tmu_2cpi_bg, " #nu_{#mu} CC2#pi^{ #pm } ", "f" );
    leg_T_bg->AddEntry( h_Tmu_other_bg, " #nu_{#mu} CCOther ", "f" );
    leg_T_bg->AddEntry( h_Tmu_n1cpi_bg, " #nu_{#mu} NC1#pi^{ #pm } ", "f" );
    leg_T_bg->AddEntry( h_Tmu_nother_bg, " #nu_{#mu} NCOther ", "f" );
    
    leg_T_sbg->AddEntry( h_Tmu_sm, " Signal ", "f" );
    leg_T_sbg->AddEntry( h_Tmu_bg, " Background ", "f" );
    
    leg_T_sbg_split->AddEntry( h_Tmu_0p_s, " #nu_{#mu} CC0#pi, 0 p ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_1p_s, " #nu_{#mu} CC0#pi, 1 p ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_2p_s, " #nu_{#mu} CC0#pi, 2 p ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_3p_s, " #nu_{#mu} CC0#pi, 3 p ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_mp_s, " #nu_{#mu} CC0#pi, >3 p ", "f" ); 
    leg_T_sbg_split->AddEntry( h_Tmu_1pip_bg, " #nu_{#mu} CC1#pi^{ + } ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_1pim_bg, " #nu_{#mu} CC1#pi^{ - } ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_1pi0_bg, " #nu_{#mu} CC1#pi^{ 0 } ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_2cpi_bg, " #nu_{#mu} CC2#pi^{ #pm } ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_other_bg, " #nu_{#mu} CCOther ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_n1cpi_bg, " #nu_{#mu} NC1#pi^{ #pm } ", "f" );
    leg_T_sbg_split->AddEntry( h_Tmu_nother_bg, " #nu_{#mu} NCOther ", "f" );
    
    leg_T_s->SetHeader("True final state");
    leg_T_bg->SetHeader("True final state");
    leg_T_sbg->SetHeader("True final state");
    leg_T_sbg_split->SetHeader("True final state");
    
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
        stringstream conv_s;
        conv_s.clear();

        string title_s;
        title_s.clear();

        char file_name_s[1024];

        conv_s << setprecision(4) << "Signal_Tmu_slice_" << low_edge_T << "_" << up_edge_T << ".png";
        title_s = conv_s.str();
        
        strcpy( file_name_s, title_s.c_str() );

        // Filenames for background histograms
        stringstream conv_bg;
        conv_bg.clear();

        string title_bg;
        title_bg.clear();

        char file_name_bg[1024];

        conv_bg << setprecision(4) << "BG_Tmu_slice_" << low_edge_T << "_" << up_edge_T << ".png";
        title_bg = conv_bg.str();
        
        strcpy( file_name_bg, title_bg.c_str() );
        
        // Filenames for signal and background histograms
        stringstream conv_sbg;
        conv_sbg.clear();

        string title_sbg;
        title_sbg.clear();

        char file_name_sbg[1024];

        conv_sbg << setprecision(4) << "SBG_Tmu_slice_" << low_edge_T << "_" << up_edge_T << ".png";
        title_sbg = conv_sbg.str();
        
        strcpy( file_name_sbg, title_sbg.c_str() );
        
        // Filenames for signal and background split histograms
        stringstream conv_sbg_split;
        conv_sbg_split.clear();

        string title_sbg_split;
        title_sbg_split.clear();

        char file_name_sbg_split[1024];

        conv_sbg_split << setprecision(4) << "Split_SBG_Tmu_slice_" << low_edge_T << "_" << up_edge_T << ".png";
        title_sbg_split = conv_sbg_split.str();
        
        strcpy( file_name_sbg_split, title_sbg_split.c_str() );
        
        // For the title of the histogram
        stringstream hist_s;
        hist_s.clear();

        char hist_name_s[1024];
        string temp_s;

        hist_s << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "_" << up_edge_T << ", Signal ";
        temp_s = hist_s.str();

        strcpy( hist_name_s, temp_s.c_str() );

        // For the title of the histogram
        stringstream hist_bg;
        hist_bg.clear();

        char hist_name_bg[1024];
        string temp_bg;

        hist_bg << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "_" << up_edge_T << ", Background ";
        temp_bg = hist_bg.str();

        strcpy( hist_name_bg, temp_bg.c_str() );

        // For the title of the histogram
        stringstream hist_sbg;
        hist_sbg.clear();

        char hist_name_sbg[1024];
        string temp_sbg;

        hist_sbg << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "_" << up_edge_T << ", Signal and Background";
        temp_sbg = hist_sbg.str();

        strcpy( hist_name_sbg, temp_sbg.c_str() );

        // For the title of the histogram
        stringstream hist_sbg_split;
        hist_sbg_split.clear();

        char hist_name_sbg_split[1024];
        string temp_sbg_split;

        hist_sbg_split << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "_" << up_edge_T << ", Signal and Background Split";
        temp_sbg_split = hist_sbg_split.str();

        strcpy( hist_name_sbg_split, temp_sbg_split.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= x_bins; ++j ){
            h_Tmu_sm->SetBinContent( j, h_smeared->GetBinContent(j, i) );
        }

        h_Tmu_bg->Reset("M");

        for( unsigned int i = 0; i < signal_h.size(); ++i ){
            signal_h[i]->Reset("M");
        }  

        for( unsigned int i = 0; i < background_h.size(); ++i ){
            background_h[i]->Reset("M");
        }  
        
        THStack *hsT_s         = new THStack("hsT_s","Signal histograms");
        THStack *hsT_bg        = new THStack("hsT_bg","Backgound histograms");
        THStack *hsT_sbg       = new THStack("hsT_sbg","Total Signal and Background histograms");
        THStack *hsT_sbg_split = new THStack("hsT_sbg_split","Signal and Background histograms");

        for ( int k = 0; k < n_entries; ++k ){
       
            info_s->GetEntry(k);

            double cth_s = l_s_cth->GetValue();
            double T_s   = l_s_T->GetValue();
            int ev_s     = l_s_ev->GetValue();
            int cc_s     = l_s_cc->GetValue();
            int nc_s     = l_s_nc->GetValue();
            int np_s     = l_s_np->GetValue();
            int npip_s   = l_s_npip->GetValue();
            int npim_s   = l_s_npim->GetValue();
            int npi0_s   = l_s_npi0->GetValue();

            if ( T_s >= low_edge_T 
                      && T_s < up_edge_T
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s == 0 ) {
             
                h_Tmu_0p_s->Fill(cth_s,1);
            }
            else if ( T_s >= low_edge_T 
                      && T_s < up_edge_T
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s == 1 ) {
             
                h_Tmu_1p_s->Fill(cth_s,1);
            }
            else if ( T_s >= low_edge_T 
                      && T_s < up_edge_T
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s == 2 ) {
             
                h_Tmu_2p_s->Fill(cth_s,1);
            }
            else if ( T_s >= low_edge_T 
                      && T_s < up_edge_T
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s == 3 ) {
             
                h_Tmu_3p_s->Fill(cth_s,1);
            }
            else if ( T_s >= low_edge_T 
                      && T_s < up_edge_T
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s > 3 ) {
             
                h_Tmu_mp_s->Fill(cth_s,1);
            }
        }

        for ( int k = 0; k < bg_entries; ++k ){
       
            info_bg->GetEntry(k);

            double cth_bg = l_bg_cth->GetValue();
            double T_bg   = l_bg_T->GetValue();
            int ev_bg     = l_bg_ev->GetValue();
            int cc_bg     = l_bg_cc->GetValue();
            int nc_bg     = l_bg_nc->GetValue();
            int np_bg     = l_bg_np->GetValue();
            int npip_bg   = l_bg_npip->GetValue();
            int npim_bg   = l_bg_npim->GetValue();
            int npi0_bg   = l_bg_npi0->GetValue();
            
            if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && cc_bg
                      && npip_bg ==  1 
                      && ( npi0_bg + npim_bg == 0 ) ) {
                
                h_Tmu_1pip_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && cc_bg
                      && npim_bg ==  1 
                      && ( npi0_bg + npip_bg == 0 ) ) {
                
                h_Tmu_1pim_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && cc_bg
                      && npi0_bg ==  1 
                      && ( npim_bg + npip_bg == 0 ) ) {
                
                h_Tmu_1pi0_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && cc_bg
                      && ( npip_bg + npim_bg ==  2  )
                      && npi0_bg == 0 ) {
                
                h_Tmu_2cpi_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && nc_bg
                      && ( npip_bg + npim_bg ==  1  )
                      && npi0_bg == 0 ) {
                
                h_Tmu_n1cpi_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && cc_bg ) {
                 
                h_Tmu_other_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && nc_bg ) {
                 
                h_Tmu_nother_bg->Fill(cth_bg,1);
            }
        }

        for ( int k = 0; k < bg_entries; ++k ){
       
            info_bg->GetEntry(k);

            double cth_bg = l_bg_cth->GetValue();
            double T_bg   = l_bg_T->GetValue();
            int ev_bg     = l_bg_ev->GetValue();
            int cc_bg     = l_bg_cc->GetValue();
            int nc_bg     = l_bg_nc->GetValue();
            int np_bg   = l_bg_npip->GetValue();
            int npip_bg   = l_bg_npip->GetValue();
            int npim_bg   = l_bg_npim->GetValue();
            int npi0_bg   = l_bg_npi0->GetValue();
            
            if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && cc_bg
                      && ( npip_bg + npim_bg ==  1  )
                      && npi0_bg == 0 ) {
                
                h_Tmu_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && cc_bg
                      && ( npip_bg + npim_bg ==  2  )
                      && npi0_bg == 0 ) {
                
                h_Tmu_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && nc_bg
                      && ( npip_bg + npim_bg ==  1  )
                      && npi0_bg == 0 ) {
                
                h_Tmu_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && cc_bg ) {
                 
                h_Tmu_bg->Fill(cth_bg,1);
            }
            else if ( T_bg >= low_edge_T 
                      && T_bg < up_edge_T
                      && nc_bg ) {
                 
                h_Tmu_bg->Fill(cth_bg,1);
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

        h_Tmu_sm->SetFillColor( kGreen + 1 );
        h_Tmu_bg->SetFillColor( kRed + 1 );
        
        h_Tmu_bg->SetFillStyle( 3354 );
        
        double norm_Tmu = h_Tmu_bg->Integral() + h_Tmu_sm->Integral();
        
        h_Tmu_sm->Scale(1/norm_Tmu);
        h_Tmu_bg->Scale(1/norm_Tmu);
        
        for ( unsigned int i = 0; i < background_h.size(); ++i ){

            int s = signal_h.size();

            background_h[i]->SetFillColor(pal[s + i]);
            background_h[i]->SetLineColor(pal[s + i]);

            background_h[i]->SetFillStyle( 3354 );
            background_h[i]->SetLineWidth(1.5);
            background_h[i]->Scale(1/norm_Tmu);

            hsT_bg->Add(background_h[i]);
            hsT_sbg_split->Add(background_h[i]);
            
        }
        for ( unsigned int i = 0; i < signal_h.size(); ++i ){

            signal_h[i]->SetFillColor(pal[i]);
            signal_h[i]->SetLineColor(pal[i]);

            signal_h[i]->SetLineWidth(1.5);
            signal_h[i]->Scale(1/norm_Tmu);

            hsT_s->Add(signal_h[i]);
            hsT_sbg_split->Add(signal_h[i]);
            
        }
        
        h_Tmu_sm->SetLineColor( kGreen + 1 );
        h_Tmu_bg->SetLineColor( kRed + 1 );

        h_Tmu_sm->SetLineWidth(1.5);     
        h_Tmu_bg->SetLineWidth(1.5);     
        
        // Fill the stacks 
        // Signal
        hsT_s->Draw();
        
        hsT_s->SetTitle(hist_name_s);
        hsT_s->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        hsT_s->GetYaxis()->SetTitle("Number of SBND events");   
        
        leg_T_s->Draw();
        c_Tmu->SaveAs(file_name_s);

        // Background 
        hsT_bg->Draw();
        
        hsT_bg->SetTitle(hist_name_bg);
        hsT_bg->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        hsT_bg->GetYaxis()->SetTitle("Number of SBND events");   
    
        leg_T_bg->Draw();
        c_Tmu->SaveAs(file_name_bg);

        // Signal and background
        hsT_sbg->Add(h_Tmu_bg);
        hsT_sbg->Add(h_Tmu_sm);
        
        hsT_sbg->Draw();
        
        hsT_sbg->SetTitle(hist_name_sbg);
        hsT_sbg->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        hsT_sbg->GetYaxis()->SetTitle("Number of SBND events");   
        
        leg_T_sbg->Draw();
        c_Tmu->SaveAs(file_name_sbg);

        // Signal and background split
        hsT_sbg_split->Draw();
        
        hsT_sbg_split->SetTitle(hist_name_sbg_split);
        hsT_sbg_split->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        hsT_sbg_split->GetYaxis()->SetTitle("Number of SBND events");   
        
        leg_T_sbg_split->Draw();
        c_Tmu->SaveAs(file_name_sbg_split);
        
        delete hsT_s;
        delete hsT_sbg;
        delete hsT_bg;
        delete hsT_sbg_split;
    } 

    delete h_Tmu_sm;
    delete h_Tmu_bg;

    for ( unsigned int i = 0; i < signal_h.size(); ++i ){
        delete signal_h[i];
    }
    
    for ( unsigned int i = 0; i < background_h.size(); ++i ){
        delete background_h[i];
    }
    
    delete c_Tmu;
    
    delete leg_T_s;
    delete leg_T_bg;
    delete leg_T_sbg;
    delete leg_T_sbg_split;
   
    // ------------------------------------------------------------------------------
    //                              Cos slices
    // ------------------------------------------------------------------------------
    
    TCanvas *c_cosmu = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_cos_s         = new TLegend( 0.72, 0.70, 0.88, 0.88 );
    TLegend *leg_cos_bg        = new TLegend( 0.72, 0.74, 0.88, 0.88 );
    TLegend *leg_cos_sbg       = new TLegend( 0.72, 0.74, 0.88, 0.88 );
    TLegend *leg_cos_sbg_split = new TLegend( 0.72, 0.60, 0.88, 0.88 ); 

    TH1D *h_cosmu_sm        = new TH1D ( "h_cosmu_sm", "", y_bins, 0, 2 );
    TH1D *h_cosmu_bg        = new TH1D ( "h_cosmu_bg", "", y_bins, 0, 2 );
 
    vector< TH1D* > signal_c_h;
    vector< TH1D* > background_c_h;
    
    TH1D *h_cosmu_0p_s     = new TH1D ( "h_cosmu_0p_s", "", y_bins, 0, 2 );
    TH1D *h_cosmu_1p_s     = new TH1D ( "h_cosmu_1p_s", "", y_bins, 0, 2 );
    TH1D *h_cosmu_2p_s     = new TH1D ( "h_cosmu_2p_s", "", y_bins, 0, 2 );
    TH1D *h_cosmu_3p_s     = new TH1D ( "h_cosmu_3p_s", "", y_bins, 0, 2 );
    TH1D *h_cosmu_mp_s     = new TH1D ( "h_cosmu_mp_s", "", y_bins, 0, 2 );
    
    TH1D *h_cosmu_1pip_bg   = new TH1D ( "h_cosmu_1pip_bg", "", y_bins, 0, 2 );
    TH1D *h_cosmu_1pim_bg   = new TH1D ( "h_cosmu_1pim_bg", "", y_bins, 0, 2 );
    TH1D *h_cosmu_1pi0_bg   = new TH1D ( "h_cosmu_1pi0_bg", "", y_bins, 0, 2 );
    TH1D *h_cosmu_2cpi_bg   = new TH1D ( "h_cosmu_2cpi_bg", "", y_bins, 0, 2 );
    TH1D *h_cosmu_n1cpi_bg  = new TH1D ( "h_cosmu_n1cpi_bg", "", y_bins, 0, 2 );
    TH1D *h_cosmu_other_bg  = new TH1D ( "h_cosmu_other_bg", "", y_bins, 0, 2 );
    TH1D *h_cosmu_nother_bg = new TH1D ( "h_cosmu_nother_bg", "", y_bins, 0, 2 );
    
    signal_c_h.push_back(h_cosmu_mp_s);
    signal_c_h.push_back(h_cosmu_3p_s);
    signal_c_h.push_back(h_cosmu_2p_s);
    signal_c_h.push_back(h_cosmu_1p_s);
    signal_c_h.push_back(h_cosmu_0p_s);

    background_c_h.push_back(h_cosmu_nother_bg);
    background_c_h.push_back(h_cosmu_n1cpi_bg);
    background_c_h.push_back(h_cosmu_other_bg);
    background_c_h.push_back(h_cosmu_2cpi_bg);
    background_c_h.push_back(h_cosmu_1pi0_bg);
    background_c_h.push_back(h_cosmu_1pim_bg);
    background_c_h.push_back(h_cosmu_1pip_bg);
    
    leg_cos_s->AddEntry( h_cosmu_0p_s, " #nu_{#mu} CC0#pi, 0 p ", "f" );
    leg_cos_s->AddEntry( h_cosmu_1p_s, " #nu_{#mu} CC0#pi, 1 p ", "f" );
    leg_cos_s->AddEntry( h_cosmu_2p_s, " #nu_{#mu} CC0#pi, 2 p ", "f" );
    leg_cos_s->AddEntry( h_cosmu_3p_s, " #nu_{#mu} CC0#pi, 3 p ", "f" );
    leg_cos_s->AddEntry( h_cosmu_mp_s, " #nu_{#mu} CC0#pi, >3 p ", "f" );
   
    leg_cos_bg->AddEntry( h_cosmu_1pip_bg, " #nu_{#mu} CC1#pi^{ + } ", "f" );
    leg_cos_bg->AddEntry( h_cosmu_1pim_bg, " #nu_{#mu} CC1#pi^{ - } ", "f" );
    leg_cos_bg->AddEntry( h_cosmu_1pi0_bg, " #nu_{#mu} CC1#pi^{ 0 } ", "f" );
    leg_cos_bg->AddEntry( h_cosmu_2cpi_bg, " #nu_{#mu} CC2#pi^{ #pm } ", "f" );
    leg_cos_bg->AddEntry( h_cosmu_other_bg, " #nu_{#mu} CCOther ", "f" );
    leg_cos_bg->AddEntry( h_cosmu_n1cpi_bg, " #nu_{#mu} NC1#pi^{ #pm } ", "f" );
    leg_cos_bg->AddEntry( h_cosmu_nother_bg, " #nu_{#mu} NCOther ", "f" );
    
    leg_cos_sbg->AddEntry( h_cosmu_sm, " Signal ", "f" );
    leg_cos_sbg->AddEntry( h_cosmu_bg, " Background ", "f" );
    
    leg_cos_sbg_split->AddEntry( h_cosmu_0p_s, " #nu_{#mu} CC0#pi, 0 p ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_1p_s, " #nu_{#mu} CC0#pi, 1 p ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_2p_s, " #nu_{#mu} CC0#pi, 2 p ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_3p_s, " #nu_{#mu} CC0#pi, 3 p ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_mp_s, " #nu_{#mu} CC0#pi, >3 p ", "f" ); 
    leg_cos_sbg_split->AddEntry( h_cosmu_1pip_bg, " #nu_{#mu} CC1#pi^{ + } ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_1pim_bg, " #nu_{#mu} CC1#pi^{ - } ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_1pi0_bg, " #nu_{#mu} CC1#pi^{ 0 } ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_2cpi_bg, " #nu_{#mu} CC2#pi^{ #pm } ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_other_bg, " #nu_{#mu} CCOther ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_n1cpi_bg, " #nu_{#mu} NC1#pi^{ #pm } ", "f" );
    leg_cos_sbg_split->AddEntry( h_cosmu_nother_bg, " #nu_{#mu} NCOther ", "f" );
    
    leg_cos_s->SetHeader("True final state");
    leg_cos_bg->SetHeader("True final state");
    leg_cos_sbg->SetHeader("True final state");
    leg_cos_sbg_split->SetHeader("True final state");

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
        stringstream conv_s_c;
        conv_s_c.clear();

        string title_s_c;
        title_s_c.clear();

        char file_name_s_c[1024];

        conv_s_c << setprecision(4) << "Signal_cosmu_slice_" << low_edge_cos << "_" << up_edge_cos << ".png";
        title_s_c = conv_s_c.str();
        
        strcpy( file_name_s_c, title_s_c.c_str() );

        // Filenames for background histograms
        stringstream conv_bg_c;
        conv_bg_c.clear();

        string title_bg_c;
        title_bg_c.clear();

        char file_name_bg_c[1024];

        conv_bg_c << setprecision(4) << "BG_cosmu_slice_" << low_edge_cos << "_" << up_edge_cos << ".png";
        title_bg_c = conv_bg_c.str();
        
        strcpy( file_name_bg_c, title_bg_c.c_str() );
        
        // Filenames for signal and background histograms
        stringstream conv_sbg_c;
        conv_sbg_c.clear();

        string title_sbg_c;
        title_sbg_c.clear();

        char file_name_sbg_c[1024];

        conv_sbg_c << setprecision(4) << "SBG_cosmu_slice_" << low_edge_cos << "_" << up_edge_cos << ".png";
        title_sbg_c = conv_sbg_c.str();
        
        strcpy( file_name_sbg_c, title_sbg_c.c_str() );
        
        // Filenames for signal and background split histograms
        stringstream conv_sbg_split_c;
        conv_sbg_split_c.clear();

        string title_sbg_split_c;
        title_sbg_split_c.clear();

        char file_name_sbg_split_c[1024];

        conv_sbg_split_c << setprecision(4) << "Split_SBG_cosmu_slice_" << low_edge_cos << "_" << up_edge_cos << ".png";
        title_sbg_split_c = conv_sbg_split_c.str();
        
        strcpy( file_name_sbg_split_c, title_sbg_split_c.c_str() );
        
        // For the title of the histogram
        stringstream hist_s_c;
        hist_s_c.clear();

        char hist_name_s_c[1024];
        string temp_s_c;

        hist_s_c << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "_" << up_edge_cos << ", Signal ";
        temp_s_c = hist_s_c.str();

        strcpy( hist_name_s_c, temp_s_c.c_str() );

        // For the title of the histogram
        stringstream hist_bg_c;
        hist_bg_c.clear();

        char hist_name_bg_c[1024];
        string temp_bg_c;

        hist_bg_c << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "_" << up_edge_cos << ", Background ";
        temp_bg_c = hist_bg_c.str();

        strcpy( hist_name_bg_c, temp_bg_c.c_str() );

        // For the title of the histogram
        stringstream hist_sbg_c;
        hist_sbg_c.clear();

        char hist_name_sbg_c[1024];
        string temp_sbg_c;

        hist_sbg_c << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "_" << up_edge_cos << ", Signal and Background";
        temp_sbg_c = hist_sbg_c.str();

        strcpy( hist_name_sbg_c, temp_sbg_c.c_str() );

        // For the title of the histogram
        stringstream hist_sbg_split_c;
        hist_sbg_split_c.clear();

        char hist_name_sbg_split_c[1024];
        string temp_sbg_split_c;

        hist_sbg_split_c << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "_" << up_edge_cos << ", Signal and Background Split";
        temp_sbg_split_c = hist_sbg_split_c.str();

        strcpy( hist_name_sbg_split_c, temp_sbg_split_c.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= y_bins; ++j ){
            h_cosmu_sm->SetBinContent( j, h_smeared->GetBinContent(i, j) );
        }

        h_cosmu_bg->Reset("M");

        for( unsigned int i = 0; i < signal_c_h.size(); ++i ){
            signal_c_h[i]->Reset("M");
        }  

        for( unsigned int i = 0; i < background_c_h.size(); ++i ){
            background_c_h[i]->Reset("M");
        }  
        
        THStack *hscos_s         = new THStack("hscos_s","Signal histograms");
        THStack *hscos_bg        = new THStack("hscos_bg","Backgound histograms");
        THStack *hscos_sbg       = new THStack("hscos_sbg","Total Signal and Background histograms");
        THStack *hscos_sbg_split = new THStack("hscos_sbg_split","Signal and Background histograms");
            
        for ( int k = 0; k < n_entries; ++k ){
       
            info_s->GetEntry(k);

            double cth_s = l_s_cth->GetValue();
            double T_s   = l_s_T->GetValue();
            int ev_s     = l_s_ev->GetValue();
            int cc_s     = l_s_cc->GetValue();
            int nc_s     = l_s_nc->GetValue();
            int np_s     = l_s_np->GetValue();
            int npip_s   = l_s_npip->GetValue();
            int npim_s   = l_s_npim->GetValue();
            int npi0_s   = l_s_npi0->GetValue();

            if ( cth_s >= low_edge_cos 
                      && cth_s < up_edge_cos
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s == 0 ) {
             
                h_cosmu_0p_s->Fill(T_s,1);
            }
            else if ( cth_s >= low_edge_cos 
                      && cth_s < up_edge_cos
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s == 1 ) {
             
                h_cosmu_1p_s->Fill(T_s,1);
            }
            else if ( cth_s >= low_edge_cos 
                      && cth_s < up_edge_cos
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s == 2 ) {
             
                h_cosmu_2p_s->Fill(T_s,1);
            }
            else if ( cth_s >= low_edge_cos 
                      && cth_s < up_edge_cos
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s == 3 ) {
             
                h_cosmu_3p_s->Fill(T_s,1);
            }
            else if ( cth_s >= low_edge_cos 
                      && cth_s < up_edge_cos
                      && cc_s
                      && ( npip_s + npim_s + npi0_s == 0 )
                      && np_s > 3 ) {
             
                h_cosmu_mp_s->Fill(T_s,1);
            }
        }

        for ( int k = 0; k < bg_entries; ++k ){
       
            info_bg->GetEntry(k);

            double cth_bg = l_bg_cth->GetValue();
            double T_bg   = l_bg_T->GetValue();
            int ev_bg     = l_bg_ev->GetValue();
            int cc_bg     = l_bg_cc->GetValue();
            int nc_bg     = l_bg_nc->GetValue();
            int np_bg     = l_bg_np->GetValue();
            int npip_bg   = l_bg_npip->GetValue();
            int npim_bg   = l_bg_npim->GetValue();
            int npi0_bg   = l_bg_npi0->GetValue();
            
            if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && cc_bg
                      && npip_bg  == 1  
                      && ( npim_bg + npi0_bg == 0 ) ){
                
                h_cosmu_1pip_bg->Fill(T_bg,1);
            }
            if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && cc_bg
                      && npim_bg  == 1  
                      && ( npip_bg + npi0_bg == 0 ) ){
                
                h_cosmu_1pim_bg->Fill(T_bg,1);
            }
            if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && cc_bg
                      && npi0_bg  == 1  
                      && ( npip_bg + npim_bg == 0 ) ){
                
                h_cosmu_1pi0_bg->Fill(T_bg,1);
            }
            else if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && cc_bg
                      && ( npip_bg + npim_bg == 2 ) 
                      && npi0_bg == 0 ){
                
                h_cosmu_2cpi_bg->Fill(T_bg,1);
            }
            else if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && nc_bg
                      && ( npip_bg + npim_bg == 1 )
                      && npi0_bg == 0 ) {
                
                h_cosmu_n1cpi_bg->Fill(T_bg,1);
            }
            else if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && cc_bg ) {
                 
                h_cosmu_other_bg->Fill(T_bg,1);
            }
            else if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && nc_bg ) {
                 
                h_cosmu_nother_bg->Fill(T_bg,1);
            }
        }

        for ( int k = 0; k < bg_entries; ++k ){
       
            info_bg->GetEntry(k);

            double cth_bg = l_bg_cth->GetValue();
            double T_bg   = l_bg_T->GetValue();
            int ev_bg     = l_bg_ev->GetValue();
            int cc_bg     = l_bg_cc->GetValue();
            int nc_bg     = l_bg_nc->GetValue();
            int npip_bg   = l_bg_npip->GetValue();
            int npim_bg   = l_bg_npim->GetValue();
            int npi0_bg   = l_bg_npi0->GetValue();
            
            if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && cc_bg
                      && ( npip_bg + npim_bg == 1 ) 
                      && npi0_bg == 0 ){
                
                h_cosmu_bg->Fill(T_bg,1);
            }
            else if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && cc_bg
                      && ( npip_bg + npim_bg == 2 ) 
                      && npi0_bg == 0 ){
                
                h_cosmu_bg->Fill(T_bg,1);
            }
            else if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && nc_bg
                      && ( npip_bg + npim_bg == 1 )
                      && npi0_bg == 0 ) {
                
                h_cosmu_bg->Fill(T_bg,1);
            }
            else if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && cc_bg ) {
                 
                h_cosmu_bg->Fill(T_bg,1);
            }
            else if ( cth_bg >= low_edge_cos 
                      && cth_bg < up_edge_cos
                      && nc_bg ) {
                 
                h_cosmu_bg->Fill(T_bg,1);
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

        h_cosmu_sm->SetFillColor( kGreen + 1 );
        h_cosmu_bg->SetFillColor( kRed + 1 );
        
        h_cosmu_bg->SetFillStyle( 3354 );
        
        double norm_cosmu = h_cosmu_bg->Integral() + h_cosmu_sm->Integral();
        
        h_cosmu_sm->Scale(1/norm_cosmu);
        h_cosmu_bg->Scale(1/norm_cosmu);
        
        for ( unsigned int i = 0; i < background_c_h.size(); ++i ){

            int s = signal_c_h.size();

            background_c_h[i]->SetFillColor(pal[s + i]);
            background_c_h[i]->SetLineColor(pal[s + i]);

            background_c_h[i]->SetFillStyle( 3354 );
            background_c_h[i]->SetLineWidth(1.5);
            background_c_h[i]->Scale(1/norm_cosmu);

            hscos_bg->Add(background_c_h[i]);
            hscos_sbg_split->Add(background_c_h[i]);
            
        }
        for ( unsigned int i = 0; i < signal_c_h.size(); ++i ){

            signal_c_h[i]->SetFillColor(pal[i]);
            signal_c_h[i]->SetLineColor(pal[i]);

            signal_c_h[i]->SetLineWidth(1.5);
            signal_c_h[i]->Scale(1/norm_cosmu);

            hscos_s->Add(signal_c_h[i]);
            hscos_sbg_split->Add(signal_c_h[i]);
            
        }
        
        h_cosmu_sm->SetLineColor( kGreen + 1 );
        h_cosmu_bg->SetLineColor( kRed + 1 );

        h_cosmu_sm->SetLineWidth(1.5);     
        h_cosmu_bg->SetLineWidth(1.5);     
        
        // Fill the stacks 
        // Signal
        hscos_s->Draw();
        
        hscos_s->SetTitle(hist_name_s_c);
        hscos_s->GetXaxis()->SetTitle("T_{#mu}");   
        hscos_s->GetYaxis()->SetTitle("Number of SBND events");   
    
        leg_cos_s->Draw();
        c_cosmu->SaveAs(file_name_s_c);

        // Background
        hscos_bg->Draw();
        
        hscos_bg->SetTitle(hist_name_bg_c);
        hscos_bg->GetXaxis()->SetTitle("T_{#mu}");   
        hscos_bg->GetYaxis()->SetTitle("Number of SBND events");   
        
        leg_cos_bg->Draw();
        c_cosmu->SaveAs(file_name_bg_c);

        // Signal and background
        hscos_sbg->Add(h_cosmu_bg);
        hscos_sbg->Add(h_cosmu_sm);
        
        hscos_sbg->Draw();
        
        hscos_sbg->SetTitle(hist_name_sbg_c);
        hscos_sbg->GetXaxis()->SetTitle("T_{#mu}");   
        hscos_sbg->GetYaxis()->SetTitle("Number of SBND events");   

        leg_cos_sbg->Draw();
        c_cosmu->SaveAs(file_name_sbg_c);

        // Signal and background split
        hscos_sbg_split->Draw();
        
        hscos_sbg_split->SetTitle(hist_name_sbg_split_c);
        hscos_sbg_split->GetXaxis()->SetTitle("T_{#mu}");   
        hscos_sbg_split->GetYaxis()->SetTitle("Number of SBND events");   
        
        leg_cos_sbg_split->Draw();
        c_cosmu->SaveAs(file_name_sbg_split_c);

        delete hscos_s;
        delete hscos_sbg;
        delete hscos_bg;
        delete hscos_sbg_split;
    } 

    delete h_cosmu_sm;
    delete h_cosmu_bg;

    for ( unsigned int i = 0; i < signal_c_h.size(); ++i ){
        delete signal_c_h[i];
    }
    
    for ( unsigned int i = 0; i < background_c_h.size(); ++i ){
        delete background_c_h[i];
    }
    
    delete c_cosmu;
    
    delete leg_cos_s;
    delete leg_cos_bg;
    delete leg_cos_sbg;
    delete leg_cos_sbg_split;
    

}

// ---------------------------------------------------------------------------
//                          Efficiency function
// ---------------------------------------------------------------------------

void Efficiency ( TH2D *h_smeared, TH2D *h_unsmeared, TH2D *h_cut, TNtuple *info_s){

    int n_entries = info_s->GetEntries();

    // Get the  signal leaves
    TLeaf *l_s_ev   = info_s->GetLeaf( "ev" );
    TLeaf *l_s_cth  = info_s->GetLeaf( "cth" );
    TLeaf *l_s_T    = info_s->GetLeaf( "T" );
    TLeaf *l_s_cc   = info_s->GetLeaf( "cc" );
    TLeaf *l_s_qel  = info_s->GetLeaf( "qel" );
    TLeaf *l_s_res  = info_s->GetLeaf( "res" );
    TLeaf *l_s_npip = info_s->GetLeaf( "npip" );
    TLeaf *l_s_npim = info_s->GetLeaf( "npim" );
    TLeaf *l_s_npi0 = info_s->GetLeaf( "npi0" );
    
    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_smeared->GetNbinsX(); // Cos theta
    int y_bins = h_smeared->GetNbinsY(); // Tmu

    TCanvas *c  = new TCanvas( "", "", 800, 600 );

    TH2D *h_eff    = new TH2D( "h_eff"," Efficiency ", x_bins, -1, 1, y_bins, 0, 2 );
    TH2D *h_pur    = new TH2D( "h_pur"," Purity ", x_bins, -1, 1, y_bins, 0, 2 );
    TH2D *h_impur  = new TH2D( "h_impur"," Impurity ", x_bins, -1, 1, y_bins, 0, 2 );
    TH2D *h_effpur = new TH2D( "h_effpur"," Efficiency * Purity ", x_bins, -1, 1, y_bins, 0, 2 );

    vector< double > T_bin_edges;
    vector< double > cos_bin_edges;
    
    for ( int i = 1; i <= y_bins; ++i ){

        // Get the lower and upper bin edges of the y axis
        T_bin_edges.push_back( h_smeared->GetYaxis()->GetBinLowEdge( i ) );
 
    }

    T_bin_edges.push_back( h_smeared->GetYaxis()->GetBinLowEdge( y_bins + 1 ) );
     
    for ( int i = 1; i <= x_bins; ++i ){

        // Get the lower and upper bin edges of the y axis
        cos_bin_edges.push_back( h_smeared->GetXaxis()->GetBinLowEdge( i ) );
 
    }

    cos_bin_edges.push_back( h_smeared->GetXaxis()->GetBinLowEdge( x_bins + 1 ) );
  
    for ( int i = 0; i < y_bins; ++i ){

        for( int j = 0; j < x_bins; ++j ){

            // Initialise some counters, these are started inside the loop
            // so that they are reset for each bin
            int ccqe_signal_count  = 0;
            int ccqe_cut_count     = 0;
            int ccqe_truth_count   = 0;
            int total_signal_count = 0;

            for ( int k = 0; k < n_entries; ++k ){
           
                info_s->GetEntry(k);
    
                int ev_s     = l_s_ev->GetValue();
                double cth_s = l_s_cth->GetValue();
                double T_s   = l_s_T->GetValue();
                int cc_s     = l_s_cc->GetValue();
                int qel_s    = l_s_qel->GetValue();
                int res_s    = l_s_res->GetValue();
                int npip_s   = l_s_npip->GetValue();
                int npim_s   = l_s_npim->GetValue();
                int npi0_s   = l_s_npi0->GetValue();
               
                if ( T_s >= T_bin_edges[i] 
                     && T_s < T_bin_edges[i+1]
                     && cth_s >= cos_bin_edges[j]
                     && cth_s < cos_bin_edges[j+1]
                     && cc_s
                     && qel_s
                     && ( npip_s + npim_s + npi0_s == 0 ) ){

                    // This is a CCQE signal event in the i,j bin
                    ccqe_signal_count++;
                 
                }
                if ( T_s >= T_bin_edges[i] 
                     && T_s < T_bin_edges[i+1]
                     && cth_s >= cos_bin_edges[j]
                     && cth_s < cos_bin_edges[j+1] ){

                    // This is a signal event in the i,j bin
                    total_signal_count++;
                 
                }
            }

            ccqe_truth_count = h_unsmeared->GetBinContent( j + 1, i + 1 );
            ccqe_cut_count   = h_cut->GetBinContent( j + 1, i + 1 );

            double purity, impurity, efficiency, effpur;
            
            if( total_signal_count != 0 ){
                
                // Fill the purity and impurity histograms
                purity   = ( ccqe_signal_count / double( total_signal_count ) );
                impurity = 1 - purity;
                
                h_pur->SetBinContent( j + 1, i + 1, purity );
                h_impur->SetBinContent( j + 1, i + 1, impurity );
            }

            if ( ccqe_truth_count != 0 ){

                // Fill the efficiency histogram
                efficiency = ( ccqe_cut_count / double( ccqe_truth_count ) );

                h_eff->SetBinContent( j + 1, i + 1,  efficiency );
            }

            if ( ccqe_truth_count != 0 && total_signal_count != 0 ){
                
                // Fill the efficiency * purity histogram
                effpur = efficiency * purity;
                h_effpur->SetBinContent( j + 1, i + 1, effpur ); 
            }
        }
    }
   
    gStyle->SetPalette(55);

    h_pur->Draw("colz");
    h_pur->SetStats(kFALSE);
    h_pur->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_pur->GetYaxis()->SetTitle("T_{#mu}");
    h_pur->GetZaxis()->SetRangeUser(0.2,1);
    c->Update();
    //c->SetLogz();
    c->SaveAs("purity.png");
    //c->SaveAs("log_purity.png");
    
    h_eff->Draw("colz");
    h_eff->SetStats(kFALSE);
    h_eff->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_eff->GetYaxis()->SetTitle("T_{#mu}");
    h_eff->GetZaxis()->SetRangeUser(0.5,1);
    c->Update();
    //c->SetLogz();
    c->SaveAs("efficiency.png");
    //c->SaveAs("log_efficiency.png");
    
    h_effpur->Draw("colz");
    h_effpur->SetStats(kFALSE);
    h_effpur->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_effpur->GetYaxis()->SetTitle("T_{#mu}");
    h_effpur->GetZaxis()->SetRangeUser(0.1,1);
    c->Update();
    //c->SetLogz();
    c->SaveAs("effpur.png");
    //c->SaveAs("log_effpur.png");
    
    h_impur->Draw("colz");
    h_impur->SetStats(kFALSE);
    h_impur->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_impur->GetYaxis()->SetTitle("T_{#mu}");
    c->Update();
    c->SetLogz();
    //c->SaveAs("impurity.png");
    c->SaveAs("log_impurity.png");
    
    delete c;

    delete h_eff;
    delete h_pur;
    delete h_impur;
    delete h_effpur;

}

