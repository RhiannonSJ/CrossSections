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

#include "smearing.h"

using namespace std; 

int smearing() {

    // -------------------------------------------------------------------------
    //                              Open event files
    // -------------------------------------------------------------------------
    
    // Open Default
    TFile f1("/hepstore/rjones/Exercises/Flavours/Default/sbnd/1M/gntp.10000.gst.root");
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
    TCanvas *c1 = new TCanvas("c1","",1000,800);
    TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} distribution before smearing ",20,-1,1,18,0,2);
    //TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} distribution before smearing, log ",20,-1,1,18,0,2);
    def_tree->Draw("( El - 0.10566 ):cthl>>h_un","fspl == 13 && cc && (nfpip + nfpim + nfpi0 == 0)","colz"); 
   
    c1->SetRightMargin(0.13);

    h_un->SetStats(kFALSE);
    h_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un->GetYaxis()->SetTitle("T_{#mu}");
    //c1->SetLogz();
    c1->SaveAs("distribution_before.png");
    //c1->SaveAs("distribution_before_log.png");
    
    // The same histogram definitions for the smeared distributions
    TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, impurity ",20,-1,1,18,0,2);
    //TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, log, impurity ",20,-1,1,18,0,2);
    
    // Take h_un and smear it
    Smear(def_tree, h_un, h_sm);

    // The the 2 2D histograms and draw slices in Tmu and cos theta mu
    // Slices( h_un, h_sm );

    // The the 2 2D histograms and draw slices in Tmu and cos theta mu
    Characterisation( h_sm, def_tree );
    
    return 0;
}

// -------------------------------------------------------------------------
//                            Function definitions
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
//                             Smearing function
// -------------------------------------------------------------------------

void Smear ( TTree *tree, 
             TH2D  *h_unsmeared,
             TH2D  *h_smeared ){

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

        double T_mu, e_mu;

        // Calculate the kinetic energy for muons
        if ( fspl == 13 ){  
         
            // Energy of the final state primary lepton
            e_mu = b_El->GetLeaf("El")->GetValue();
            T_mu = e_mu - m_mu;

            T_mu_vect.push_back(T_mu);

        }
        // If the final state primary is a lepton, push back a number that will
        // be removed in the cuts later
        else if ( fspl == 11 ){
            T_mu_vect.push_back(-99999);
        }
        else{
            T_mu_vect.push_back(-99999);
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
    TCanvas *c2 = new TCanvas("c2","",1000,800);


    // Fill h_smeared normally and with the impurities and THEN clone it
    for ( int i = 0; i < n_values; ++i ){        
        tree->GetEntry(i);

        double cc    = b_cc->GetLeaf( "cc" )->GetValue(); 
        double fspl  = b_fspl->GetLeaf( "fspl" )->GetValue(); 
        double nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        double nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        double nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();

        // Set the energy and impurity cuts
        if ( fspl == 13 
            && cc != 0 
            && (nfpip + nfpim + nfpi0 == 0) 
            && T_mu_vect[i] > 0.05 ){
                
            h_smeared->Fill(cos_mu_prime[i], T_mu_prime[i]);
        }
        else if ( Impurity[i] == 1 
            && T_pi_vect[i] > 0.05 ){

            h_smeared->Fill(cos_pi_vect[i], T_pi_vect[i]);
        }
    }

    double nX = h_smeared->GetNbinsX();
    double nY = h_smeared->GetNbinsY();
 
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
    //c2->SetLogz();
    c2->SaveAs("distribution_after_impure.png");
    //c2->SaveAs("distribution_after_log_impure.png");


}


void Slices ( TH2D *h_unsmeared, TH2D *h_smeared ){

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_unsmeared->GetNbinsX(); // Cos theta
    int y_bins = h_unsmeared->GetNbinsY(); // Tmu

    TCanvas *c_Tmu   = new TCanvas ( "c_Tmu", "", 800, 600 );
   
    TLegend *leg_T = new TLegend( 0.12, 0.78, 0.28, 0.88 );

    TH1D *h_Tmu      = new TH1D ( "h_Tmu", "", x_bins, -1, 1 );
    TH1D *h_Tmu_sm   = new TH1D ( "h_Tmu_sm", "", x_bins, -1, 1 );
    
    leg_T->AddEntry( h_Tmu, " Unsmeared ", "l" );
    leg_T->AddEntry( h_Tmu_sm, " Smeared ", "l" );
    
    for ( int i = 1; i <= y_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_T;
        double upp_edge_T;

        low_edge_T = h_unsmeared->GetYaxis()->GetBinLowEdge(i);
        upp_edge_T = h_unsmeared->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv;
        conv.clear();

        string title;
        title.clear();

        char file_name[1024];

        conv << setprecision(4) << "Tmu_slice_" << low_edge_T << ";" << upp_edge_T << ".png";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        stringstream hist;
        hist.clear();

        char hist_name[1024];
        string temp1;

        hist << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "-" << upp_edge_T;
        temp1 = hist.str();

        strcpy( hist_name, temp1.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= x_bins; ++j ){
            h_Tmu->SetBinContent( j, h_unsmeared->GetBinContent(j, i) );
            h_Tmu_sm->SetBinContent( j, h_smeared->GetBinContent(j, i) );
        }

        h_Tmu->Draw();
        h_Tmu_sm->Draw("same");
        h_Tmu->SetTitle(hist_name);
        h_Tmu->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        h_Tmu->GetYaxis()->SetTitle("Number of SBND events");   
        h_Tmu->SetLineColor( kRed + 2 );
        h_Tmu_sm->SetLineColor( kGreen + 2 );
        h_Tmu->SetTitleOffset(1.5, "Y");
        h_Tmu->SetStats(kFALSE);

        leg_T->Draw();
        c_Tmu->SaveAs(file_name);

    } 
   
    delete h_Tmu;
    delete h_Tmu_sm;

    delete c_Tmu;
    
    delete leg_T;
    
    TCanvas *c_cosmu = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_c = new TLegend( 0.72, 0.78, 0.88, 0.88 );

    TH1D *h_cosmu    = new TH1D ( "h_cosmu", "", y_bins, 0, 2 );
    TH1D *h_cosmu_sm = new TH1D ( "h_cosmu_sm", "", y_bins, 0, 2 );
     
    leg_c->AddEntry( h_cosmu, " Unsmeared ", "l" );
    leg_c->AddEntry( h_cosmu_sm, " Smeared ", "l" );
    
    // Cos theta mu slices
    for ( int i = 1; i <= x_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_cos;
        double upp_edge_cos;

        low_edge_cos = h_unsmeared->GetXaxis()->GetBinLowEdge(i);
        upp_edge_cos = h_unsmeared->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv1;
        conv1.clear();

        string title1;
        title1.clear();

        char file_name1[1024];

        conv1 << setprecision(4) << "cos_thetamu_slice_" << low_edge_cos << ";" << upp_edge_cos << ".png";
        title1 = conv1.str();
        
        strcpy( file_name1, title1.c_str() );

        // For the title of the histogram
        stringstream hist1;
        hist1.clear();

        char hist_name1[1024];
        string temp2;

        hist1 << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "," << upp_edge_cos;
        temp2 = hist1.str();

        strcpy( hist_name1, temp2.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= y_bins; ++j ){
            h_cosmu->SetBinContent( j, h_unsmeared->GetBinContent(i, j) );
            h_cosmu_sm->SetBinContent( j, h_smeared->GetBinContent(i, j) );
        }

        h_cosmu->Draw();
        h_cosmu_sm->Draw("same");
        h_cosmu->SetTitle(hist_name1);
        h_cosmu->GetXaxis()->SetTitle("T_{#mu}");   
        h_cosmu->GetYaxis()->SetTitle("Number of SBND events");   
        h_cosmu->SetLineColor( kRed + 2 );
        h_cosmu_sm->SetLineColor( kGreen + 2 );
        h_cosmu->SetTitleOffset(1.5, "Y");
        h_cosmu->SetStats(kFALSE);

        leg_c->Draw();
        c_cosmu->SaveAs(file_name1);

    } 
    
    delete h_cosmu;
    delete h_cosmu_sm;

    delete c_cosmu;

    delete leg_c;

}

void Characterisation ( TH2D *h_smeared, TTree *tree ){

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
    TBranch *b_nfpip = tree->GetBranch( "nfpip" );
    TBranch *b_nfpim = tree->GetBranch( "nfpim" );
    TBranch *b_nfpi0 = tree->GetBranch( "nfpi0" );

    // Number of events in the TTree
    int n_values = tree->GetEntries();
    //int n_values = 1000; // Using 100 events for debugging purposes
    
    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_smeared->GetNbinsX(); // Cos theta
    int y_bins = h_smeared->GetNbinsY(); // Tmu

    TCanvas *c_Tmu   = new TCanvas ( "c_Tmu", "", 800, 600 );
   
    TLegend *leg_T = new TLegend( 0.12, 0.78, 0.28, 0.88 );

    TH1D *h_Tmu_sm   = new TH1D ( "h_Tmu_sm", "", x_bins, -1, 1 );
    TH1D *h_Tmu_cc0pi   = new TH1D ( "h_Tmu_cc0pi", "", x_bins, -1, 1 );
    TH1D *h_Tmu_other   = new TH1D ( "h_Tmu_other", "", x_bins, -1, 1 );
    
    leg_T->AddEntry( h_Tmu_sm, " Total signal ", "l" );
    leg_T->AddEntry( h_Tmu_cc0pi, " cc0pi ", "l" );
    leg_T->AddEntry( h_Tmu_other, " other ", "l" );
    
    for ( int i = 1; i <= y_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_T;
        double upp_edge_T;

        low_edge_T = h_smeared->GetYaxis()->GetBinLowEdge(i);
        upp_edge_T = h_smeared->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv;
        conv.clear();

        string title;
        title.clear();

        char file_name[1024];

        conv << setprecision(4) << "Tmu_slice_" << low_edge_T << ";" << upp_edge_T << ".png";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        stringstream hist;
        hist.clear();

        char hist_name[1024];
        string temp1;

        hist << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "-" << upp_edge_T;
        temp1 = hist.str();

        strcpy( hist_name, temp1.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= x_bins; ++j ){
            h_Tmu_sm->SetBinContent( j, h_smeared->GetBinContent(j, i) );
        }

        h_Tmu_sm->Draw();
        h_Tmu_sm->SetTitle(hist_name);
        h_Tmu_sm->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        h_Tmu_sm->GetYaxis()->SetTitle("Number of SBND events");   
        //h_Tmu->SetLineColor( kRed + 2 );
        h_Tmu_sm->SetLineColor( kGreen + 2 );
        h_Tmu_sm->SetTitleOffset(1.5, "Y");
        h_Tmu_sm->SetStats(kFALSE);

        leg_T->Draw();
        c_Tmu->SaveAs(file_name);

    } 
   
    delete h_Tmu_sm;

    delete c_Tmu;
    
    delete leg_T;
    
    TCanvas *c_cosmu = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_c = new TLegend( 0.72, 0.78, 0.88, 0.88 );

    TH1D *h_cosmu_sm = new TH1D ( "h_cosmu_sm", "", y_bins, 0, 2 );
     
    leg_c->AddEntry( h_cosmu_sm, " Total signal ", "l" );
    
    // Cos theta mu slices
    for ( int i = 1; i <= x_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_cos;
        double upp_edge_cos;

        low_edge_cos = h_smeared->GetXaxis()->GetBinLowEdge(i);
        upp_edge_cos = h_smeared->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv1;
        conv1.clear();

        string title1;
        title1.clear();

        char file_name1[1024];

        conv1 << setprecision(4) << "cos_thetamu_slice_" << low_edge_cos << ";" << upp_edge_cos << ".png";
        title1 = conv1.str();
        
        strcpy( file_name1, title1.c_str() );

        // For the title of the histogram
        stringstream hist1;
        hist1.clear();

        char hist_name1[1024];
        string temp2;

        hist1 << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "," << upp_edge_cos;
        temp2 = hist1.str();

        strcpy( hist_name1, temp2.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= y_bins; ++j ){
            h_cosmu_sm->SetBinContent( j, h_smeared->GetBinContent(i, j) );
        }

        h_cosmu_sm->Draw();
        h_cosmu_sm->SetTitle(hist_name1);
        h_cosmu_sm->GetXaxis()->SetTitle("T_{#mu}");   
        h_cosmu_sm->GetYaxis()->SetTitle("Number of SBND events");   
        //h_cosmu->SetLineColor( kRed + 2 );
        h_cosmu_sm->SetLineColor( kGreen + 2 );
        h_cosmu_sm->SetTitleOffset(1.5, "Y");
        h_cosmu_sm->SetStats(kFALSE);

        leg_c->Draw();
        c_cosmu->SaveAs(file_name1);

    } 
    
    delete h_cosmu_sm;

    delete c_cosmu;

    delete leg_c;
}

