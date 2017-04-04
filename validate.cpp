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

#include "validate.h"

using namespace std; 

int validate() {

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
    
    // Number of MB bins (data) 
    // costhetamu : 20
    // Tmu        : 18
    // 
    // Amount to smear 
    // costhetamu : 5%
    // Tmu        : 10%
    //
    // Kinetic energy cuts
    // mu         : > 50MeV
    // charged pi : > 50MeV
   
    // -------------------------------------------------------------------------
    //                  Get the tree, normalisation and the branches
    // -------------------------------------------------------------------------
    TTree *def_tree = (TTree*) f1.Get( "gst" );
    
    TCanvas *c1 = new TCanvas("c1","",1000,800);
    //TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} distribution before smearing ",20,-1,1,18,0,2);
    TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} distribution before smearing, log ",20,-1,1,18,0,2);
    def_tree->Draw("pl:cthl>>h_un","fspl == 13 && cc && (nfpip + nfpim + nfpi0 == 0)","colz"); 
   
    c1->SetRightMargin(0.13);

    h_un->SetStats(kFALSE);
    h_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un->GetYaxis()->SetTitle("T_{#mu}");
    c1->SetLogz();
    //c1->SaveAs("distribution_before.png");
    c1->SaveAs("distribution_before_log.png");
    
    //TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, impurity ",20,-1,1,18,0,2);
    TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, log, impurity ",20,-1,1,18,0,2);
    
    // Take h_un and smear it
    Smear(def_tree, h_un, h_sm, 0.05, 0.1, 0.05, 0.05, 0.2);

    // The the 2 2D histograms and draw slices in Tmu and cos theta mu
    Slices( h_un, h_sm );

    return 0;
}

// -------------------------------------------------------------------------
//                            Function definitions
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
//                             Smearing function
// -------------------------------------------------------------------------

void Smear ( TTree *tree, 
             TH2D *h_unsmeared,
             TH2D *h_smeared,
             double x_smear, 
             double y_smear, 
             double E_mu_cut,
             double E_pi_cut,
             double mu_pi_switch ){

    // File for writing lots to
    ofstream file;
    file.open("out_file.txt"); 

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
    TBranch *b_cthl  = tree->GetBranch( "cthl" );
    TBranch *b_El    = tree->GetBranch( "El" );
    TBranch *b_pl    = tree->GetBranch( "pl" );
    TBranch *b_fspl  = tree->GetBranch( "fspl" );
    TBranch *b_cthf  = tree->GetBranch( "cthf" );
    TBranch *b_Ef    = tree->GetBranch( "Ef" );
    TBranch *b_pf    = tree->GetBranch( "pf" );
    TBranch *b_pdgf  = tree->GetBranch( "pdgf" );
    TBranch *b_cc    = tree->GetBranch( "cc" );
    TBranch *b_nc    = tree->GetBranch( "nc" );
    TBranch *b_nfpip = tree->GetBranch( "nfpip" );
    TBranch *b_nfpim = tree->GetBranch( "nfpim" );
    TBranch *b_nfpi0 = tree->GetBranch( "nfpi0" );

    int n_values = tree->GetEntries();

    int e_count = 0;
    int other_count = 0;
    int pion_count = 0;

    vector< double > T_mu_vect;
    vector< double > T_pi_vect;
    vector< int > Impurity;
   
    int impurity_count = 0;
    int muon_without_cuts = 0;
    int all_cuts_count = 0;

    // Loop over the branches and calculate the kinetic energy
    // print a few of them, along with the momentum to find out the difference
    for ( int i = 0; i < n_values; ++i ){
        
        tree->GetEntry(i);
     
        double T_mu, T_pi, e_mu, e_pi;
        double m_mu = 0.10566; // Muon mass, GeV
        double m_pi = 0.13957; // Charged pion mass, GeV

        int pdgf = b_pdgf->GetLeaf("pdgf")->GetValue();
        int fspl = b_fspl->GetLeaf("fspl")->GetValue();

        // Calculate the kinetic energy for muons
        if ( fspl == 13 ){  
         
            e_mu = b_El->GetLeaf("El")->GetValue();
            T_mu = sqrt( pow(e_mu,2) - pow(m_mu,2) );

            T_mu_vect.push_back(T_mu);

            other_count++;
        }
        else if ( fspl == 11 ){
            T_mu_vect.push_back(-99999);
            e_count++;
        }
        else{
            T_mu_vect.push_back(-99999);
            other_count++;
        }

        // Calculate the kinetic energy of the pions
        if ( pdgf == 211 || pdgf == -211 ){

            e_pi = b_Ef->GetLeaf("Ef")->GetValue();
            T_pi = sqrt( pow(e_pi,2) - pow(m_pi,2) );

            T_pi_vect.push_back(T_pi);

            pion_count++;   
        }
        else{
            T_pi_vect.push_back(-99999); 
        }

        // Apply the impurity cut
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 && (b_nfpip->GetLeaf( "nfpip" )->GetValue() + b_nfpim->GetLeaf( "nfpim" )->GetValue() == 1 )){
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
       
    // Implement the smearing factors and plot the histograms both before and after to observe the effect
    // Let's start with Gaussian smearing for now and move onto log normal smearing later
    // Now try the smearing
    TCanvas *c2 = new TCanvas("c2","",1000,800);


    // Fill h_smeared normally and with the impurities and THEN clone it
    for ( int i = 0; i < n_values; ++i ){        
        tree->GetEntry(i);

        double fspl  = b_fspl->GetLeaf( "fspl" )->GetValue(); 
        double pdgf  = b_pdgf->GetLeaf( "pdgf" )->GetValue(); 
        double cc    = b_cc->GetLeaf( "cc" )->GetValue(); 
        double nc    = b_nc->GetLeaf( "nc" )->GetValue(); 
        double nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        double nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        double nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();
        double Ef    = b_Ef->GetLeaf( "Ef" )->GetValue();
        double pf    = b_pf->GetLeaf( "pf" )->GetValue();
        double cthf  = b_cthf->GetLeaf( "cthf" )->GetValue();
        double El    = b_El->GetLeaf( "El" )->GetValue();
        double pl    = b_pl->GetLeaf( "pl" )->GetValue();
        double cthl  = b_cthl->GetLeaf( "cthl" )->GetValue();

        // Set the energy and impurity cuts
        if ( fspl == 13 
            && cc != 0 
            && (nfpip + nfpim + nfpi0 == 0) 
            && T_mu_vect[i] > 0.05 ){
                
            h_smeared->Fill(cthl, pl);
        }
        else if ( ( pdgf == 211 || pdgf == -211 ) 
            && Impurity[i] == 1 
            && nc != 0
            && (nfpip + nfpim == 1)
            && T_pi_vect[i] > 0.05 ){

            h_smeared->Fill(cthf, pf);
        }
    }

    double x, y, w, z;

    double sX = h_unsmeared->GetStdDev(1) * ( 1 + x_smear ); // Bin width plus 5%
    double sY = h_unsmeared->GetStdDev(2) * ( 1 + y_smear ); // Bin width plus 10%

    TH2D *h3 = (TH2D*) h_smeared->Clone();
    h_smeared->Reset();

    double nX = h_smeared->GetNbinsX();
    double nY = h_smeared->GetNbinsY();
  
    // Loop over x and y bins to fill with the smeared distribution
    for ( int j = 1; j <= nX; ++j ){
        x = h_smeared->GetXaxis()->GetBinCenter(j);            
                
        for ( int k = 1; k <= nY; ++k ){
            y = h_smeared->GetYaxis()->GetBinCenter(k);            
            double content = 0;
                    
            for ( int l = 1; l <= nX; ++l ){
                w = h3->GetXaxis()->GetBinCenter(l);

                for ( int m = 1; m <= nY; m++ ){
                    z = h3->GetYaxis()->GetBinCenter(m);
                    content += h3->GetBinContent(j,k) * TMath::Gaus(x,w,sX,false) * TMath::Gaus(y,z,sY,false);
                }
            }
            content = TMath::Floor( content * h_smeared->GetXaxis()->GetBinWidth(j) * h_smeared->GetYaxis()->GetBinWidth(k) / 2. / TMath::Pi() / sX / sY );
            h_smeared->SetBinContent(j,k, content);
        }
    }

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

    /*
    // A histogram for the differences
    TCanvas *c4 = new TCanvas("c4","",1000,800);
    TH2D *h4 = new TH2D("h4"," The absolute differences between the smeared and unsmeared bins ",20,-1,1,18,0,2);
  
    //gStyle->SetPalette(kBird);

    double d1, d2, d4 = 0;

    for ( int i = 1; i <= nX; ++i){
        for ( int j = 1; j <= nY; ++j ){
            if ( !(  h_unsmeared->GetBinContent(i,j) != h_unsmeared->GetBinContent(i,j) ) ){
                h4->SetBinContent(i, j, TMath::Abs( ( h_unsmeared->GetBinContent(i,j) - h_smeared->GetBinContent(i,j) ) / h_unsmeared->GetBinContent(i,j) ) );            
            }
            else{
                h4->SetBinContent(i,j,0);    
            }
        }
    }

    h4->Draw("colz");
    h4->SetStats(kFALSE);
    h4->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h4->GetYaxis()->SetTitle("T_{#mu}");
    c4->SaveAs("bin_differences.png");

    */

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
    
    TLegend *leg_c = new TLegend( 0.12, 0.78, 0.28, 0.88 );

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

