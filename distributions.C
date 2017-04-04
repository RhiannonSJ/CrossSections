// A quick macro to test the root LogNormal function
// Using 100 values
// Random number generation

#include "smearing.h"

int distributions(){

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();    
    _random_gen->SetSeed( time( NULL ) ); 

    TCanvas *c = new TCanvas( "c", "", 800, 600 );
    TH1D *h    = new TH1D( "h", "Log-Normal distribution example", 50, 0, 50 );

    for ( int i = 0; i < 100000; ++i ){
        
        // Calculate the mean and sigma for the LogNormal function
        // zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
        // sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
        // m     = expectation value
        // var   = variance = s.d.^2 = ( expectation value * 0.1 ) ^ 2

        double  m     = 10;
        double  var   = pow( 10 * 0.5, 2 ); 
        
        double sigma = TMath::Sqrt( TMath::Log( 1 + ( var / TMath::Power( m, 2 ) ) ) );
        double zeta  = TMath::Log( m * ( 1 / TMath::Sqrt( 1 + ( var / TMath::Power( m, 2 ) ) ) ) );

        double lognorm = _random_gen->LogNormal( zeta, sigma );
        
        h->Fill(lognorm);

    }   

    h->Draw("l");
    h->SetStats(kFALSE);    
    c->SaveAs( "log_norm.png" );

    TCanvas *c1 = new TCanvas( "c1", "", 800, 600 );
    TH1D *h1    = new TH1D( "h1", "Gaussian distribution example", 50, -20, 30 );

    for ( int i = 0; i < 100000; ++i ){
        
        // Calculate the mean and sigma for the LogNormal function
        double  mean = 10;
        double  sd   = 10 * 0.5; 
        
        double gaus  = mean + _random_gen->Gaussian( sd );
        
        h1->Fill(gaus);

    }   

    h1->Draw("l");
    h1->SetStats(kFALSE);    
    c1->SaveAs( "gaus.png" );

    return 0;

}
