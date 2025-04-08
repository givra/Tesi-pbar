#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>    
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TText.h"
#include "TStyle.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TF1.h"
#include "TF2.h"
#include "TString.h"
#include "TPad.h"

using namespace std;

// __________________________________________
//
// DA CAMBIARE: nome tree e file @ 197-198
//              colore @ 296: kViolet+2
//              num label @ 275: %8
//              media sul tempo @ 209: 9000
//              titolo asse Y + u.m @ 293: Still pressure piezo [mbar]
//              range asse Y @ 295: 0.65, 0.75
//              titolo canvas @ 260: Still pressure piezo vs Time GMT+00 (2.5-hour interval)
//___________________________________________

// Create a custom formatter function for your timestamps
TString formatTimeStamp(double unixTime) {
    TTimeStamp ts(unixTime);
    
    // Use UInt_t instead of int
    UInt_t year, month, day, hour, min, sec;
    ts.GetDate(kFALSE, 0, &year, &month, &day);
    ts.GetTime(kFALSE, 0, &hour, &min, &sec);
    
    // Format as a shorter string
    return TString::Format("%u-%02u-%02u %02u:%02u", 
                          year, month, day, hour, min);
}

void DCS_PV8piezo(){

    TFile* eventFile = new TFile("histDCS-newdata.root", "READ");
    TTree* eventTree = (TTree*)eventFile->Get("OutTree");

    // Variables to read from 2023 eventTree
    int eventNum, nVertices, nBeam;
    double unixTime, timeInSpill;

    eventTree->SetBranchAddress("Nevent", &eventNum);
    eventTree->SetBranchAddress("UnixSeconds", &unixTime);
    eventTree->SetBranchAddress("TiS", &timeInSpill);
    eventTree->SetBranchAddress("NVrtx", &nVertices);
    eventTree->SetBranchAddress("Nbeam", &nBeam);

    const int timeWindow = 9000;  // 2.5 hours
    int nEvents = eventTree->GetEntries();

   // cout << "n entries in eventTree " << nEvents << endl;
    
    
    // Get the first entry to initialize min/max
    eventTree->GetEntry(0);
    double minTimeEV = unixTime;
    eventTree->GetEntry(nEvents-1);
    double maxTimeEV = unixTime;
    
    
    // Now we know our time range, calculate number of windows
    const int numWindowsEV = int((maxTimeEV - minTimeEV) / timeWindow) + 1;
    
    // Create arrays once we know the size
    int* sumVertices = new int[numWindowsEV]();  
    int* sumBeam = new int[numWindowsEV]();
    bool* windowHasRatioData = new bool[numWindowsEV]();
    double* meanVrtx = new double[numWindowsEV]();
    double* meanBeam = new double[numWindowsEV]();
    double* meanRatio = new double[numWindowsEV]();
    double* windowTime = new double[numWindowsEV]();
    vector<TString> labels(numWindowsEV); // Use TString to properly manage memory

    // counters for number of vrtx/beam in each time window
    int* VerticesPerWindow = new int[numWindowsEV](); 
    int* BeamPerWindow = new int[numWindowsEV](); 
    
    
    // Initialize window times and labels (only once)
    for (int i = 0; i < numWindowsEV; i++) {
        windowTime[i] = minTimeEV + (i * timeWindow) + (timeWindow / 2);
        labels[i] = TTimeStamp(windowTime[i]).AsString();
    }
    
    // Now collect the data in a single pass through the events
    for (int i = 0; i < nEvents; i++) {
        eventTree->GetEntry(i);
        
        int windowIndex = int((unixTime - minTimeEV) / timeWindow);
        
        if (windowIndex >= 0 && windowIndex < numWindowsEV && nBeam > 0) {
            sumVertices[windowIndex] += nVertices;
            if(nVertices > 0 )VerticesPerWindow[windowIndex] ++;
            sumBeam[windowIndex] += nBeam;
            if(nBeam>0) BeamPerWindow[windowIndex] ++;
            windowHasRatioData[windowIndex] = true;

        }
        
    }
    
    // Calculate ratios
    int validPoints = 0;
    for (int i = 0; i < numWindowsEV; i++) {
            meanVrtx[i] = (double)sumVertices[i]/VerticesPerWindow[i];
            meanBeam[i] = (double)sumBeam[i]/BeamPerWindow[i];
            meanRatio[i] = (meanVrtx[i] /meanBeam[i]);
            validPoints++;

            // cout << "Window " << i << ": sum vrtx " << sumVertices[i] << " num vrtx " << VerticesPerWindow[i] << " sum beam " << sumBeam[i] << " num beam " << BeamPerWindow[i] << endl;
            // cout << " mean ratio " << meanRatio[i] << endl;
    }
    
    // Create graph with the correct number of points
    TGraph* gRatio = new TGraph(numWindowsEV, windowTime, meanRatio);
    gRatio->SetName("gRatio");
    gRatio->SetTitle("Ratio nVertices/nBeam vs Time GMT+00 (2.5-hour interval)");
    
    TCanvas* c1 = new TCanvas("c1", "PV8piezo vs Time", 1200, 800);
    c1->Divide(1, 2);
    c1->cd(2);
    
    // Set axis properties
    gRatio->GetXaxis()->SetTimeDisplay(1);
    
    // Apply labels (much more efficient approach)
    TAxis* xaxis = gRatio->GetXaxis();
    for (int i = 0; i < numWindowsEV; i++) {
        // Find the bin that contains this x value
        int bin = xaxis->FindBin(windowTime[i]);

        if (i % 8 == 0) { // Only show every 8th label
            xaxis->SetBinLabel(bin, formatTimeStamp(windowTime[i]).Data());
        }
        else {
            xaxis->SetBinLabel(bin, ""); // Empty string for other bins
        }

        //xaxis->SetBinLabel(bin, labels[i].Data());
    }
    
    gRatio->GetXaxis()->SetLabelSize(0.05);
    gRatio->GetXaxis()->SetLabelOffset(0.01);
    gRatio->GetXaxis()->LabelsOption("h");
    gRatio->GetXaxis()->SetNdivisions(510);
    
    gRatio->SetMarkerStyle(20);
    gRatio->SetMarkerSize(0.9);
    gRatio->GetYaxis()->SetRangeUser(0.9985, 1.0001);
    gRatio->GetXaxis()->SetTitle("Date");
    gRatio->GetYaxis()->SetTitle("Ratio");
    gRatio->SetMarkerColor(kRed);

   // gRatio->GetYaxis()->SetRangeUser(0.95, 1);
    gPad->SetGrid(1,1);
    gRatio->Draw("AP");
    c1->Update();


    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_PV8piezo_vs_time.png");

    // Clean up
    delete[] sumVertices;
    delete[] sumBeam;
    delete[] windowHasRatioData;
    delete[] meanVrtx;
    delete[] meanBeam;
    delete[] meanRatio;
    delete[] windowTime;
    delete[] VerticesPerWindow;
    delete[] BeamPerWindow;
    eventFile->Close();

    // __________________________________________________________ PV8piezo _____________________________________________________________________


    TFile *hfilePV8piezo  = new TFile("/eos/home-g/gmeinard/root_MCfiles/DCS_data/Still_pressure_piezo_(PV8_piezo)_2023-05-19_00_00_00_2023-06-26_00_00_00.root", "READ");
    TTree *treePV8piezo = (TTree*)hfilePV8piezo->Get("Still pressure piezo (PV8 piezo)");
     
    // Variables to read from DCS TTrees

    float PV8piezo = 0;
    TTimeStamp* time_PV8piezo = nullptr;

    // trees branches
    treePV8piezo->SetBranchAddress("timestamp", &time_PV8piezo);
    treePV8piezo->SetBranchAddress("fValue",&PV8piezo);

    const int timeWindow_PV8piezo = 9000;  // 2.5 hours
    int nPV8piezo = treePV8piezo->GetEntries();
    // cout << " n entries in PV8piezo " << nPV8piezo << endl;

    // Get the first entry to initialize min/max
    treePV8piezo->GetEntry(0);
    double minTimePV8piezo = time_PV8piezo->AsDouble();
    // treePV8piezo->GetEntry(nPV8piezo-1);
    // double maxTimePV8piezo = time_PV8piezo->AsDouble();

    // Now we know our time range, calculate number of windows
    const int numWindowsPV8piezo = int((maxTimeEV - minTimePV8piezo) / timeWindow_PV8piezo) + 1;
    
    // Create arrays once we know the size
    float* sumPV8piezo = new float[numWindowsPV8piezo](); 
    double* windowTime_PV8piezo = new double[numWindowsPV8piezo]();
    vector<TString> labels_PV8piezo(numWindowsPV8piezo); // Use TString to properly manage memory
    // counters for PV8piezo values in each time window
    int* PV8piezoPerWindow = new int[numWindowsPV8piezo]();
    double* meanPV8piezo = new double[numWindowsPV8piezo]();
    
    // Initialize window times and labels (only once)
    for (int i = 0; i < numWindowsPV8piezo; i++) {
        windowTime_PV8piezo[i] = minTimePV8piezo + (i * timeWindow_PV8piezo) + (timeWindow_PV8piezo / 2);
        labels_PV8piezo[i] = TTimeStamp(windowTime_PV8piezo[i]).AsString();                    // time_PV8piezo->AsString();
    }
    
    // Now collect the data in a single pass through the events
    for (int i = 0; i < nPV8piezo; i++) {
        treePV8piezo->GetEntry(i);
        
        int windowIndex_PV8piezo = int((time_PV8piezo->AsDouble() - minTimePV8piezo) / timeWindow_PV8piezo);
        
        if (windowIndex_PV8piezo >= 0 && windowIndex_PV8piezo < numWindowsPV8piezo) {
            sumPV8piezo[windowIndex_PV8piezo] += PV8piezo;
            PV8piezoPerWindow[windowIndex_PV8piezo] ++;
            // cout << " windows index " << windowIndex_PV8piezo <<" tpo4 " << sumPV8piezo[windowIndex_PV8piezo] << endl;
        }
    }
    
    // Calculate mean
    int validPointsPV8piezo = 0;
    for (int i = 0; i < numWindowsPV8piezo; i++) {
            meanPV8piezo[i] = (sumPV8piezo[i]/(double)(PV8piezoPerWindow[i]));
            validPointsPV8piezo++;
    }


    // Create graph with the correct number of points
    TGraph* gPV8piezo = new TGraph(numWindowsPV8piezo, windowTime_PV8piezo, meanPV8piezo);
    gPV8piezo->SetName("gPV8piezo");
    gPV8piezo->SetTitle("Still pressure piezo vs Time GMT+00 (2.5-hour interval)");
    
    // TCanvas* c1 = new TCanvas("c1", "PV8piezo vs Time", 1200, 800);
   // c1->Divide(1, 2);
    c1->cd(1);
    
    // Set axis properties
    gPV8piezo->GetXaxis()->SetTimeDisplay(1);
    
    // Apply labels (much more efficient approach)
    TAxis* xaxis_PV8piezo = gPV8piezo->GetXaxis();
    for (int i = 0; i < numWindowsPV8piezo; i++) {
        // Find the bin that contains this x value
        int bin_PV8piezo = xaxis_PV8piezo->FindBin(windowTime_PV8piezo[i]);

        if (i % 8 == 0) { // Only show every 8th label
            xaxis_PV8piezo->SetBinLabel(bin_PV8piezo, formatTimeStamp(windowTime_PV8piezo[i]).Data());
        }
        else {
            xaxis_PV8piezo->SetBinLabel(bin_PV8piezo, ""); // Empty string for other bins
        }

        //xaxis->SetBinLabel(bin, labels[i].Data());
    }
    
    gPV8piezo->GetXaxis()->SetLabelSize(0.05);
    gPV8piezo->GetXaxis()->SetLabelOffset(0.01);
    gPV8piezo->GetXaxis()->LabelsOption("h");
    gPV8piezo->GetXaxis()->SetNdivisions(510);
    
    gPV8piezo->SetMarkerStyle(20);
    gPV8piezo->SetMarkerSize(0.9);
    gPV8piezo->GetXaxis()->SetTitle("Date");
    gPV8piezo->GetYaxis()->SetTitle("Still pressure piezo [mbar]");
    gPV8piezo->GetYaxis()->LabelsOption("h");
    gPV8piezo->GetYaxis()->SetRangeUser(-0.05, 0.9);
    gPV8piezo->SetMarkerColor(kTeal+3);

    gPad->SetGrid(1,1);
    gPV8piezo->Draw("AP");

    // Clean up
    delete[] sumPV8piezo;
    delete[] meanPV8piezo;
    delete[] PV8piezoPerWindow;
    delete[] windowTime_PV8piezo;
    
    c1->Update();
    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_PV8piezo_vs_time.png");
    
    hfilePV8piezo->Close();
    
}
