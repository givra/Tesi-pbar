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
//              colore @ 296: kOrange+7
//              num label @ 275: %8
//              media sul tempo @ 209: 9000
//              titolo asse Y + u.m @ 293: TTH4 temperature [K]
//              range asse Y @ 295: 0.95, 1.05
//              titolo canvas @ 260: TTH4 temperature vs Time GMT+00 (2.5-hour interval)
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

void DCS_TTH4(){

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
    
    TCanvas* c1 = new TCanvas("c1", "TTH4 vs Time", 1200, 800);
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


    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_TTH4_vs_time.png");

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

    // __________________________________________________________ TTH4 _____________________________________________________________________


    TFile *hfileTTH4  = new TFile("/eos/home-g/gmeinard/root_MCfiles/DCS_data/Mxc_upstream_end_(TTH4_temperature)_2023-05-19_00_00_00_2023-06-26_00_00_00.root", "READ");
    TTree *treeTTH4 = (TTree*)hfileTTH4->Get("Mxc upstream end (TTH4 temperature)");
     
    // Variables to read from DCS TTrees

    float TTH4 = 0;
    TTimeStamp* time_TTH4 = nullptr;

    // trees branches
    treeTTH4->SetBranchAddress("timestamp", &time_TTH4);
    treeTTH4->SetBranchAddress("fValue",&TTH4);

    const int timeWindow_TTH4 = 9000;  // 2.5 hours
    int nTTH4 = treeTTH4->GetEntries();
    // cout << " n entries in TTH4 " << nTTH4 << endl;

    // Get the first entry to initialize min/max
    treeTTH4->GetEntry(0);
    double minTimeTTH4 = time_TTH4->AsDouble();
    // treeTTH4->GetEntry(nTTH4-1);
    // double maxTimeTTH4 = time_TTH4->AsDouble();

    // Now we know our time range, calculate number of windows
    const int numWindowsTTH4 = int((maxTimeEV - minTimeTTH4) / timeWindow_TTH4) + 1;
    
    // Create arrays once we know the size
    float* sumTTH4 = new float[numWindowsTTH4](); 
    double* windowTime_TTH4 = new double[numWindowsTTH4]();
    vector<TString> labels_TTH4(numWindowsTTH4); // Use TString to properly manage memory
    // counters for TTH4 values in each time window
    int* TTH4PerWindow = new int[numWindowsTTH4]();
    double* meanTTH4 = new double[numWindowsTTH4]();
    
    // Initialize window times and labels (only once)
    for (int i = 0; i < numWindowsTTH4; i++) {
        windowTime_TTH4[i] = minTimeTTH4 + (i * timeWindow_TTH4) + (timeWindow_TTH4 / 2);
        labels_TTH4[i] = TTimeStamp(windowTime_TTH4[i]).AsString();                    // time_TTH4->AsString();
    }
    
    // Now collect the data in a single pass through the events
    for (int i = 0; i < nTTH4; i++) {
        treeTTH4->GetEntry(i);
        
        int windowIndex_TTH4 = int((time_TTH4->AsDouble() - minTimeTTH4) / timeWindow_TTH4);
        
        if (windowIndex_TTH4 >= 0 && windowIndex_TTH4 < numWindowsTTH4) {
            sumTTH4[windowIndex_TTH4] += TTH4;
            TTH4PerWindow[windowIndex_TTH4] ++;
            // cout << " windows index " << windowIndex_TTH4 <<" tpo4 " << sumTTH4[windowIndex_TTH4] << endl;
        }
    }
    
    // Calculate mean
    int validPointsTTH4 = 0;
    for (int i = 0; i < numWindowsTTH4; i++) {
            meanTTH4[i] = (sumTTH4[i]/(double)(TTH4PerWindow[i]));
            validPointsTTH4++;
    }


    // Create graph with the correct number of points
    TGraph* gTTH4 = new TGraph(numWindowsTTH4, windowTime_TTH4, meanTTH4);
    gTTH4->SetName("gTTH4");
    gTTH4->SetTitle("TTH4 temperature vs Time GMT+00 (2.5-hour interval)");
    
    // TCanvas* c1 = new TCanvas("c1", "TTH4 vs Time", 1200, 800);
   // c1->Divide(1, 2);
    c1->cd(1);
    
    // Set axis properties
    gTTH4->GetXaxis()->SetTimeDisplay(1);
    
    // Apply labels (much more efficient approach)
    TAxis* xaxis_TTH4 = gTTH4->GetXaxis();
    for (int i = 0; i < numWindowsTTH4; i++) {
        // Find the bin that contains this x value
        int bin_TTH4 = xaxis_TTH4->FindBin(windowTime_TTH4[i]);

        if (i % 8 == 0) { // Only show every 8th label
            xaxis_TTH4->SetBinLabel(bin_TTH4, formatTimeStamp(windowTime_TTH4[i]).Data());
        }
        else {
            xaxis_TTH4->SetBinLabel(bin_TTH4, ""); // Empty string for other bins
        }

        //xaxis->SetBinLabel(bin, labels[i].Data());
    }
    
    gTTH4->GetXaxis()->SetLabelSize(0.05);
    gTTH4->GetXaxis()->SetLabelOffset(0.01);
    gTTH4->GetXaxis()->LabelsOption("h");
    gTTH4->GetXaxis()->SetNdivisions(510);
    
    gTTH4->SetMarkerStyle(20);
    gTTH4->SetMarkerSize(0.9);
    gTTH4->GetXaxis()->SetTitle("Date");
    gTTH4->GetYaxis()->SetTitle("TTH4 temperature [K]");
    gTTH4->GetYaxis()->LabelsOption("h");
    gTTH4->GetYaxis()->SetRangeUser(0.95, 1.05);
    gTTH4->SetMarkerColor(kOrange+7);

    gPad->SetGrid(1,1);
    gTTH4->Draw("AP");

    // Clean up
    delete[] sumTTH4;
    delete[] meanTTH4;
    delete[] TTH4PerWindow;
    delete[] windowTime_TTH4;
    
    c1->Update();
    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_TTH4_vs_time.png");
    
    hfileTTH4->Close();
    
}
