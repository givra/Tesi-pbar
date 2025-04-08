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
//              colore @ 296: kAzure+7
//              num label @ 275: %8
//              media sul tempo @ 209: 9000
//              titolo asse Y + u.m @ 293: He4 flow rate [mmol/s]
//              range asse Y @ 295: 90, 105
//              titolo canvas @ 260: He4 flow rate vs Time GMT+00 (1.5-hour interval)
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

void DCS_FM60(){

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
    
    TCanvas* c1 = new TCanvas("c1", "FM60 vs Time", 1200, 800);
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


    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_FM60_vs_time.png");

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

    // __________________________________________________________ FM60 _____________________________________________________________________


    TFile *hfileFM60  = new TFile("/eos/home-g/gmeinard/root_MCfiles/DCS_data/He4_flow_rate_(FM60)_2023-05-19_00_00_00_2023-06-26_00_00_00.root", "READ");
    TTree *treeFM60 = (TTree*)hfileFM60->Get("He4 flow rate (FM60)");
     
    // Variables to read from DCS TTrees

    float FM60 = 0;
    TTimeStamp* time_FM60 = nullptr;

    // trees branches
    treeFM60->SetBranchAddress("timestamp", &time_FM60);
    treeFM60->SetBranchAddress("fValue",&FM60);

    const int timeWindow_FM60 = 9000;  // 2.5 hours
    int nFM60 = treeFM60->GetEntries();
    // cout << " n entries in FM60 " << nFM60 << endl;

    // Get the first entry to initialize min/max
    treeFM60->GetEntry(0);
    double minTimeFM60 = time_FM60->AsDouble();
    // treeFM60->GetEntry(nFM60-1);
    // double maxTimeFM60 = time_FM60->AsDouble();

    // Now we know our time range, calculate number of windows
    const int numWindowsFM60 = int((maxTimeEV - minTimeFM60) / timeWindow_FM60) + 1;
    
    // Create arrays once we know the size
    float* sumFM60 = new float[numWindowsFM60](); 
    double* windowTime_FM60 = new double[numWindowsFM60]();
    vector<TString> labels_FM60(numWindowsFM60); // Use TString to properly manage memory
    // counters for FM60 values in each time window
    int* FM60PerWindow = new int[numWindowsFM60]();
    double* meanFM60 = new double[numWindowsFM60]();
    
    // Initialize window times and labels (only once)
    for (int i = 0; i < numWindowsFM60; i++) {
        windowTime_FM60[i] = minTimeFM60 + (i * timeWindow_FM60) + (timeWindow_FM60 / 2);
        labels_FM60[i] = TTimeStamp(windowTime_FM60[i]).AsString();                    // time_FM60->AsString();
    }
    
    // Now collect the data in a single pass through the events
    for (int i = 0; i < nFM60; i++) {
        treeFM60->GetEntry(i);
        
        int windowIndex_FM60 = int((time_FM60->AsDouble() - minTimeFM60) / timeWindow_FM60);
        
        if (windowIndex_FM60 >= 0 && windowIndex_FM60 < numWindowsFM60) {
            sumFM60[windowIndex_FM60] += FM60;
            FM60PerWindow[windowIndex_FM60] ++;
            // cout << " windows index " << windowIndex_FM60 <<" tpo4 " << sumFM60[windowIndex_FM60] << endl;
        }
    }
    
    // Calculate mean
    int validPointsFM60 = 0;
    for (int i = 0; i < numWindowsFM60; i++) {
            meanFM60[i] = (sumFM60[i]/(double)(FM60PerWindow[i]));
            validPointsFM60++;
    }


    // Create graph with the correct number of points
    TGraph* gFM60 = new TGraph(numWindowsFM60, windowTime_FM60, meanFM60);
    gFM60->SetName("gFM60");
    gFM60->SetTitle("He4 flow rate vs Time GMT+00 (2.5-hour interval)");
    
    // TCanvas* c1 = new TCanvas("c1", "FM60 vs Time", 1200, 800);
   // c1->Divide(1, 2);
    c1->cd(1);
    
    // Set axis properties
    gFM60->GetXaxis()->SetTimeDisplay(1);
    
    // Apply labels (much more efficient approach)
    TAxis* xaxis_FM60 = gFM60->GetXaxis();
    for (int i = 0; i < numWindowsFM60; i++) {
        // Find the bin that contains this x value
        int bin_FM60 = xaxis_FM60->FindBin(windowTime_FM60[i]);

        if (i % 8 == 0) { // Only show every 8th label
            xaxis_FM60->SetBinLabel(bin_FM60, formatTimeStamp(windowTime_FM60[i]).Data());
        }
        else {
            xaxis_FM60->SetBinLabel(bin_FM60, ""); // Empty string for other bins
        }

        //xaxis->SetBinLabel(bin, labels[i].Data());
    }
    
    gFM60->GetXaxis()->SetLabelSize(0.05);
    gFM60->GetXaxis()->SetLabelOffset(0.01);
    gFM60->GetXaxis()->LabelsOption("h");
    gFM60->GetXaxis()->SetNdivisions(510);
    
    gFM60->SetMarkerStyle(20);
    gFM60->SetMarkerSize(0.9);
    gFM60->GetXaxis()->SetTitle("Date");
    gFM60->GetYaxis()->SetTitle("He4 flow rate [mmol/s]");
    gFM60->GetYaxis()->LabelsOption("h");
    gFM60->GetYaxis()->SetRangeUser(90, 105);
    gFM60->SetMarkerColor(kAzure+7);

    gPad->SetGrid(1,1);
    gFM60->Draw("AP");

    // Clean up
    delete[] sumFM60;
    delete[] meanFM60;
    delete[] FM60PerWindow;
    delete[] windowTime_FM60;
    
    c1->Update();
    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_FM60_vs_time.png");
    
    hfileFM60->Close();
    
}
