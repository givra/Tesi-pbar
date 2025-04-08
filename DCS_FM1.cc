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

using namespace std;

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

void DCS_FM1(){

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
    
    TCanvas* c1 = new TCanvas("c1", "FM1 vs Time", 1200, 800);
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


    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_FM1_vs_time.png");

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

    // __________________________________________________________ FM1 _____________________________________________________________________


    TFile *hfileFM1  = new TFile("/eos/home-g/gmeinard/root_MCfiles/DCS_data/He3_flow_rate_(FM1)_2023-05-19_00_00_00_2023-06-26_00_00_00.root", "READ");
    TTree *treeFM1 = (TTree*)hfileFM1->Get("He3 flow rate (FM1)");
     
    // Variables to read from DCS TTrees

    float FM1 = 0;
    TTimeStamp* time_FM1 = nullptr;

    // trees branches
    treeFM1->SetBranchAddress("timestamp", &time_FM1);
    treeFM1->SetBranchAddress("fValue",&FM1);

    const int timeWindow_FM1 = 9000;  // 2.5 hours
    int nFM1 = treeFM1->GetEntries();
    // cout << " n entries in FM1 " << nFM1 << endl;

    // Get the first entry to initialize min/max
    treeFM1->GetEntry(0);
    double minTimeFM1 = time_FM1->AsDouble();
    // treeFM1->GetEntry(nFM1-1);
    // double maxTimeFM1 = time_FM1->AsDouble();

    // Now we know our time range, calculate number of windows
    const int numWindowsFM1 = int((maxTimeEV - minTimeFM1) / timeWindow_FM1) + 1;
    
    // Create arrays once we know the size
    float* sumFM1 = new float[numWindowsFM1](); 
    double* windowTime_FM1 = new double[numWindowsFM1]();
    vector<TString> labels_FM1(numWindowsFM1); // Use TString to properly manage memory
    // counters for FM1 values in each time window
    int* FM1PerWindow = new int[numWindowsFM1]();
    double* meanFM1 = new double[numWindowsFM1]();
    
    // Initialize window times and labels (only once)
    for (int i = 0; i < numWindowsFM1; i++) {
        windowTime_FM1[i] = minTimeFM1 + (i * timeWindow_FM1) + (timeWindow_FM1 / 2);
        labels_FM1[i] = TTimeStamp(windowTime_FM1[i]).AsString();                    // time_FM1->AsString();
    }
    
    // Now collect the data in a single pass through the events
    for (int i = 0; i < nFM1; i++) {
        treeFM1->GetEntry(i);
        
        int windowIndex_FM1 = int((time_FM1->AsDouble() - minTimeFM1) / timeWindow_FM1);
        
        if (windowIndex_FM1 >= 0 && windowIndex_FM1 < numWindowsFM1) {
            sumFM1[windowIndex_FM1] += FM1;
            FM1PerWindow[windowIndex_FM1] ++;
        }
    }
    
    // Calculate mean
    int validPointsFM1 = 0;
    for (int i = 0; i < numWindowsFM1; i++) {
            meanFM1[i] = (sumFM1[i]/(double)(FM1PerWindow[i]));
            validPointsFM1++;
    }


    // Create graph with the correct number of points
    TGraph* gFM1 = new TGraph(numWindowsFM1, windowTime_FM1, meanFM1);
    gFM1->SetName("gFM1");
    gFM1->SetTitle("He3 flow rate vs Time GMT+00 (2.5-hour interval)");
    
    // TCanvas* c1 = new TCanvas("c1", "FM1 vs Time", 1200, 800);
   // c1->Divide(1, 2);
    c1->cd(1);
    
    // Set axis properties
    gFM1->GetXaxis()->SetTimeDisplay(1);
    
    
    // Apply labels (much more efficient approach)
    TAxis* xaxis_FM1 = gFM1->GetXaxis();
    for (int i = 0; i < numWindowsFM1; i++) {
        // Find the bin that contains this x value
        int bin_FM1 = xaxis_FM1->FindBin(windowTime_FM1[i]);

        if (i % 8 == 0) { // Only show every 8th label
            xaxis_FM1->SetBinLabel(bin_FM1, formatTimeStamp(windowTime_FM1[i]).Data());
        }
        else {
            xaxis_FM1->SetBinLabel(bin_FM1, ""); // Empty string for other bins
        }

        //xaxis->SetBinLabel(bin, labels[i].Data());
    }
    
    gFM1->GetXaxis()->SetLabelSize(0.05);
    gFM1->GetXaxis()->SetLabelOffset(0.01);
    gFM1->GetXaxis()->LabelsOption("h");
    gFM1->GetXaxis()->SetNdivisions(510);
    
    gFM1->SetMarkerStyle(20);
    gFM1->SetMarkerSize(0.9);
    gFM1->GetXaxis()->SetTitle("Date");
    gFM1->GetYaxis()->SetTitle("He3 flow rate [SLPM]");
    gFM1->GetYaxis()->LabelsOption("h");
    gFM1->GetYaxis()->SetRangeUser(8, 11.5);
    gFM1->SetMarkerColor(kBlue);

    gPad->SetGrid(1,1);
    gFM1->Draw("AP");

    // Clean up
    delete[] sumFM1;
    delete[] meanFM1;
    delete[] FM1PerWindow;
    delete[] windowTime_FM1;
    
    c1->Update();
    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_FM1_vs_time.png");
    
    hfileFM1->Close();
    
}
