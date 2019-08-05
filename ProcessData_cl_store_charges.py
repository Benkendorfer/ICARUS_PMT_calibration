# explanation of the wfm format from the begining
import struct
import numpy as np
from array import array
from ROOT import TCanvas, TGraph, TH1F, TFile, TH2F
import peakutils

import sys # for command-line input

#some input
tempStr = "filelist_"
filehead = sys.argv[1] # ask user for data directory input
filelistStr = tempStr + filehead + ".txt"
rtfileoutputStr = filehead + "_result.root"
# dataDirectory = "/home/azhang/ICARUS/PMT/Data201905/" + filehead + "/"
# print("Analyzing data from " + filehead + ". Files listed in " + filelistStr)


def isLEDOn(directoryName):
    if directoryName.find("On") > -1:
        return True
    else:
        return False

# record whether LED was on or off
ledStatus = isLEDOn(filehead)

#some constants
NSamples = 12500 # number of data points in one waveform
tpt = 1.6 # time interval between sample points, in ns.
NCH = 4  # number of PMTs
QFactor = tpt/50.0*1000.0 # convert V*ns/50Ohm to charge in pC
Nsigma = 10

Nwaves = len(open(filelistStr).readlines(  ))/4 # divided by 4 because 4 PMTs in group
# print "number of waveforms in this run: ", Nwaves

def decode_wfm(filename):
    #awave1 = np.zeros(NSamples)
    awave1=[]
    #filename = "w-w-2nd-ch1-1200v-longCable-2ms_ch1_20190502101025236"
    file_wfm = filename
    #file_root = filename +".root"
    with open(file_wfm, "rb") as bfile:
        #read the first 2 bytes
        bfile.seek(0)
        pCode = struct.unpack('H',bfile.read(2))
        #read the next 8 bytes - version number - char[8]
        vNumber = bfile.read(8)
        #read the next 1 byte - number of digits in byte count
        ndbc = struct.unpack('c',bfile.read(1))
        #read the next 4 bytes - number of bytes to the end of the file
        nbEOF = struct.unpack('i', bfile.read(4))
        #read the next 1 byte - number of bytes per point
        nbPt = struct.unpack('c',bfile.read(1))
        #print "number of bytes per point", nbPt
        #read the next 4 bytes - byte offset to begining of curve buffer
        bOffset1 = struct.unpack('i',bfile.read(4))
        #print "byte offset to begining of curve buffer ", bOffset1[0]
        #read the next 4 bytes - horizontal zoom scale factor 
        horScalF = struct.unpack('i',bfile.read(4))
        #print "horizontal scale zoom factor ", horScalF[0]
        #read the next 4 bytes - horizontal zoom position
        horZoomP = struct.unpack('f', bfile.read(4))
        #print "horizontal zoom position ", horZoomP[0]
        #read the next 8 bytes - vertical zoom scale factor
        verScaF = struct.unpack('d', bfile.read(8))
        #print "vertical scale zoom factor ", verScaF[0]
        #read the next 4 bytes - vertical zoom position
        verZoomP = struct.unpack('f', bfile.read(4))
        #print "vertical zoom position", verZoomP[0]
        #read the next 32 bytes - waveform label - char[32]
        wfmLabel = bfile.read(32)
        #print "Waveform label (user defined or ref waveform)", wfmLabel
        #read the next 4 bytes - number of FastFrams-1
        nbFF = struct.unpack('i', bfile.read(4))
        #print "Number of FastFrams-1", nbFF[0]
        #read the next 2 bytes - size of the waveform header
        sizeWfmHead = struct.unpack('H', bfile.read(2))
        #print "Size of the waveform header", sizeWfmHead[0]

        #print ""
    
        #print "waveform header:"
        #Reference file data
        typeWfmSet = struct.unpack('i',bfile.read(4))
        #print "Type of waveform set (0:single waveform, 1:FastFrame)", typeWfmSet[0]
        WfmCnt = struct.unpack('i',bfile.read(4))
        #print "Number of waveforms in the set", WfmCnt[0]
        AcqCnt = struct.unpack('l',bfile.read(8))
        #print "Acquisition counter", AcqCnt[0]
        TransCnt = struct.unpack('l',bfile.read(8))
        #print "Transaction counter", TransCnt[0]
        SlotID = struct.unpack('i',bfile.read(4)) #Not for use
        IsStaticFlag = struct.unpack('i',bfile.read(4))
        #print "IsStaticFlag, static or live 0-reference; 1-waveform or math", IsStaticFlag[0]
        WfmUpdateSpecCnt = struct.unpack('i',bfile.read(4))
        #print "WfmUpdateSpecCnt", WfmUpdateSpecCnt[0]
        ImpDimRefCnt = struct.unpack('i',bfile.read(4))
        #print "ImpDimRefCnt",ImpDimRefCnt[0]
        ExpDimRefCnt = struct.unpack('i',bfile.read(4))
        #print "ExpDimRefCnt",ExpDimRefCnt[0]
        DataType = struct.unpack('i',bfile.read(4))
        #print "DataType",DataType[0]
        GenPurposeCnt = struct.unpack('l',bfile.read(8))
        #print "GenPurposeCnt",GenPurposeCnt[0]
        AccuWfmCnt = struct.unpack('i',bfile.read(4))
        #print "AccuWfmCnt",AccuWfmCnt[0]
        TarAccuCnt = struct.unpack('i',bfile.read(4))
        #print "TarAccuCnt",TarAccuCnt[0]
        CurveCnt = struct.unpack('i',bfile.read(4))
        #print "CurveCnt",CurveCnt[0]
        nbReqFastFrame = struct.unpack('i',bfile.read(4))
        #print "nbReqFastFrame",nbReqFastFrame[0]
        nbAcqFastFrame = struct.unpack('i',bfile.read(4))
        #print "Number of frames", nbAcqFastFrame[0]
        SummFrame = struct.unpack('h',bfile.read(2))
        #print "SummFrame",SummFrame[0]
        PixMapDisFormat = struct.unpack('i',bfile.read(4))
        #print "PixMapDisFormat",PixMapDisFormat[0]
        PixMapMax = struct.unpack('l',bfile.read(8))
        #print "PixMapMax",PixMapMax[0]

        DimScale = struct.unpack('d',bfile.read(8))
        #print "DimScale",DimScale[0]
        DimOffset = struct.unpack('d',bfile.read(8))
        #print "DimOffset",DimOffset[0]
        DimSize = struct.unpack('i',bfile.read(4))
        #print "DimSize",DimSize[0]
        ExtDimUnits = bfile.read(20)
        #print "ExtDimUnits",ExtDimUnits[0]
        DimExtentMin = struct.unpack('d',bfile.read(8))
        #print "DimExtentMin",DimExtentMin[0]
        DimExtentMax = struct.unpack('d',bfile.read(8))
        #print "DimExtentMax",DimExtentMax[0]
        DimResolution = struct.unpack('d',bfile.read(8))
        #print "DimResolution",DimResolution[0]
        DimRefPoint = struct.unpack('d',bfile.read(8))
        #print "DimRefPoint",DimRefPoint[0]
        ExpDimFormat = struct.unpack('i',bfile.read(4))
        #print "ExpDimFormat",ExpDimFormat[0]
        ExpDimStorageType = struct.unpack('i',bfile.read(4))
        #print "ExpDimStorageType",ExpDimStorageType[0]
        ExpDimNvalue = struct.unpack('i',bfile.read(4))
        #print "ExpDimNvalue",ExpDimNvalue[0]
        ExpDimOverRange = struct.unpack('i',bfile.read(4))
        #print "ExpDimOverRange",ExpDimOverRange[0]
        ExpDimUnderRange = struct.unpack('i',bfile.read(4))
        #print "ExpDimUnderRange",ExpDimUnderRange[0]
        ExpDimHighRange = struct.unpack('i',bfile.read(4))
        #print "ExpDimHighRange",ExpDimHighRange[0]
        ExpDimLowRange = struct.unpack('i',bfile.read(4))
        #print "ExpDimLowRange",ExpDimLowRange[0]
        ExpDimUserScale = struct.unpack('d',bfile.read(8))
        #print "ExpDimUserScale",ExpDimUserScale[0]
        ExpDimUserUnits = bfile.read(20)
        #print "ExpDimUserUnits",ExpDimUserUnits[0]
        ExpDimUserOffset = struct.unpack('d',bfile.read(8))
        #print "ExpDimUserOffset",ExpDimUserOffset[0]
        ExpDimPointDensity = struct.unpack('d',bfile.read(8))
        #print "ExpDimPointDensity",ExpDimPointDensity[0]
        ExpDimHRef = struct.unpack('d',bfile.read(8))
        #print "ExpDimHRef",ExpDimHRef[0]
        ExpDimTrigDelay = struct.unpack('d',bfile.read(8))
        #print "ExpDimTrigDelay",ExpDimTrigDelay[0]
        
        ImpDimScale = struct.unpack('d',bfile.read(8))
        #print "ImpDimScale",ImpDimScale[0]
        ImpDimOffset = struct.unpack('d',bfile.read(8))
        #print "ImpDimOffset",ImpDimOffset[0]
        ImpDimSize = struct.unpack('i',bfile.read(4))
        #print "ImpDimSize",ImpDimSize[0]
        ImpDimUnits = bfile.read(20)
        #print "ImpDimUnits",ImpDimUnits[0]
        ImpDimExtentMin = struct.unpack('d',bfile.read(8))
        #print "ImpDimExtentMin",ImpDimExtentMin[0]
        ImpDimExtentMax = struct.unpack('d',bfile.read(8))
        #print "ImpDimExtentMax",ImpDimExtentMax[0]
        ImpDimRes = struct.unpack('d',bfile.read(8))
        #print "ImpDimRes",ImpDimRes[0]
        ImpDimRefPoint = struct.unpack('d',bfile.read(8))
        #print "ImpDimRefPoint",ImpDimRefPoint[0]
        ImpDimSpacing = struct.unpack('i',bfile.read(4))
        #print "ImpDimSpacing",ImpDimSpacing[0]
        ImpDimUserScale = struct.unpack('d',bfile.read(8))
        #print "ImpDimUserScale",ImpDimUserScale[0]
        ImpDimUserScale1 = struct.unpack('d',bfile.read(8))
        #print "ImpDimUserScale1 (caution)",ImpDimUserScale1[0]
        ImpDimUserUnits = bfile.read(20)
        #print "ImpDimUserUnits",ImpDimUserUnits[0]
        ImpDimUserOffset = struct.unpack('d',bfile.read(8))
        #print "ImpDimUserOffset",ImpDimUserOffset[0]
        ImpDimPointDensity = struct.unpack('d',bfile.read(8))
        #print "ImpDimPointDensity",ImpDimPointDensity[0]
        ImpDimPointDensity1 = struct.unpack('d',bfile.read(8))
        #print "ImpDimPointDensity1 (caution)",ImpDimPointDensity1[0]
        ImpDimHRef = struct.unpack('d',bfile.read(8))
        #print "ImpDimHRef",ImpDimHRef[0]
        ImpDimTrigDelay = struct.unpack('d',bfile.read(8))
        #print "ImpDimTrigDelay",ImpDimTrigDelay[0]

        TimeBase_RealPointSpacing = struct.unpack('i',bfile.read(4))
        #print "TimeBase_RealPointSpacing",TimeBase_RealPointSpacing[0]
        TimeBase_Sweep = struct.unpack('i',bfile.read(4))
        #print "TimeBase_Sweep",TimeBase_Sweep[0]
        TimeBase_TypeBase = struct.unpack('I',bfile.read(4))
        #print "TimeBase_TypeBase",TimeBase_TypeBase[0]

        WfmUpdaSpec_RealPointOffset = struct.unpack('i',bfile.read(4))
        #print "WfmUpdaSpec_RealPointOffset",WfmUpdaSpec_RealPointOffset[0]
        WfmUpdaSpec_TTOffset = struct.unpack('d',bfile.read(8))
        #print "WfmUpdaSpec_TTOffset",WfmUpdaSpec_TTOffset[0]
        WfmUpdaSpec_FracSec = struct.unpack('d',bfile.read(8))
        #print "WfmUpdaSpec_FracSec",WfmUpdaSpec_FracSec[0]
        WfmUpdaSpec_GmtSec = struct.unpack('i',bfile.read(4))
        #print "WfmUpdaSpec_GmtSec",WfmUpdaSpec_GmtSec[0]
    
        WfmCurveInfo_StateFlag = struct.unpack('i',bfile.read(4))
        #print "WfmCurveInfo_StateFlag",WfmCurveInfo_StateFlag[0]
        WfmCurveInfo_TypeCheckSum = struct.unpack('i',bfile.read(4))
        #print "WfmCurveInfo_TypeCheckSum",WfmCurveInfo_TypeCheckSum[0]    
        WfmCurveInfo_CheckSum = struct.unpack('H',bfile.read(2))
        #print "WfmCurveInfo_CheckSum",WfmCurveInfo_CheckSum[0]
        WfmCurveInfo_preQstartOffset = struct.unpack('i',bfile.read(4))    
        #print "WfmCurveInfo_preQstartOffset",WfmCurveInfo_preQstartOffset[0]
        WfmCurveInfo_dataStartOffset = struct.unpack('i',bfile.read(4))
        #print "WfmCurveInfo_dataStartOffset",WfmCurveInfo_dataStartOffset[0]
        WfmCurveInfo_postQstartOffset = struct.unpack('i',bfile.read(4))
        #print "WfmCurveInfo_postQstartOffset",WfmCurveInfo_postQstartOffset[0]
        WfmCurveInfo_postQStopOffset = struct.unpack('i',bfile.read(4))
        #print "WfmCurveInfo_postQStopOffset",WfmCurveInfo_postQStopOffset[0]
        WfmCurveInfo_EOCoffset = struct.unpack('i',bfile.read(4))
        #print "WfmCurveInfo_EOCoffset",WfmCurveInfo_EOCoffset[0]
        bfile.seek(bOffset1[0])
        temp = bfile.read(74) #74
        for i in range(NSamples):
            v = struct.unpack('h',bfile.read(2))
            #awave1[i] = v[0]*DimScale[0]+DimOffset[0]
            awave1.append(v[0]*DimScale[0]+DimOffset[0])
            #print awave[i]
            #gWave.SetBinContent(i+1,v[0]*DimScale[0]+DimOffset[0])
        #WfmFileCheckSum - 8 bytes
        WfmFileCheckSum = struct.unpack('l',bfile.read(8))
        #print "Wfm file check sum",WfmFileCheckSum[0]
        return awave1
        
# main function, process the file the 2nd time
def main():

    f = open(filelistStr,"r")
    rtfileoutput = TFile(rtfileoutputStr,"recreate")
    waveDir = rtfileoutput.mkdir("Waveforms")
    resultsDir = rtfileoutput.mkdir("Results")

    # prepare for histgorams
    hAmplitude_list = []
    hAmplitudeBin_list = []
    hWave_list = []
    hPedMean_list = []
    hPedWidth_list = []
    hFinalCharge_list = []
    hNbOfPulses_list = []
    hPulseStartTime_list = [] # start time of the pulse
    hPulseTimeDist_list = [] # time distribution of the pulses that are not due to fiber trigger (ie., dark pulses or else)
    hPulseWidth_list = [] # time distribution of the pulses that are not due to fiber trigger (ie., dark pulses or else)
    hPulseAmplitudeVsTime_list = [] # 2D histogram: pulse amplitude vs. pulse time bin
    hPulseAmplitudeVsWidth_list = [] # 2D histogram: pulse amplitude vs. pulse time bin
    hWaveAvg_list = [] # average of raw waveforms
    #special test on time difference between the first channel
    hTimeDiff = TH1F("hTimeDiff","Time diference",100,-10,10) # unit in ns
    hTimeDiff.SetXTitle("Time difference (ns)");
    hTimeDiff.SetYTitle("N");
    for i in range(0,NCH,1):
        # will keep a few waveforms 
        name = "Wave_" + str(i)
        hist = TH1F(name,"",NSamples,0,NSamples)
        hist.SetXTitle("Sample number")
        hist.SetYTitle("Amplitude (mV)")
        hist.SetLineColor(i+1)
        hWave_list.append(hist)
        # will keep average waveforms 
        name = "WaveAvg_" + str(i)
        hist = TH1F(name,"",NSamples,0,NSamples)
        hist.SetXTitle("Sample number")
        hist.SetYTitle("Amplitude (mV)")
        hist.SetLineColor(i+1)
        hWaveAvg_list.append(hist)
        # histogram of pulse baseline 
        name = "PedMean_" + str(i)
        hist = TH1F(name,"",1000,-5,5)
        hist.SetXTitle("Amplitude (mv)")
        hist.SetYTitle("Counts")
        hist.SetLineColor(i+1)
        hPedMean_list.append(hist)
        # histogram of pulse baseline width
        name = "PedWidth_" + str(i)
        hist = TH1F(name,"",1000,-5,5)
        hist.SetXTitle("Amplitude (mV)")
        hist.SetYTitle("Counts")
        hist.SetLineColor(i+1)
        hPedWidth_list.append(hist)
        # pulse charge due to fiber triggers, unit converted to fC
        name = "FinalCharge_" + str(i)
        # would be best to have different settings for LED on vs LED off
        if ledStatus == False:
            hist = TH1F(name,"",1000,-5,20)
        else:
            hist = TH1F(name,"",2000,-5,200)
        hist.SetXTitle("Charge (pC)")
        hist.SetYTitle("Counts")
        hist.SetLineColor(i+1)
        hFinalCharge_list.append(hist)
        # pulse amplitude, or nimimum (because of negative pulse)
        name = "PulseAmplitude_" + str(i)
        hist = TH1F(name,"",500,0,500)
        hist.SetXTitle("Amplitude (mV)")
        hist.SetYTitle("Counts")
        hist.SetLineColor(i+1)
        hAmplitude_list.append(hist)
        # pulse width, or nimimum (because of negative pulse)
        name = "PulseWidth_" + str(i)
        hist = TH1F(name,"",100,0,100)
        hist.SetXTitle("Pulse width (ns)")
        hist.SetYTitle("Counts")
        #hist.SetLineColor(i+1)
        hPulseWidth_list.append(hist)
        # time bin where pulse amplitude is found, used to determine charge integration region
        name = "PulseAmplitudeBin_" + str(i)
        hist = TH1F(name,"",200,100,300)
        hist.SetXTitle("Amplitude bin number")
        hist.SetYTitle("Counts")
        hist.SetLineColor(i+1)
        hAmplitudeBin_list.append(hist)
        # number of pulses found in the waveform not in the fiber trigger region
        name = "NbOfPulses_" + str(i)
        hist = TH1F(name,"",20,0,20)
        hist.SetXTitle("Number of pulses in a 20 #mus window")
        hist.SetYTitle("Counts")
        hist.SetLineColor(i+1)
        hNbOfPulses_list.append(hist)
        # pulse start time, only choose those pulses above threshold
        name = "PulseStartTime_" + str(i)
        hist = TH1F(name,"",200,100,300)
        hist.SetXTitle("Pulse start time bin (1.6 ns/bin)")
        hist.SetYTitle("Counts")
        #hist.SetLineColor(i+1)
        hPulseStartTime_list.append(hist)
        # pulse time distributions for those pulses are 100 ns later than the fiber triggered pulse
        name = "PulseTimeDist_" + str(i)
        hist = TH1F(name,"",NSamples/10,0,NSamples)
        hist.SetXTitle("Pulse time bin (16 ns/bin)")
        hist.SetYTitle("Counts")
        hist.SetLineColor(i+1)
        hPulseTimeDist_list.append(hist)
        # pulse time distributions for those pulses are 100 ns later than the fiber triggered pulse
        name = "PulseAmpVsTimeBin_" + str(i)
        hist = TH2F(name,"",NSamples/10,0,NSamples,500,0,500)
        hist.SetXTitle("Pulse time bin")
        hist.SetYTitle("Pulse Amplitude (mV)")
        hist.SetLineColor(i+1)
        hPulseAmplitudeVsTime_list.append(hist)
        # pulse time distributions for those pulses are 100 ns later than the fiber triggered pulse
        name = "PulseAmpVsWidth_" + str(i)
        hist = TH2F(name,"",100,0,100,500,0,500)
        hist.SetXTitle("Pulse time bin")
        hist.SetYTitle("Pulse Amplitude (mV)")
        hPulseAmplitudeVsWidth_list.append(hist)

    #process the waveforms
    baseline_mean = 0.0
    baseline_width = 0.0
    threshold = 0.0
    sumwave = [[0 for x in range(NSamples)] for y in range(NCH)]
    #special test for time diff
    flag1 = [ False for x in range(Nwaves)]
    time1 = [ 0 for x in range(Nwaves)]
    flag2 = [ False for x in range(Nwaves)]
    time2 = [ 0 for x in range(Nwaves)]
    for ch in range(NCH):
        for waveNb in range(Nwaves):
            afilename = f.readline().rstrip()
            #print afilename
            #if waveNb>20:
            #    continue
            awave = np.asarray(decode_wfm(afilename))

            baseline_mean = np.average(awave[NSamples-1000:NSamples])
            baseline_width = np.std(awave[NSamples-1000:NSamples])
            threshold = baseline_mean - Nsigma*baseline_width
            hPedMean_list[ch].Fill(baseline_mean)
            hPedWidth_list[ch].Fill(baseline_width)
            TimeBinOfAmplitude = np.argmin(awave[70:800]) + 70

            if ch == 0 and awave[TimeBinOfAmplitude]<threshold:
                flag1[waveNb] = True
                time1[waveNb] = TimeBinOfAmplitude
            if ch == 1 and awave[TimeBinOfAmplitude]<threshold:
                flag2[waveNb] = True
                time2[waveNb] = TimeBinOfAmplitude

            # get charge integration about the amplitude: 20.8 ns before to 36.8 ns after the amplitude
            #sumcharge = 0.0
            #if TimeBinOfAmplitude>=168 and TimeBinOfAmplitude<=180:
            sumcharge = np.sum(awave[TimeBinOfAmplitude-13:TimeBinOfAmplitude+23])
            pulsestartbin = TimeBinOfAmplitude
            while awave[pulsestartbin]<threshold and pulsestartbin>=TimeBinOfAmplitude-20:
                pulsestartbin -= 1
            hPulseStartTime_list[ch].Fill(pulsestartbin-1)                
            #else:
            #    sumcharge = np.sum(awave[172-13:172+23])
            fC = (sumcharge-baseline_mean*(23+13+1))*QFactor # convert charge involving the 50 Ohm impedance
            #if awave[TimeBinOfAmplitude] < threshold:
            hAmplitude_list[ch].Fill(-1000.0*(awave[TimeBinOfAmplitude]-baseline_mean))
            hAmplitudeBin_list[ch].Fill(TimeBinOfAmplitude)
            #print TimeBinOfAmplitude, awave[TimeBinOfAmplitude]
            hFinalCharge_list[ch].Fill(-1.0*fC)
            if ch == 0:
                print -1.0*fC

            # search for number of pulses 
            peakindex = peakutils.peak.indexes(-1.0*awave[0:NSamples],thres=-1.0*threshold, min_dist=31, thres_abs=True)
            for pulse_i in range(len(peakindex)):

                # skip pulse if it is below threshold in magnitude
                if awave[peakindex[pulse_i]] > threshold:
                    continue

                hPulseTimeDist_list[ch].Fill(peakindex[pulse_i]) #+TimeBinOfAmplitude+200
                hPulseAmplitudeVsTime_list[ch].Fill(peakindex[pulse_i], -1000.0*(awave[peakindex[pulse_i]]-baseline_mean))
                #determine the pulse width
                thisindex = peakindex[pulse_i]
                thiswidth = 0
                while awave[thisindex]<=threshold and thisindex>=2:
                    #print pulse_i, ch, thiswaveform[ch][thisindex], threshold[ch]
                    thiswidth += 1
                    thisindex -= 1
                    #if thisindex<=1:
                    #    print eventNb, ch, pulse_i, thiswaveform[ch][thisindex]
                thisindex = peakindex[pulse_i]
                while awave[thisindex+1]<=threshold and thisindex<=NSamples-3:
                    thiswidth += 1
                    thisindex += 1
                    #    if thisindex>=NSamples-2:
                    #        print eventNb, ch, pulse_i, thiswaveform[ch][thisindex]
                hPulseWidth_list[ch].Fill((thiswidth+2)*2)
                hPulseAmplitudeVsWidth_list[ch].Fill((thiswidth+2)*2, -1000.0*(awave[peakindex[pulse_i]]-baseline_mean))

            # if waveNb%100 ==0:
            #     # print ch, waveNb, afilename
            #     for bin in range(0,NSamples,1):
            #         hWave_list[ch].SetBinContent(bin+1,-1000.0*(awave[bin]-baseline_mean))
            #         #hWave_list[ch].SetBinContent(bin+1,awave[bin])
            #     waveDir.cd()
            #     lenname = int(len(afilename))
            #     newName = "Wave_Ch_" + str(ch) + "_" + str(waveNb)+afilename[lenname-4-17:lenname-4]
            #     titleName = "Evt "+str(waveNb)+"_Ch"+str(ch)+"_PulseBins_"
            #     for pulse_i in range(len(peakindex)):
            #         titleName += str(peakindex[pulse_i])
            #         titleName += "_"
            #     hWave_list[ch].SetName(newName)
            #     hWave_list[ch].SetTitle(titleName)
            #     hWave_list[ch].SetEntries(NSamples)
            #     hWave_list[ch].Write()

            hNbOfPulses_list[ch].Fill( len(peakindex) )
            #average of raw waveforms for each PMT
            for bin in range(0,NSamples,1):
                sumwave[ch][bin] += 1000.0*(awave[bin]-baseline_mean)
    for ch in range(NCH):
        for bin in range(NSamples):
            hWaveAvg_list[ch].SetBinContent(bin+1,sumwave[ch][bin]/Nwaves)
        waveDir.cd()
        hWaveAvg_list[ch].Write()
    for waveNb in range(Nwaves):
        if flag1[waveNb] is True and flag2[waveNb] is True:
            hTimeDiff.Fill(time1[waveNb]-time2[waveNb])
    f.close()
    resultsDir.cd()
    for i in range(0,NCH,1):
        # fit the charge distributions with Poisson distribution plus exponential background
        hFinalCharge_list[i].Write()
    for i in range(0,NCH,1):
        hAmplitude_list[i].Write()
    for i in range(0,NCH,1):
        hAmplitudeBin_list[i].Write()
    for i in range(0,NCH,1):
        hPedMean_list[i].Write()
    for i in range(0,NCH,1):
        hPulseStartTime_list[i].Write()
    for i in range(0,NCH,1):
        hPedWidth_list[i].Write()
    for i in range(0,NCH,1):
        hPulseTimeDist_list[i].Write()
        #print "number of pulses: ", hPulseTimeDist_list[i].GetEntries()
    for i in range(0,NCH,1):
        hPulseWidth_list[i].Write()
    for i in range(0,NCH,1):
        hPulseAmplitudeVsTime_list[i].Write()
    for i in range(0,NCH,1):
        hPulseAmplitudeVsWidth_list[i].Write()
    for i in range(0,NCH,1):
        #darkrate = 0.0
        #for j in range(2,21):
        #    if hNbOfPulses_list[i].GetBinContent(j)>0:
        #        darkrate += hNbOfPulses_list[i].GetBinContent(j)*(j-1)
        #darkrate = darkrate/(10.0*1.0e-6*hNbOfPulses_list[i].GetEntries())
        #title = "dark rate on average " + str(darkrate) + " Hz"
        #hNbOfPulses_list[i].SetTitle(title)
        hNbOfPulses_list[i].Write()
    hTimeDiff.Write()
    rtfileoutput.Close()

if __name__ == "__main__":
    main()
