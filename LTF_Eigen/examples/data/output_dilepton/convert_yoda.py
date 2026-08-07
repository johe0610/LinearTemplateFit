#! /usr/bin/env python


############################################################################################################
#                                                                                                          #
#  python script to convert yoda file to root                                                              #
#  Takes the (compressed) yoda file as input                                                               #
#  Script modified from                                                                                    #
#  https://gitlab.cern.ch/atlas/athena/-/blob/main/Generators/Rivet_i/examples/convert2root                #
#  To run the conversion do                                                                                #
#  $ python convert_yoda.py WbWb_singlelepton.yoda.gz                                                      #
#                                                                                                          #
############################################################################################################

from array import array
import ROOT as rt
import yoda, sys


fName = str(sys.argv[1])
yodaAOs = yoda.read(fName)
rtFile = rt.TFile(fName[:fName.find('.yoda.gz')] + '.root', 'recreate')
#rtFile = rt.TFile(fName[:fName.find('.yoda')] + '.root', 'recreate')
print("Convert "+fName+" to root")
for name in yodaAOs:
  print("Now processing "+name)
  yodaAO = yodaAOs[name];  rtAO = None
  if 'Histo1D' in str(yodaAO):
    rtAO = rt.TH1D(name, '', yodaAO.numBins(), array('d', yodaAO.xEdges()))
    rtAO.Sumw2(); rtErrs = rtAO.GetSumw2()
    for i in range(rtAO.GetNbinsX()):
      rtAO.SetBinContent(i + 1, yodaAO.bin(i+1).sumW())
      rtErrs.AddAt(yodaAO.bin(i+1).sumW2(), i+1)
  elif 'Scatter2D' in str(yodaAO):
    rtAO = rt.TGraphAsymmErrors(yodaAO.numPoints())
    for i in range(yodaAO.numPoints()):
      x = yodaAO.point(i).x(); y = yodaAO.point(i).y()
      xLo, xHi = yodaAO.point(i).xErrs()
      yLo, yHi = yodaAO.point(i).yErrs()
      rtAO.SetPoint(i, x, y)
      rtAO.SetPointError(i, xLo, xHi, yLo, yHi)
  else:
    continue
  rtAO.Write(str(name.replace("/RAW/WbWb_singlelepton/","")))
  print("Write out histogram " + name.replace("/RAW/WbWb_singlelepton/",""))
rtFile.Close()
