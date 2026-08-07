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
import math as m

fName = str(sys.argv[1])
yodaAOs = yoda.read(fName) # creates dictionary holding all the hists

rtFile = rt.TFile(fName[:fName.find('.yoda')] + '.root', 'recreate')

for name in yodaAOs:
  print("Now processing "+name+" for yodaAO "+str(yodaAOs[name]))
  yodaAO = yodaAOs[name];  # gets the histogram
  rtAO = None
  if 'Estimate1D' in str(yodaAO):
    if 'BEAMPZ' in str(yodaAO):
    #if type(yodaAO.xEdges()[0]) is str:
      continue
    rtAO = rt.TH1D(name, '', yodaAO.numBins(), array('d', yodaAO.xEdges()))
    values = yodaAO.vals()
    matrix = yodaAO.covarianceMatrix()
    rtAO.Sumw2(); rtErrs = rtAO.GetSumw2()

    for i in range(rtAO.GetNbinsX()):
      rtAO.SetBinContent(i + 1, values[i])
      rtAO.SetBinError(i + 1, m.sqrt(matrix[i][i]))
      #print("Value and error: ",values[i]," ",m.sqrt(matrix[i][i]))
  elif 'Scatter2D' in str(yodaAO):
    rtAO = rt.TGraphAsymmErrors(yodaAO.numPoints())
    for i in range(yodaAO.numPoints()):
      x = yodaAO.point(i).x(); y = yodaAO.point(i).y()
      xLo, xHi = yodaAO.point(i).xErrs()
      yLo, yHi = yodaAO.point(i).yErrs()
      rtAO.SetPoint(i, x, y)
      rtAO.SetPointError(i, xLo, xHi, yLo, yHi)
  elif 'Estimate2D' in str(yodaAO):
    if type(yodaAO.xEdges()[0]) is str:
      continue
    rtAO = rt.TH2D(name, '', yodaAO.numBinsX(), array('d', yodaAO.xEdges()), yodaAO.numBinsY(), array('d', yodaAO.yEdges()))
    values = yodaAO.vals()
    for i in range(rtAO.GetNbinsX()):
      for j in range(rtAO.GetNbinsY()):
        rtAO.SetBinContent(i + 1, j+1, values[j*rtAO.GetNbinsX()+i])
  else:
    continue
  rtAO.Write(str(name.replace("/WbWb_singlelepton/","")))
  print("Write out histogram " + name.replace("/WbWb_singlelepton/",""))
rtFile.Close()
