import ROOT

root_file = "EPPI0plots.root"

cuts_text = """Active event cuts:

nPhotons >= 2, nElectrons >= 1, nProtons >= 1
channel cuts: Q2 >= 1, W > 2, y <= 0.8

Electrons:
    p >= 1
    vertex
    DC + ECAL fiducial cuts
    ePCAL > 0.07 + sampling fraction cuts
Protons:
    p >= 0.3
    vertex 
    DC + CVT fiducial cuts
    Energy loss corrections
Photons:
    FD only
    p1, p2 >= 0.4
    0.9 <= beta1, beta2 <= 1.1
    Edep1, Edep2 >= 0.15
    ECAL fiducial cuts

SKIM18: TotalCharge = 1.76953e+07 nC, EventsProcessed = 33752755, Fills = 640598"""

# Open ROOT file in UPDATE mode
f = ROOT.TFile(root_file, "UPDATE")

# Create a TPaveText object
ptext = ROOT.TPaveText(0.1, 0.1, 0.9, 0.9, "NDC")  # NDC = normalized coords
ptext.SetName("cuts_info")
ptext.SetBorderSize(1)
ptext.SetFillColor(0)
ptext.SetTextAlign(12)
ptext.SetTextFont(42)
ptext.SetTextSize(0.03)

# Add each line of the cuts text
for line in cuts_text.split("\n"):
    ptext.AddText(line)

# Write it to the ROOT file
ptext.Write()

f.Close()
print("Cuts metadata saved as a TPaveText")
