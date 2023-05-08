"""small example in a container"""
print("Importing libraries")
import sys
sys.path.append("./CPlantBox/")
import plantbox as pb
print("libraries loaded")

rootsystem = pb.RootSystem()
# Open plant and root parameter from a file
path = "./CPlantBox/modelparameter/rootsystem/"
output = "12_rootsystem"
name = "FP"

print("Read parameter xml file")
rootsystem.readParameters(path + name + ".xml")

# Create and set geometry
# Creates a soil container
print("Create soil container")
rhizotron = pb.SDF_PlantBox(900, 900, 900)

# Pick 1, or 2
print("Set geometry")
rootsystem.setGeometry(rhizotron)  # soilcore, or rhizotron

# Initialize
print("Initialize")
rootsystem.initialize()

# Simulate
print("Run -- CPlantBox -- ")
rootsystem.simulate(12)  # days
print("analyse")
ana = pb.SegmentAnalyser(rootsystem)

print("write rootsytem")
ana.write("{}.txt".format(str(output)))
  
