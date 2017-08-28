# Reconstruction workflow
* Export the skull to a .mesh format
* Scale the skull to [0,1]
* Align the skull to the template (super4PCS then ICP)
* Warp the skull with the template ellipsoid to clean it
* Take measures on the skull
* Run the ACP on the skull
* Extract the skin
* Run the morphing

# Database creation workflow
* For each scan:
  * Run mmgs on the skull and face
  * Warp the skull 
  
# Database addition workflow
* Align the scan
