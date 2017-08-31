# Usage
1. Cranio-facial database creation
2. Muscles database creation
3. Reconstruction workflow

## 1 - Cranio-facial database creation
For each skull:
1. Segmentation
2. Scale by 0.0035 and translate to [0.5, 0.5, 0.5]
3. Uniform remeshing with two Hausdorf distances (rough and fine)
4. Align the skulls to template_skull.mesh
5. Generate the warping shell
6. Warp the skull
7. Create the computational domain for signed distance computation
8. Compute the signed distance with mshdist
9. Fill the wrapped surface with sphere and tets
10. Morph the reference mesh onto the skull
11. Extract surface displacement scalar field
12. Generate the shell of the morphed skull and face

## 2 -  Muscles database creation
For each mandible + muscles pair:
1. Segmentation of muscles
2. Cut the mandible in two
3. Symetrize the left part of the mandible, along with its associated muscle
4. Scale to a [0.2, 0.8] box and translate to [0.5,0.5,0.5]
5. Uniform remeshing with two Hausdorf distances (rough and fine)
6. Align the half mandible and the muscle to template_mandible_masseter.mesh
7. Generate the warping shell for the half mandible
6. Warp it
7. Create the computational domain for signed distance computation
8. Compute the signed distance with mshdist
9. Fill the wrapped surface with sphere and tets
10. Morph the reference mesh onto the skull
11. Extract surface displacement scalar field

Global step:
1. Run the PCA on the half mandibles, and extract the coefficients
2. Run the PCA on the masseters, and extract the coefficients
3. Train the Neural Network from:
..* PCA coefficients of the muscles
..* PCA coefficients or measurements of the half mandible


## 3 - Reconstruction workflow
From an input skull, how to reconstruct the corresponding face?
1. Measure the distances on the skull
2. Scale by 0.0035 and translate to [0.5, 0.5, 0.5]
3. (Optionnal) Reconstruct the masseter muscles from the Neural Network, giving measurements as inputs
4. Replace the muscles in their location
5. Align the skull (and muscles) to the reference skull
6. Generate the warping shell
7. Warp the skull or the skull + muscles
8. Fill the warped surface with sphere and tets
9. Morph the reference template to the warped object
10. Link the unknown skull to each mask of the database
11. Compute average shapes
