#--------------------------------------------------------------------------------------------------------------------------------------------#
#								     DESIGNN manual							     #
#--------------------------------------------------------------------------------------------------------------------------------------------#

Basic Requirements for database-preparation , training and prediction using DESIGNN (Versions) :
i.	Python 3.9
ii.	Numpy (1.20.3)
iii.    Pytorch (2.2.2)
iv. 	Pytorch_Geometric (2.5.2)
v.      Fortran 90


#--------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------DATABASE PREPARATION FOR THE DESIGNN APPROACH------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------------------------#


#---> 1. Generation of database for training the DESIGNN model

	 Prepare an input file of the structure in the following format with name structure-01


	 X	0.000000	0.000000	0.000000
	 X	0.350000	-0.94567	1.005666
 	 .
 	 .
         .
	 Mg	3.67778		6.334451	-1.990003
	 Mg	-4.55699	1.332367	 2.334567
	 .
	 .
	 .

	which contains the coordinates of the minimum first followed by the coordinates of the metal atoms.


#---> 2. Prepare the atomic description files as follows: 


	(a) RSF_Parameter_file : contains the parameters to evaluate the radial feature descriptor (Behler’s ACSF) in the format of
	     η (Gaussian width)  , Rs (Shift Parameter), Rc (cutoff-radii).

	    Example :

	    0.00500 0.00000 4.20000
	    0.01500 0.00000 4.20000
	    0.12500 0.00000 4.20000
	    ...
	    ...
	    ...

	(b) ASF_Parameter_file : contains the parameters to evaluate the angular feature descriptor (Behler’s ACSF) in the format of
	    η (Gaussian width) , Rs (Shift Parameter), ζ(angular resolution parameter), λ(cosine maximum shift parameter), Rc (cutoff-radii).

	    Example:

	    0.00015 0.00000 3.00000 -1.00000 4.20000
	    0.00025 0.00000 2.00000 -1.00000 4.20000
	    ...
	    ...
	    ...

#---> 3. Prepare a list of the all the different structure files along with its path in the following format :

	/path/to/the/structure/file/structure-01
	/path/to/the/structure/file/structure-02
	....
	....
	....


#---> 4. Compile the codes in Fortran90 as follows :

	(a) gfortran asc_order_dist.f90 feature_calculator.f90 file_reading.f90 line_counter.f90 nsf_radial.f90 pair_indices.f90
            sf_angular.f90 topo-face-association.f90 -o FilePreparation.out

	(b) gfortran transform.f90 line_counter.f90 -o transformation.out 

#---> 5. Run the Fortran executable as : ./FilePreparation.out


#---> 6. Run the programs in the [Database] folder to obtain the following files for each minimum in the structures of the database.

	(a) nAtoms.txt : contains the number of atoms of each structure in the database
	(b) atomic-features.txt: contains the feature of each atom of each structure
	(c) face-index.txt : contains the information of the polyhedral face with which the topography minima is associated.
	(d) Pairs.txt : contains all the pair informations of the database.
	(e) nPairs-list.txt : contains the number of pairs in each structure.
	(f) r-value : the projected value of r ,i.e., the  property that is to be trained.
	(g) t-value : the angle after projection θ, i.e., the property that is to be trained.



#---------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------ TRAINING THE DESIGNN MODEL------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------#


#---> 7. The [Training] folder contains the files which carries out the following functions for training:

	(a) trainDESIGNN.py : considers the inputs from the user, i.e. number of epochs, batch size, input dimension, number of
	    message passing steps and the hyperparameters, such as : optimizer, learning rate, loss function.
	(b) GraphInformation.py : reads the database for training
	(c) RNETWORK.py : trains the r-networks within the given architecture.
	(d) TNETWORK.py : train the θ-networks within the given architecture.

#---> 8. To run the training for r-value copy the database (as mentioned in 6.) along with the codes in 7. and run :
         python3.x Training/R-Property/trainDESIGNN.py

#---> 9. To run the training for t-value copy the database (as mentioned in 6.) along with the codes in 7. and run :
         python3.x Training/T-Property/trainDESIGNN.py

#---> 10. The final outcomes will yield two trained models : mpnn-DESIGNN_model-R.pth and mpnn-DESIGNN_model-T.pth 



#-------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------PREDICTION VIA THE DESIGNN APPROACH-------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------# 


#---> 11. Prepare the input file of the structure in the following format and name the file as “Full” along with a copy of the file named “Shell”:
	  Exmaple :


	  Mg 3.67778 6.334451 -1.990003
	  Mg -4.55699 1.332367 2.334567
	  .
	  .


#---> 12. Prepare the atomic description files similar to the database preparation(6.) :

	(a) RSF_Parameter_file
	(b) ASF_Parameter_file

#---> 13. Compile the codes in Fortran90 as follows :

          gfortran designn_mainroutine.f90 adjacency_matrix.f90 asc_order_dist.f90 centroidcheck.f90 faceindexfinder.f90 simplex_three.f90
          simplex_four.f90 simplex_five.f90 simplex_six.f90 feature_calculator.f90 file_reading.f90 line_counter.f90 mapping.f90 nsf_radial.f90
          pair_indices.f90 sf_angular.f90 vecdecision.f90 vectorcheck.f90 -o FinalPrediction.out

#---> 14. Copy the trained models(.pth) (from 10.) to the same path as the FinalPrediction.out executable, along with the prediction modules : 

	(a) PredictionGraphInformation.py
	(b) prediction-model-DESIGNN-R.py
	(c) prediction-model-DESIGNN-T.py
	(d) RNETWORKDESIGNN.py
	(e) TNETWORKDESIGNN.py


#---> 15. To run the prediction for r and theta properties : ./Prediction/FinalPrediction.out
          where, Prediction folder contains structure of Mg15, with the relevant files required for carrying out the prediction of the property.
          The predicted values are stored in files "predicted_r" and "predicted_t" 

#---> 16. The final structure file generated containing the predicted topography minima:
    
	   ESP_Framework_Predicted.xyz


#---------------------------------------------------------------------------------------------------------------------------------------------#
