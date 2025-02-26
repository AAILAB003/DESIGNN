#--------------------------------------------------------------------------------------------------------------------------------------------#
#								     DESIGNN manual							     #
#--------------------------------------------------------------------------------------------------------------------------------------------#

Basic Requirements for database-preparation , training and prediction using DESIGNN (Versions) :
i.	Python 3.9
ii.	Numpy (1.20.3)
iii.    Pytorch (2.2.2)
iv. 	Pytorch_Geometric (2.5.2)
v.      Fortran 90


#---------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------ TRAINING THE DESIGNN MODEL------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------#


#---> 1. The [Training] folder contains the files which carries out the following functions for training:

	(a) trainDESIGNN.py : considers the inputs from the user, i.e. number of epochs, batch size, input dimension, number of
	    message passing steps and the hyperparameters, such as : optimizer, learning rate, loss function.
	(b) GraphInformation.py : reads the database for training
	(c) RNETWORK.py : trains the r-networks within the given architecture.
	(d) TNETWORK.py : train the θ-networks within the given architecture.

#---> 2. To run the training copy the database files along with the codes in 1. and run :
         python3.x Training/R-Property/trainDESIGNN.py

#---> 3. To run the training copy the database files along with the codes in 1. and run :
         python3.x Training/T-Property/trainDESIGNN.py

#---> 4. The final outcomes will yield two trained models : mpnn-DESIGNN_model-R.pth and mpnn-DESIGNN_model-T.pth 



#-------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------PREDICTION VIA THE DESIGNN APPROACH-------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------# 


#---> 5. Copy the trained models(.pth) (from 4.) along with the prediction modules : 

	(a) PredictionGraphInformation.py
	(b) prediction-model-DESIGNN-R.py
	(c) prediction-model-DESIGNN-T.py
	(d) RNETWORKDESIGNN.py
	(e) TNETWORKDESIGNN.py


#---> 6. To run the prediction for r and theta properties : 
	 python3.x Training/R-Property/prediction-model-DESIGNN-R.py
	 python3.x Training/R-Property/prediction-model-DESIGNN-T.py

          where, Prediction folder contains structure of Mg15, with the relevant files required for carrying out the prediction of the property.
          The predicted values can be stored in files "predicted_r" and "predicted_t" 




#---------------------------------------------------------------------------------------------------------------------------------------------#
