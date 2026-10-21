# Readme for python workflow


## Step 1 - Learning the latent representation

Execute either from a terminal:
For a CVAE model with 3 latent dimensions 
`python train_traj_vae.py --nz 3 --model "cvae"`

For a VAE model with 4 latent dimensions 
`python train_traj_vae.py --nz 4 --model "vae"`

or line by line with an editor such as Visual Studio Code

You need first to put "xgv_100_14-18.csv" obtained by executing step 2 of the R workflow in the `data` folder.

The train_trak_vae.py file will call the cnn-vae architecture from the `cvae_module.py` file.

The model parameters as well as the monitoring data will be saved in the save folder

## Step2 - Analysis and extraction with python Notebook

`comparaison_monitor.ipynb` contains the code to plot the monitoring data and compare the ELBO, DKL or reconstruction error of several models.

`cvae_analysis.ipynb` exports the latent space parameters into a .csv file, generates synthetic data if needed, and propose various plots for post-learning analysis. 
The relative path of the saved data is `"../R_workflow/pydata/"`.
