import numpy as np 
import nibabel as nib 

from sklearn.decomposition import PCA
from nilearn.input_data import NiftiMasker



class DTM():

    def __init__(self,
                n_conditions,
                dt=0.01,
                poly_degree=2
                )
        self.n_conditions = n_conditions
        self.dt = dt 
        self.poly_degree=2
    
    
    def fit(self,
            nifti,
            onsets,
            mask=None,
            pca=True, 
            n_components=10,
            parcellation=None
            ):

        if not isinstance(nifti, nib.Nifti1Image):
            try:
                nifti = nib.load(nifti)
            except:
                raise ValueError('nifti must be a string or nibabel image')
        if nifti.get_data().shape != (91,109,91):
            raise NotImplementedError('nifti files must be in MNI152 space')
        if mask is None:
            mask = nib.load(MNI) # fix this
        else:
            try:
                mask = nib.load(mask)
            except:
                raise ValueError('nifti must be a string or nibabel image')
        self.mask = mask
        data = NiftiMasker(mask_img=mask).fit_transform(nifti) 

        #if parcellation is not None:
        #    try:
        #        parc = parcellation
        
        if pca:
            pca_model = PCA(n_components=n_components).fit(data.T)
            self.pca_model = pca_model
            data = pca_model.transform(data.T).T

        # sindy model here



        def predict(self,nifti,onsets):
            """
            integrates the system using RK4 to generate predictions
            given new onset files and initial conditions
            """
            data = 


        def plot_components(self):
            """
            render pca components to a nifti and visualize them
            """


        

        



    




