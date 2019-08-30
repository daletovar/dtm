import numpy as np 
import nibabel as nib 

from pySINDy import SINDy 
from sklearn.decomposition import PCA
from nilearn.input_data import NiftiMasker



class DTM(object):

    def __init__(self,
                n_conditions,
                dt=0.01,
                poly_degree=2):
        self.n_conditions = n_conditions
        self.dt = dt 
        self.poly_degree=2
    
    
    def fit(self,
            nifti,
            onsets,
            mask=None,
            pca=True, 
            n_components=10):

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
        
        if pca:
            pca_model = PCA(n_components=n_components).fit(data.T)
            self.pca_model = pca_model
            data = pca_model.transform(data.T).T

        # read onset files
        
        self.n_state_vars = data.shape[0]
        data = np.vstack((data,onsets))

        # sindy model here
        self.sindy_model = SINDy()
        self.sindy_model.fit(data,self.dt,poly_degree=self.poly_degree)
        


        def equations(self):
            """
            show equations
            """
            coef = self.sindy_model.coefficients
            desc = np.array(self.sindy_model.descriptions)
            for i in range(coef.shape[1] - self.n_conditions):
                eq = 'u{}/dt = '.format(i)
                idx = np.nonzero(coef[:,i])[0]
                for j in idx: 
                    eq += str(coef[j,i]) + desc[j] + ' '
                    if j != idx[-1]:
                        eq += '+ '
                print(eq)

        def desc_to_functions(self):
            """
            5.668(u0 * u1)**2
            """
        
        def RK4(n,current_vals,f,step=0.2):
            """
            n: number of state variables
            """

            t1 = t2 = t3 = np.empty(n,dtype=np.float)
            k1 = k2 = k3 = k4 = np.empty(n,dtype=np.float)
            
            for i in range(n): t1[i] = current_vals[i]+0.5*(k1[i]=step*f(x, y, i))
            for i in range(n): t1[i] = current_vals[i]+0.5*(k1[i]=step*f(x, y, i))
            for i in range(n): t1[i] = current_vals[i]+0.5*(k1[i]=step*f(x, y, i))
            for i in range(n): t1[i] = current_vals[i]+0.5*(k1[i]=step*f(x, y, i))

            for i in range(n):

            k1 = stepsize * f(current_vals)
            k2 = stepsize * f(current_vals + 0.5 * stepsize)
            k3 = stepsize * f()

        def predict(self,nifti,onsets):
            """
            integrates the system using RK4 to generate predictions
            given new onset files and initial conditions
            """
            
            def f(model,vals,i):
                exponents = model.get_ordered_poly_exponents(self.n_state_vars + self.n_conditions, 
                            self.poly_degree)[:,i]
                coef = self.sindy_model._coef[:,i]
                nnz = np.nonzero(coef)[0]
                total = 0
                for idx in nnz:




            # somehow get data
            data = self.pca_model.transform(NiftiMasker(self.mask).fit_transform(nifti))

            # read onset files

            init_conditions = data[:,0]
            vals = np.zeros((data.shape[0]+self.n_conditions, data.shape[1]))
            vals[:self.n_state_vars,0] =  data[:,0]
            vals[self.n_state_vars:,:] =  onsets

            for step in range(1,data.shape[1]):
                vals[:,step] = RK4(model,vals[:,step-1])

            return vals, data



            


        def plot_components(self):
            """
            render pca components to a nifti and visualize them
            """


        

        

