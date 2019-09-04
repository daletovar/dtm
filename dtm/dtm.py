import numpy as np 
import nibabel as nib 
import sympy
from pySINDy import SINDy 
from sklearn.decomposition import PCA,FastICA
from nilearn.input_data import NiftiMasker
from collections.abc import Iterable


class DTM(object):

    def __init__(self,
                n_conditions,
                condition_names=None,
                dt=0.01,
                poly_degree=2):
        """
        n_conditions: the number of conditions - 'EV' in fsl - integer
        conditions_names: the name of each of the conditions - list,array
        dt: the size of the time step for the SINDy differentiation 
            and the later RK4 integration - float
        poly_degree: the degree of the polynomials for the SINDy theta library - integer
        """

        self.n_conditions = n_conditions
        self.dt = dt 
        self.poly_degree=2
        if condition_names is not None:
            if not isinstance(condition_names,Iterable):
                raise NotImplementedError('condition_names must be Iterable')
        self.condition_names = condition_names
    
    
    def fit(self,
            nifti,
            onset_files,
            mask=None,
            decomp='pca', 
            n_components=3):

        """
        nifti: a 4D nifti file containing the fmri time series
        onset_files: a list of paths to the onset files for each EV
        mask: an roi mask in the same space as the time series
        decomp: decomposition method to apply - pca or ica
        n_components: number of components to use in decomposition
        """

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
        
        
        decomp_model = {'pca': PCA,
                        'ica': FastICA
                        }[decomp](n_components=n_components).fit(data.T)
        self.decomp_model = decomp_model
        data = decomp_model.transform(data.T).T

        # read onset files
        onsets = np.zeros((self.n_conditions,data.shape[1]))
        min = 100
        max = 0
        for i,file in enumerate(onset_files):
            try:
                idx = (np.loadtxt(file)[0]/2 + 2).astype(np.int_p)
                onsets[i,idx] = 1
                task_min = idx.min()
                task_max = idx.max()
                if task_min < min:
                    min = task_min
                if task_max > max:
                    max = task_max
            except:
                raise IOError('cannot read {}'.format(file))


        self.n_state_vars = data.shape[0]
        data = np.vstack((data,onsets))
        data = data[:,min:max+1]

        # sindy model here
        self.sindy_model = SINDy()
        self.sindy_model.fit(data,self.dt,poly_degree=self.poly_degree)
        self.desc_to_functions()

    @property
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
        desc = self.sindy_model.descriptions
        var_dict = {}
        for n in range(self.n_state_vars + self.n_conditions):
            var_dict['u{}'.format(n)] = sympy.Symbol('u{}'.format(n))
        
        expressions = [1]
        
        for d in desc[1:]:
            # single vars first
            if d in var_dict:
                expressions.append(var_dict[d])
            elif '^' not in d:
                splits = d.split('u')[1:] # remove empty string
                expr = var_dict['u{}'.format(splits[0])]
                for s in range(1,len(splits)):
                    expr *= var_dict['u{}'.format(splits[s])]
                    expressions.append(expr)
                var_dict[d] = expr
            else:
                loc = d.find('^')
                expr = var_dict[d[:loc]]**int(d[loc+2])
                splits = d.split('}')[1:]
                if splits[0] != '':
                    expr *= var_dict[splits[0]]
                expressions.append(expr)
        equations = []
        functions = []
        coefs = self.sindy_model.coefficients
        state_vars = []
        for n in range(self.n_state_vars + self.n_conditions): 
            state_vars.append(sympy.Symbol('u{}'.format(n)))

        for i in range(coefs.shape[1]):
            idx = np.nonzero(coefs[:,i])[0]
            expr = expressions[idx[0]] * coefs[idx[0],i]
            for j in idx[1:]:
                expr += expressions[idx[0]] * coefs[j,i]
            equations.append(expr)               
            functions.append(sympy.lambdify(state_vars,expr))
        self.equations = equations 
        self.functions = functions
    
    def predict(self,nifti,onset_files):
        """
        integrates the system using RK4 to generate predictions
        given new onset files and initial conditions

        nifti: a 4D nifti file containing the fmri time series
        onset_files: a list of paths to the onset files for each EV
        """

        # somehow get data
        data = self.decomp_model.transform(NiftiMasker(self.mask).fit_transform(nifti))
        

        # read onset files
        onsets = np.zeros((self.n_conditions,data.shape[1]))
        min = 100
        max = 0
        for i,file in enumerate(onset_files):
            try:
                idx = (np.loadtxt(file)[0]/2 + 2).astype(np.int_p)
                onsets[i,idx] = 1
                task_min = idx.min()
                task_max = idx.max()
                if task_min < min:
                    min = task_min
                if task_max > max:
                    max = task_max
            except:
                raise IOError('cannot read {}'.format(file))

        data = np.vstack((data,onsets))
        data = data[:,min:max+1]

        vals = np.zeros(data.shape)
        vals[:self.n_state_vars,0] =  data[:,0] # initial conditions
        vals[self.n_state_vars:,:] =  data[self.n_state_vars:,:] # copy onsets
        for step in range(1,data.shape[1]):
            vals[:self.n_state_vars,step] = RK4(
                        self.n_state_vars,vals[:,step-1],self.functions,self.dt)
        return vals, data


    def plot_components(self):
        """
        render pca components to a nifti and visualize them
        """

def RK4(n,y,f,step=0.2):
    """
    n: number of state variables
    y: values of each of the state variables at the current time step
    f: the derivative function - different for each state variable
    step: delta t
    """
    k1 = np.empty(n,dtype=np.float)
    k2 = np.empty(n,dtype=np.float)
    k3 = np.empty(n,dtype=np.float)
    k4 = np.empty(n,dtype=np.float)
    output = np.empty(n,dtype=np.float)
    for i in range(n):
        k1[i] = step * f[i](*y)
        k2[i] = step * f[i](*y + (0.5 * k1[i]))
        k3[i] = step * f[i](*y + (0.5 * k2[i]))
        k4[i] = step * f[i](*y + k3[i])
        output[i] = y[i] + (1.0/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
    return output
        

        

