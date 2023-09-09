from subprocess import call, check_output, CalledProcessError
import subprocess
from tempfile import TemporaryDirectory
import os
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import ElasticNetCV, ElasticNet, Lasso, Ridge, BayesianRidge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import cross_validate, KFold
from sklearn.decomposition import PCA
import pandas_plink

class pyrrBLUP:
    """A class to run rrBLUP through Python
    
    Predicts phenotypes given genotype data.
    Requires companion pyrrBLUP.r script.
    Requires r enviornment.
    """
    
    def __init__(self, geno_path, pheno_path, pyrrBLUP_path, r_env_name):
        """
        Parameters
        ----------
        geno_path : str
            Path to plink genotype files (.bed, .bim, .fam)
        pheno_path : str
            Path to phenotype csv file
        pyrrBLUP_Path : str
            Path to pyrrBLUP.r script
        r_env_name : str
            Name of user's conda r enviornment 
        """
        
        self.geno = geno_path
        self.pheno = pheno_path
        self.pyrrBLUP = pyrrBLUP_path
        self.r_env = r_env_name
    
    def fit_predict(self, train_rats, test_rats, rfid_col, pheno_col): 
        """Fits rrBLUP model and returns trait predictions
        
        Fits model on training rats data and returns
        predictions for all rats in a DataFrame. 
        
        Parameters
        ----------
        train_rats : numpy.array
            Array of train animal rfids. 
            Intended input: X_train.index
        test_rats : numpy.array
            Array of test animal rfids
            Inteded input: X_test.index
        rfid_col : str
            Name of rfid column in self.pheno
        pheno_col : str
            Name of phenotype trait column in self.pheno

        Raises
        ------
        CalledProcessError
            If there is an error in pyrrBLUP.r script. 
            Most likely due to incorrect paths or column names.
        """
        
        temp_outdir = TemporaryDirectory()
        
        # set up rat data files
        with open(os.path.join(temp_outdir.name, "test_rats.txt"), 'w') as f:
            f.writelines(["%s\n" %l for l in test_rats])
        
        with open(os.path.join(temp_outdir.name, "all_rats.txt"), 'w') as f:
            f.writelines(["0 %s\n" %l for l in np.concatenate((test_rats,train_rats))])
   
        # call rrBLUP
        cmd = ['Rscript', self.pyrrBLUP]
        cmd_env = ['conda', 'run', '-n', self.r_env] 
        
        try:
            subprocess.check_call(cmd_env + cmd + [self.geno, self.pheno, temp_outdir.name, rfid_col, 
                                    pheno_col], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # read in predictions
            preds = pd.read_csv(os.path.join(temp_outdir.name, 'preds.csv'))
            preds.rename(columns={preds.columns[0]: 'rfid'}, inplace=True)
            
        except CalledProcessError:
            print("Something went wrong in getting your predictions.\n" \
                    "Check file paths and column names are accurate.")
            preds = None

        # clean up temporary directory/in-between files
        temp_outdir.cleanup()
        return preds
    
    
# df has rfid, obs, and pred columns 
def score_clf(df, test_rfids, verbose=True):
    all_rats = r2_score(df.obs, df.pred)
    test_rats = r2_score(df[df.rfid.isin(test_rfids)].obs, df[df.rfid.isin(test_rfids)].pred)
    train_rats = r2_score(df[~df.rfid.isin(test_rfids)].obs, df[~df.rfid.isin(test_rfids)].pred)
    if verbose:
        print(f'all rats r2: {all_rats}')
        print(f'test rats r2: {test_rats}')
        print(f'train rats r2: {train_rats}')
    return [all_rats, test_rats, train_rats]

def decompose_grm(grm_path, n_comp = 50, verbose = True):
    (grmxr, ar) =  pandas_plink.read_grm(grm_path)
    mass_rats = pheno.rfid.to_numpy()[np.isin(pheno.rfid.to_numpy(), grmxr['sample_0'].values)]
    grmxr = grmxr.sel(sample_0=mass_rats, sample_1=mass_rats)
    print(len(grmxr.sample_0))
    eigval, eigvec = np.linalg.eig(grmxr)
    eig_pairs = [(np.abs(eigval[i]), eigvec[:,i]) for i in range(len(eigval))]
    eig_pairs.sort(key=lambda x: x[0], reverse=True)
    explained_var = np.array([eig_pairs[i][0] for i in range(n_comp)])/sum(list(zip(*eig_pairs))[0])
    if verbose: print(f'explained_variances:{explained_var}\
          \n total explained var:{sum(explained_var)}' )
    return pd.DataFrame(np.vstack((eig_pairs[i][1] for i in range(n_comp))).T,
             columns = [f'GRM_PC{x}' for x in range(1, n_comp+1)],
             index= grmxr.sample_0.astype(str))

def decompose_grm_pca(grm_path, n_comp = 5, verbose = True):
    (grmxr, ar) =  pandas_plink.read_grm(grm_path)
    pca = PCA(n_components=n_comp)
    decomp = pca.fit_transform(grmxr)
    if verbose: print(f'explained_variances:{pca.explained_variance_ratio_}\
          \ntotal explained var:{sum(pca.explained_variance_ratio_)}' )
    return pd.DataFrame(decomp, columns = [f'GRM_PC{x}' for x in range(1, decomp.shape[1]+1)],
             index= grmxr.sample_0.astype(str))

def plink2pddf(plinkpath,rfids = 0, c = 0, pos_start = 0, pos_end = 0, snplist = 0):
    snps, iid, gen = pandas_plink.read_plink(plinkpath)
    
    # snps.chrom = snps.chrom.astype(int)
    snps.pos = snps.pos.astype(int)
    isnps = snps.set_index(['chrom', 'pos'])
    iiid = iid.set_index('iid')
    if not snplist:
        if (pos_start == pos_start == 0 ):
            if not c: index = isnps
            else: index = isnps.loc[(slice(c, c)), :]
        index = isnps if (not c and not pos_start and not pos_end) else isnps.loc[(slice(c, c),slice(pos_start, pos_end) ), :]
    else:
        index = isnps.set_index('snp').loc[snplist].reset_index()
    col = iiid  if not rfids else iiid.loc[rfids]
    return pd.DataFrame(gen.astype(np.float16)[index.i.values ][:, col.i].T, index = col.index.values.astype(str), columns = index.snp.values )

pheno = pd.read_csv('/projects/ps-palmer/bbjohnson/rattaca/nida_poster/all_pheno_merged.csv')
pheno = pheno[~pheno.mass.isna()][['rfid', 'mass']]
pheno.head()
pca = decompose_grm('/projects/ps-palmer/tsanches/gwaspipeline/gwas/allratsGRM/grm/AllchrGRM.grm.bin', 400)

file = open('../snp_sets/test_clumps.clumped', 'r')
lines = file.readlines()
file.close()
clump = pd.DataFrame([line.strip().split() for line in lines[1:]], columns = lines[0].strip().split())
clump.SP2 = clump.SP2.apply(lambda x: [snp for snp in x.split(',')])
clump.CHR = clump.CHR.apply(lambda x: 'X' if x=='23' else 'MT' if x=='26' else x)
clump[['F', 'BP', 'P', 'TOTAL', 'NSIG', 'S05', 'S01', 'S001', 'S0001']] = clump[['F', 'BP', 'P', 'TOTAL', 'NSIG', 'S05', 'S01', 'S001', 'S0001']].apply(pd.to_numeric)
snps = [i.split('(')[0] for i in sum(clump[clump['P'] < 1e-20].SP2.to_numpy(), [])]
snps = np.array([i for i in snps if i != "NONE"])
isnps = clump[clump['P'] < 1e-12].SNP.to_numpy()
allsnps = np.concatenate((isnps, snps))

# get snps that have low p values on the gwas on mass
gen_window = plink2pddf( '/projects/ps-palmer/hs_rats/round10/round10',
         snplist= isnps.tolist())
gen_window = gen_window.loc[:, :]
gen_window = gen_window.T.drop_duplicates().T
gen_window

df = pd.concat([pca, gen_window, pheno.set_index('rfid')[['mass']]], axis = 1).dropna(subset=['mass', 'GRM_PC1'])
df.columns = df.columns.str.replace(':', ".")
df.head()

X_train, X_test, y_train, y_test = train_test_split(df.drop('mass', axis=1), df.mass)

print('--rrBLUP--')
# set up rrblup object
blupclf = pyrrBLUP('/projects/ps-palmer/hs_rats/round10/round10', '/projects/ps-palmer/bbjohnson/rattaca/nida_poster/all_pheno_merged.csv', '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/scripts/pyrrBLUP/pyrrBLUP.r', 'r_env')
# fit and predict
bluppreds = blupclf.fit_predict(X_train.index, X_test.index, 'rfid', 'mass')
if type(bluppreds) == pd.core.frame.DataFrame: score_clf(bluppreds, X_test.index, verbose=True)
