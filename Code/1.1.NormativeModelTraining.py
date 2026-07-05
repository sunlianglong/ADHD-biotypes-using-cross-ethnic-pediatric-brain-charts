
import pandas as pd
import os
import sys
import pickle
import pcntoolkit as pcn
from sklearn.model_selection import train_test_split

current_dir = os.path.dirname(os.path.abspath(__file__)) # Gets the directory where the current script is located
sys.path.append(current_dir) # Add the script directory to the Python module search path
from plot_curves_bspline import plot_curves_bspline

def fun_Step1_TrainingModel(TrainDataPath, TestDataPath, Feature, ResultDir):

    ############################# Step 2: Prepare features
    # TrainDataPath = '/public1/home/xxxxx/ADHD_project/TrainSet_CT_Cammoun2012_res125.csv'
    # TestDataPath = '/public1/home/xxxxx/ADHD_project/TestSet_CT_Cammoun2012_res125.csv'

    # Feature = 'L_lateralorbitofrontal_1'

    # ResultDir = '/public1/home/xxxxx/ADHD_project/Results/CT/' + Feature

    TrainSet_data = pd.read_csv(TrainDataPath)
    TestSet_data = pd.read_csv(TestDataPath)


    # select training cov
    x_train = TrainSet_data[['age']]
    # batch effect for train set
    batch_effects_train = TrainSet_data[['sex', 'site']]
    # select training imaging features
    y_train = TrainSet_data.loc[:, Feature]

    # preparing testing dataset
    x_test = TestSet_data[['age']]
    # batch effect for test set
    batch_effects_test = TestSet_data[['sex', 'site']]
    # imaging features for test set
    y_test = TestSet_data.loc[:, Feature]

    # Step 4: Create binary files for model estimation
    os.makedirs(ResultDir, exist_ok=True)
    with open(ResultDir + '/x_train.pkl', 'wb') as file:
        pickle.dump(x_train, file)
    with open(ResultDir + '/y_train.pkl', 'wb') as file:
        pickle.dump(y_train, file)
    with open(ResultDir + '/trbefile.pkl', 'wb') as file:
        pickle.dump(batch_effects_train, file)

    with open(ResultDir + '/x_test.pkl', 'wb') as file:
        pickle.dump(x_test, file)
    with open(ResultDir + '/y_test.pkl', 'wb') as file:
        pickle.dump(y_test, file)
    with open(ResultDir + '/tsbefile.pkl', 'wb') as file:
        pickle.dump(batch_effects_test, file)

    # Step 5: Set modeling choices
    output_path = os.path.join(ResultDir, 'Models')
    log_dir = os.path.join(ResultDir, 'log')
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)


    likelihood = 'Normal'
    model_type = 'bspline'
    linear_mu = 'True'
    linear_sigma = 'True'
    random_slope_mu = 'False'
    random_intercept_mu = 'True'
    random_intercept_sigma = 'False'
    inscaler_type = 'None'
    outscaler_type = 'None'

    n_cores = '1'
    alg = 'hbr'
    saveoutput = 'True'
    savemodel = 'True'
    binary = 'True'
    outputsuffix = '_ct'

    os.chdir(ResultDir)
    # Step 6: Estimate and apply the normative model to patients and hc_test
    if not os.path.exists(os.path.join(output_path, 'NM_0_0_ct.pkl')):
        pcn.normative.estimate(covfile=os.path.join(ResultDir, 'x_train.pkl'),
                            respfile=os.path.join(ResultDir, 'y_train.pkl'),
                            trbefile=os.path.join(ResultDir, 'trbefile.pkl'),
                            testcov=os.path.join(ResultDir, 'x_test.pkl'),
                            testresp=os.path.join(ResultDir, 'y_test.pkl'),
                            tsbefile=os.path.join(ResultDir, 'tsbefile.pkl'),
                            log_path=log_dir,
                            saveoutput=saveoutput,
                            output_path=output_path,
                            savemodel=savemodel,
                            binary=binary,
                            outputsuffix=outputsuffix,
                            alg=alg,
                            cores=n_cores,
                            inscaler=inscaler_type,
                            outscaler=outscaler_type,
                            likelihood=likelihood,
                            model_type=model_type,
                            linear_mu=linear_mu,
                            random_slope_mu=random_slope_mu,
                            random_intercept_mu=random_intercept_mu,
                            linear_sigma=linear_sigma,
                            random_intercept_sigma=random_intercept_sigma)

        # Step 7 Check and save results
        mod = pd.read_pickle(os.path.join(output_path, 'NM_0_0_ct.pkl'))
        Z_ct = pd.read_pickle(os.path.join(ResultDir, 'Z_ct.pkl'))
        # x_test = pd.read_pickle('x_test.pkl')
        Z_ct.columns = [Feature]
        Z_ct.index = x_test.index
        required_columns = ['id', 'sex', 'age', 'site', 'group']
        Z_ct = pd.merge(TestSet_data[required_columns], Z_ct, how='inner', left_index=True, right_index=True)
        Z_ct.to_csv(ResultDir + '/Z_ct.csv', header=True, index=True)
    else:
        print(f"File {os.path.join(output_path, 'NM_0_0_ct.pkl')} already exists. Skipping estimation and result saving steps.")

    # Step 8: Plot growth curves
    plot_curves_bspline(TrainDataPath, ResultDir, Feature)

    # Step 9: 10-fold within training set to evaluate the robustness of the models
    os.makedirs(ResultDir + '/10foldCV', exist_ok=True)
    os.chdir(ResultDir + '/10foldCV')
    # ResultDir = os.getcwd()

    pcn.normative.estimate(covfile=os.path.join(ResultDir, 'x_train.pkl'),
                           respfile=os.path.join(ResultDir, 'y_train.pkl'),
                           trbefile=os.path.join(ResultDir, 'trbefile.pkl'),
                           cvfolds=10,
                           alg='hbr',
                           log_path=os.path.join(ResultDir, 'log/'),
                           output_path=os.path.join(ResultDir, 'Models/'),
                           outputsuffix='_10fold',
                           savemodel=False,
                           saveoutput=True,
                           cores=n_cores,
                           inscaler=inscaler_type,
                           outscaler=outscaler_type,
                           likelihood=likelihood,
                           model_type=model_type,
                           linear_mu=linear_mu,
                           random_slope_mu=random_slope_mu,
                           random_intercept_mu=random_intercept_mu,
                           linear_sigma=linear_sigma,
                           random_intercept_sigma=random_intercept_sigma,
                           binary=binary)

    EV_ct = pd.read_pickle('EXPV_ct.pkl')
    MSLL_ct = pd.read_pickle('MSLL_ct.pkl')
    NLL_ct = pd.read_pickle('NLL_ct.pkl')
    pRho_ct = pd.read_pickle('pRho_ct.pkl')
    Rho_ct = pd.read_pickle('Rho_hbrct.pkl')
    RMSE_ct = pd.read_pickle('RMSE_hbrct.pkl')
    SMSE_ct = pd.read_pickle('SMSE_hbrct.pkl')

    # ---------- [ADD] save model performance to csv ----------
    def _to_scalar(x):
        """Convert pcn pickles (could be scalar/np array/pd obj) to a python scalar if possible."""
        if isinstance(x, (pd.Series, pd.DataFrame)):
            x = x.values
        try:
            import numpy as np
            if isinstance(x, np.ndarray):
                if x.size == 1:
                    return float(x.reshape(-1)[0])
                # if multiple values (e.g., per-fold), store mean
                return float(np.nanmean(x))
        except Exception:
            pass
        # python scalar
        try:
            return float(x)
        except Exception:
            return x

    perf = pd.DataFrame([{
        "EV":   _to_scalar(EV_ct),
        "MSLL": _to_scalar(MSLL_ct),
        "NLL":  _to_scalar(NLL_ct),
        "pRho": _to_scalar(pRho_ct),
        "Rho":  _to_scalar(Rho_ct),
        "RMSE": _to_scalar(RMSE_ct),
        "SMSE": _to_scalar(SMSE_ct),
    }])

    perf.to_csv("ModelPerformance.csv", index=False)
    # ---------------------------------------------------------


if __name__ == '__main__':
    import sys
    TrainDataPath = sys.argv[1]
    TestDataPath = sys.argv[2]
    Feature = sys.argv[3]
    ResultDir = sys.argv[4]
    
    fun_Step1_TrainingModel(TrainDataPath, TestDataPath, Feature, ResultDir)