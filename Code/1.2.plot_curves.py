import os
import pandas as pd
import numpy as np
import pcntoolkit as pcn
from pcntoolkit.model.hbr import HBR
from matplotlib import pyplot as plt

def plot_curves_bspline(CSVname, Outdir, Feature):

    # plot growth curves
    os.chdir(Outdir)
    processing_dir = os.getcwd()

    # Feature='L_precentral_5'
    # CSVname='/public1/home/sch10024/ADHD_project/Code_All/data/CT_Final_Train_80.csv'
    # Outdir='/public1/home/sch10024/ADHD_project/Code_All/Results/CT_bspline_80Training/' + Feature

    cov_train = pd.read_csv(CSVname)
    train_site = cov_train['site'].unique()

    age_point = 131
    site_num = len(train_site)

    ################################# create dummy data
    ages = list(np.arange(6, 19.1, 0.1)) * site_num * 2 *4  # 每个中心15个年龄,58个中心，2个性别（男/女）
    x_dummy = pd.DataFrame({
        'age': ages,
    })
    sexes = [1] * (age_point * site_num)*4 + [0] * (age_point * site_num)*4  # 生成 sex 数据: 每个中心15个数据，前58个中心是男性（1），后58个中心是女性（0）
    sites = [] # 生成 site 数据

    for i in train_site:
        sites += [i] * age_point # 每个中心15个数据
    sites = (sites * 4) * 2  # 乘2以包括男性和女性

    # races = ( [1] * (age_point * site_num) + [2] * (age_point * site_num) + [3] * (age_point * site_num) + [4] * (age_point * site_num) )*2

    batch_effects_dummy = pd.DataFrame({
        'sex': sexes,
        'site': sites#,
        #'race': races
    })
    batch_effects_dummy = np.array(batch_effects_dummy)

    print(x_dummy.shape)
    print(batch_effects_dummy.shape)

    # load data
    x_train = pd.read_pickle('x_train.pkl')
    y_train = pd.read_pickle('y_train.pkl')
    cov_train = pd.read_csv(CSVname)
    roi_names = [Feature]

    # Z_ct = pd.read_pickle('Z_ct.pkl')
    # yhat_ct = pd.read_pickle('yhat_ct.pkl')
    # ys2_ct = pd.read_pickle('ys2_ct.pkl')
    # y_adj = Z_ct.to_numpy() * ys2_ct.to_numpy() - yhat_ct.to_numpy()
    # y_adj = y_adj *(-1)
    # x_test = pd.read_pickle('x_test.pkl')
    # cov_test = pd.read_csv('/public1/home/sch10024/ADHD_project/Code_All/data/CT_Final_Test_20.csv')


    # fit model for visualization
    yhat_dummy = pd.DataFrame()
    ys2_dummy = pd.DataFrame()
    for i, roi in enumerate(roi_names):
        model = pd.read_pickle('./Models/NM_0_' + str(i) + '_ct.pkl')

        # ## add by llsun: to get the scatter of training set
        # batch_effects_train = cov_train[['sex', 'site', 'Race']]
        # batch_effects_train = np.array(batch_effects_train)
        # yhat, ys2 = model.hbr.predict(X=np.array(cov_train['age']),
        #                                             batch_effects=batch_effects_train,
        #                                             batch_effects_maps=model.batch_effects_maps)
        
        yhat_dummy1, ys2_dummy1 = model.hbr.predict(X=np.array(x_dummy),
                                                    batch_effects=batch_effects_dummy,
                                                    batch_effects_maps=model.batch_effects_maps)
        yhat_dummy1 = pd.DataFrame(yhat_dummy1)
        ys2_dummy1 = pd.DataFrame(ys2_dummy1)
        yhat_dummy = pd.concat([yhat_dummy, yhat_dummy1], axis=1)
        ys2_dummy = pd.concat([ys2_dummy, ys2_dummy1], axis=1)
    del(yhat_dummy1, ys2_dummy1, model, i, roi)
    yhat_dummy.columns = roi_names
    ys2_dummy.columns = roi_names


    print("yhat_dummy:", yhat_dummy.shape)
    print("ys2_dummy:", ys2_dummy.shape)

    yhat_dummy_boy = yhat_dummy.iloc[0:age_point*site_num*4, :]
    ys2_dummy_boy = ys2_dummy.iloc[0:age_point*site_num*4, :]
    yhat_dummy_girl = yhat_dummy.iloc[age_point*site_num*4:age_point*site_num*2*4, :]
    ys2_dummy_girl = ys2_dummy.iloc[age_point*site_num*4:age_point*site_num*2*4, :]

    print("yhat_dummy_boy:", yhat_dummy_boy.shape)
    print("yhat_dummy_girl:", yhat_dummy_girl.shape)

    # average across sites
    # for boy
    yhat_dummy_mean = pd.DataFrame()
    for j in range(0, age_point):
        indices = [j + k*age_point for k in range(site_num)]
        mean = yhat_dummy_boy.iloc[indices, :].mean(axis=0)
        yhat_dummy_mean = pd.concat([yhat_dummy_mean, mean], axis=1)
    # yhat_dummy_mean.columns = ['6y', '7y', '8y', '9y', '10y', '11y', '12y', '13y', '14y', '15y', '16y', '17y', '18y', '19y', '20y']
    yhat_dummy_mean = yhat_dummy_mean.T


    ys2_dummy_mean = pd.DataFrame()
    for j in range(0, age_point):
        indices = [j + k*age_point for k in range(site_num)]
        mean = ys2_dummy_boy.iloc[indices, :].mean(axis=0)
        ys2_dummy_mean = pd.concat([ys2_dummy_mean, mean], axis=1)
    # ys2_dummy_mean.columns = ['6y', '7y', '8y', '9y', '10y', '11y', '12y', '13y', '14y', '15y', '16y', '17y', '18y', '19y', '20y']
    ys2_dummy_mean = ys2_dummy_mean.T

    yhat_dummy_mean.to_csv('yhat_mean_boy.csv', header=True, index=True)
    ys2_dummy_mean.to_csv('ys2_mean_boy.csv', header=True, index=True)

    yhat_dummy_mean = pd.DataFrame()
    for j in range(0, age_point):
        indices = [j + k*age_point for k in range(site_num)]
        mean = yhat_dummy_girl.iloc[indices, :].mean(axis=0)
        yhat_dummy_mean = pd.concat([yhat_dummy_mean, mean], axis=1)
    # yhat_dummy_mean.columns = ['6y', '7y', '8y', '9y', '10y', '11y', '12y', '13y', '14y', '15y', '16y', '17y', '18y', '19y', '20y']
    yhat_dummy_mean = yhat_dummy_mean.T
    # for girl
    ys2_dummy_mean = pd.DataFrame()
    for j in range(0, age_point):
        indices = [j + k*age_point for k in range(site_num)]
        mean = ys2_dummy_girl.iloc[indices, :].mean(axis=0)
        ys2_dummy_mean = pd.concat([ys2_dummy_mean, mean], axis=1)
    # ys2_dummy_mean.columns = ['6y', '7y', '8y', '9y', '10y', '11y', '12y', '13y', '14y', '15y', '16y', '17y', '18y', '19y', '20y']
    ys2_dummy_mean = ys2_dummy_mean.T

    yhat_dummy_mean.to_csv('yhat_mean_girl.csv', header=True, index=True)
    ys2_dummy_mean.to_csv('ys2_mean_girl.csv', header=True, index=True)
    del(j , mean)




    # select the data point belonging to male/female (1/0)
    idx_train = cov_train.index[cov_train['sex'] == 0]
    yhat_dummy_mean = pd.read_csv('yhat_mean_girl.csv', index_col=0)
    ys2_dummy_mean = pd.read_csv('ys2_mean_girl.csv', index_col=0)
    for i, roi in enumerate(roi_names):
        ax = plt.axes()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #plt.figure()
        plt.plot(list(np.arange(6, 19.1, 0.1)), yhat_dummy_mean[roi], linewidth=3, c='k')
        plt.plot(list(np.arange(6, 19.1, 0.1)), yhat_dummy_mean[roi] + 2*np.sqrt(ys2_dummy_mean[roi]), linewidth=1, linestyle='--', c='k')
        plt.plot(list(np.arange(6, 19.1, 0.1)), yhat_dummy_mean[roi] - 2*np.sqrt(ys2_dummy_mean[roi]), linewidth=1, linestyle='--', c='k')
        plt.scatter(x_train.loc[idx_train, 'age'], y_train[idx_train], c='indianred', label='Girl')  # boys: steelblue, girls:indianred
        plt.legend(loc='upper right', fontsize=16)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        #plt.title(roi)
        plt.savefig(roi + '_girl.png')
        plt.close()
    del(i, roi, ax)

    # select the data point belonging to male/female (1/0)
    idx_train = cov_train.index[cov_train['sex'] == 1]
    yhat_dummy_mean = pd.read_csv('yhat_mean_boy.csv', index_col=0)
    ys2_dummy_mean = pd.read_csv('ys2_mean_boy.csv', index_col=0)
    for i, roi in enumerate(roi_names):
        ax = plt.axes()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        #plt.figure()
        plt.plot(list(np.arange(6, 19.1, 0.1)), yhat_dummy_mean[roi], linewidth=3, c='k')
        plt.plot(list(np.arange(6, 19.1, 0.1)), yhat_dummy_mean[roi] + 2*np.sqrt(ys2_dummy_mean[roi]), linewidth=1, linestyle='--', c='k')
        plt.plot(list(np.arange(6, 19.1, 0.1)), yhat_dummy_mean[roi] - 2*np.sqrt(ys2_dummy_mean[roi]), linewidth=1, linestyle='--', c='k')
        plt.scatter(x_train.loc[idx_train, 'age'], y_train.loc[idx_train], c='steelblue', label='boy')  # boys: steelblue, girls:indianred
        plt.legend(loc='upper right', fontsize=16)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        #plt.title(roi)
        plt.savefig(roi + '_boy.png')
        plt.close()
    del(i, roi, ax)


if __name__ == '__main__':
    import sys
    CSVname = sys.argv[1]
    Outdir = sys.argv[2]
    Feature = sys.argv[3]
    
    plot_curves_bspline( CSVname, Outdir, Feature)

    