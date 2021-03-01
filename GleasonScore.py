from sklearn import svm
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.model_selection import cross_val_score
from imblearn.over_sampling import SMOTE
from collections import Counter
from imblearn.over_sampling import RandomOverSampler
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition
from sklearn import datasets
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from imblearn.combine import SMOTEENN
from sklearn.neighbors import KNeighborsClassifier
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition
url_class = "prad_tcga_clinical_data.xlsx"
url_features = "prad_tcga_genes.xlsx"
#code to read the datasets
dataset=pd.read_excel(url_features)
dataset2=pd.read_excel(url_class)

#code to delete all the 0 values from the gene rows
data = dataset [~( dataset.loc[:,dataset.columns!='ID'] == 0).all(axis=1)]
X=np.transpose(data) #to transpose the matrix
y_gleason=dataset2['GLEASON_SCORE'] # to read the Gleason score column
Counter(y_gleason)

#code to perform classifcation without feature selection

#SVM-RBF 10 fold cross validation classification
model = svm.SVC(kernel='rbf',gamma='auto')
scores=cross_val_score(model, X.iloc[1:,1:], y_gleason, scoring='accuracy', cv=10)
print("Accuracy scores for each fold for SVM_RBF on the original data are :" , scores)
print("Mean Accuracy score :" , scores.mean())

#Random Forest Classifier 10 fold cross validation classification
model = RandomForestClassifier(n_estimators=1000)
scores=cross_val_score(model, X.iloc[1:,1:], y_gleason, scoring='accuracy', cv=10)
print("Accuracy scores for each fold for Random Forest algorithm on the original data are :" , scores)
print("Mean Accuracy score :" , scores.mean())

#KNN 10 fold classification
k_scores = []

# Calculating best values for K values between 1 and 20
for i in range(1, 20):
    knn = KNeighborsClassifier(n_neighbors=i)
    scores = cross_val_score(knn, X.iloc[1:,1:], y_gleason, cv=10, scoring='accuracy')
    k_scores.append(scores.mean())

plt.figure(figsize=(12, 6))
plt.plot(range(1, 20), k_scores, color='red', linestyle='dashed', marker='o',
         markerfacecolor='blue', markersize=10)
plt.title('Accuracy rate K value after running K-NN algorithm on the original Data')
plt.xlabel('K Value')
plt.ylabel('Cross-validated accuracy')

#we found best value of K = 11 so we used that for Classification
classifier = KNeighborsClassifier(n_neighbors=11)
scores = cross_val_score(classifier, X.iloc[1:,1:], y_gleason, cv=10, scoring='accuracy')
print("Accuracy scores for each fold for K-NN algorithm on original data are :" , scores)
print("Mean Accuracy score :" , scores.mean())

#code to plot 3-D graph using PCA on the orignal dataset
fig = plt.figure(1, figsize=(10, 10))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

plt.cla()
pca = decomposition.PCA(n_components=3)
pca.fit(X.iloc[1:,1:])
X_pca = pca.transform(X.iloc[1:,1:])

for i in range(len(y_gleason)):
    if y_gleason[i]==6:
        six=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='y',label=y_gleason[i])
    elif y_gleason[i]==7:
        seven=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='g',label=y_gleason[i])
    elif y_gleason[i]==8:
        eight=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='b',label=y_gleason[i])
    elif y_gleason[i]==9:
        nine=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='r',label=y_gleason[i])
    elif y_gleason[i]==10:
        ten=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='m',label=y_gleason[i])
plt.legend((six,seven,eight,nine,ten),
           ('Six', 'Seven', 'Eight','Nine','Ten'),
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=20)
plt.title('3-D Visualization of the original data without Feature Selection')
plt.show()


#code to perform feature selection

clf = RandomForestClassifier(n_estimators=1)
clf = clf.fit(X.iloc[1:,1:], y_gleason)
clf.feature_importances_

model = SelectFromModel(clf, prefit=True)
X_new_gleason = model.transform(X.iloc[1:,1:])
X_new_gleason.shape

#Code to Over Sample the data using SMOTE algorithm
sm = SMOTE(random_state=42,k_neighbors=3, m_neighbors='deprecated')
X_res, y_res = sm.fit_resample(X_new_gleason, y_gleason)
print('Resampled dataset shape after running SMOTE algorithm for Data Oversampling %s' % Counter(y_res))

#Code to Over Sample the data using Random Sampler algorithm
#ros = RandomOverSampler(random_state=42)
#X_res, y_res = ros.fit_resample(X_new, y)
#print('Resampled dataset shape %s' % Counter(y_res))
#X_res.shape
#y_res.shape
#code to plot 3-D graph using PCA after doing over sampling
fig = plt.figure(1, figsize=(10, 10))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

plt.cla()
pca = decomposition.PCA(n_components=3)
pca.fit(X_res)
X_pca = pca.transform(X_res)

for i in range(len(y_res)):
    if y_res[i]==6:
        six=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='y',label=y_res[i])
    elif y_res[i]==7:
        seven=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='g',label=y_res[i])
    elif y_res[i]==8:
        eight=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='b',label=y_res[i])
    elif y_res[i]==9:
        nine=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='r',label=y_res[i])
    elif y_res[i]==10:
        ten=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='m',label=y_res[i])
plt.legend((six,seven,eight,nine,ten),
           ('Six', 'Seven', 'Eight','Nine','Ten'),
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=20)
plt.title('3-D Visualization of the Oversampled Data')
plt.show()

#code to perform combination of over and under sampling using SMOTEENN
sme = SMOTEENN(random_state=42,smote=SMOTE(random_state=42, k_neighbors=3, m_neighbors='deprecated'))
X_res_smoteenn, y_res_smoteenn = sme.fit_resample(X_new_gleason, y_gleason)
print('Resampled dataset shape after running SMOTEENN algorithm for combination of Data Oversampling and Undersampling %s' % Counter(y_res_smoteenn))

#code to plot 3-D graph using PCA after doing over sampling and under sampling using SMOTEENN
fig = plt.figure(1, figsize=(10, 10))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

plt.cla()
pca = decomposition.PCA(n_components=3)
pca.fit(X_res_smoteenn)
X_pca = pca.transform(X_res_smoteenn)

for i in range(len(y_res_smoteenn)):
    if y_res_smoteenn[i]==6:
        six=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='y',label=y_res_smoteenn[i])
    elif y_res_smoteenn[i]==7:
        seven=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='g',label=y_res_smoteenn[i])
    elif y_res_smoteenn[i]==8:
        eight=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='b',label=y_res_smoteenn[i])
    elif y_res_smoteenn[i]==9:
        nine=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='r',label=y_res_smoteenn[i])
    elif y_res_smoteenn[i]==10:
        ten=ax.scatter(X_pca[i][0], X_pca[i][1], X_pca[i][2], c='m',label=y_res_smoteenn[i])
plt.legend((six,seven,eight,nine,ten),
           ('Six', 'Seven', 'Eight','Nine','Ten'),
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=20)
plt.title('3-D Visualization after running SMOTEENN algorithm to do combination of Data Oversampling and Undersampling')
plt.show()

#code to perform classification of over sampled data using SMOTE

#SVM-RBF 10 fold cross validation classification
model = svm.SVC(kernel='rbf',gamma='auto')
scores=cross_val_score(model, X_res, y_res, scoring='accuracy', cv=10)
print("Accuracy scores for each fold for SVM-RBF algorithm after doing Data Oversampling using SMOTE algorithm are :" , scores)
print("Mean Accuracy score :" , scores.mean())

#Random Forest Classifier 10 fold cross validation classification
model = RandomForestClassifier(n_estimators=1000)
scores=cross_val_score(model, X_res, y_res, scoring='accuracy', cv=10)
print("Accuracy scores for each fold for Random Forest algorithm after doing Data Oversampling using SMOTE algorithm are :" , scores)
print("Mean Accuracy score :" , scores.mean())

#KNN 10 fold classification
k_scores = []

# Calculating best values for K values between 1 and 20
for i in range(1, 20):
    knn = KNeighborsClassifier(n_neighbors=i)
    scores = cross_val_score(knn, X_res, y_res, cv=10, scoring='accuracy')
    k_scores.append(scores.mean())

plt.figure(figsize=(12, 6))
plt.plot(range(1, 20), k_scores, color='red', linestyle='dashed', marker='o',
         markerfacecolor='blue', markersize=10)
plt.title('Accuracy rate K value of K-NN algorithm after running Combination of Data Oversampling and Undersampling')
plt.xlabel('K Value')
plt.ylabel('Cross-validated accuracy')

#we found best alue of K = 1 so we used that for Classification
classifier = KNeighborsClassifier(n_neighbors=1)
scores = cross_val_score(classifier, X_res, y_res, cv=10, scoring='accuracy')
print("Accuracy scores for each fold for K-NN algorithm after doing Data Oversampling using SMOTE algorithmare :" , scores)
print("Mean Accuracy score :" , scores.mean())

#code to perform classification of over sampled data using SMOTEENN

#SVM-RBF 10 fold cross validation classification
model = svm.SVC(kernel='rbf',gamma='auto')
scores=cross_val_score(model, X_res_smoteenn, y_res_smoteenn, scoring='accuracy', cv=10)
print("Accuracy scores for each fold for SVM-RBF algorithm after doing Combination of Data Oversampling and Undersampling using SMOTEENN algorithm are :" , scores)
print("Mean Accuracy score :" , scores.mean())

#Random Forest Classifier 10 fold cross validation classification
model = RandomForestClassifier(n_estimators=1000)
scores=cross_val_score(model, X_res_smoteenn, y_res_smoteenn, scoring='accuracy', cv=10)
print("Accuracy scores for each fold for Random Forest algorithm after doing Combination of Data Oversampling and Undersampling using SMOTEENN algorithm are :" , scores)
print("Mean Accuracy score :" , scores.mean())

#KNN 10 fold classification
k_scores = []

# Calculating best values for K values between 1 and 40
for i in range(1, 20):
    knn = KNeighborsClassifier(n_neighbors=i)
    scores = cross_val_score(knn, X_res_smoteenn, y_res_smoteenn, cv=10, scoring='accuracy')
    k_scores.append(scores.mean())

plt.figure(figsize=(12, 6))
plt.plot(range(1, 20), k_scores, color='red', linestyle='dashed', marker='o',
         markerfacecolor='blue', markersize=10)
plt.title('Accuracy rate K value after doing Combination of Data Oversampling and Undersampling using SMOTEENN algorithm ')
plt.xlabel('K Value')
plt.ylabel('Cross-validated accuracy')

#we found best value of K = 1 so we used that for Classification
classifier = KNeighborsClassifier(n_neighbors=1)
scores = cross_val_score(classifier, X_res_smoteenn, y_res_smoteenn, cv=10, scoring='accuracy')
print("Accuracy scores for each fold for K-NN algorithm after doing Combination of Data Oversampling and Undersampling using SMOTEENN algorithm are :" , scores)
print("Mean Accuracy score :" , scores.mean())

#Code to give list of features
feat_labels=X.iloc[0]
# Print the name and gini importance of each feature
for feature in zip(feat_labels, clf.feature_importances_):
    if(feature[1]!=0):
        print(feature)
