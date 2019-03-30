import numpy as np
from sklearn.neighbors import KNeighborsClassifier as KNN

X_train = np.genfromtxt('Spectogram_train_1.csv',delimiter = ',')
X_test = np.genfromtxt('Spectogram_test_1.csv',delimiter = ',')
#X_train = np.genfromtxt('Spectogram_train_2.csv',delimiter = ',')
#X_test = np.genfromtxt('Spectogram_test_2.csv',delimiter = ',')
#X_train = np.genfromtxt('Spectogram_train_3.csv',delimiter = ',')
#X_test = np.genfromtxt('Spectogram_test_3.csv',delimiter = ',')
print('Done loading data')

y_train = np.zeros((270,1));
y_test = np.zeros((27,1));

for i in range(3):
	for j in range(9):
		if j < 3:
			y_test[9*i + j,0] = 1
		elif j < 6:
			y_test[9*i + j,0] = 2
		else:
			y_test[9*i + j,0] = 3


for j in range(270):
	if j < 90:
		y_train[j,0] = 1
	elif j < 180:
		y_train[j,0] = 2
	else:
		y_train[j,0] = 3
	

classifier = KNN(n_neighbors = 5, weights = 'distance', metric = 'euclidean')

classifier.fit(X_train, y_train.ravel())
train_accuracy = classifier.score(X_train, y_train.ravel())
print('Train accuracy: ' + str(train_accuracy))

predictions_for_test = classifier.predict(X_test)
print(predictions_for_test)
test_accuracy = classifier.score(X_test, y_test.ravel())
print('Test accuracy: ' + str(test_accuracy))
