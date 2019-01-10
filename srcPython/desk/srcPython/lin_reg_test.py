import numpy as np

y = [1,2,3,4,3,4,5,4,5,5,4,5,4,5,4,5,6,5,4,5,4,3,4]
x =     [[4,2,3,4,5,4,5,6,7,4,8,9,8,8,6,6,5,5,5,5,5,5,5],
     [4,1,2,3,4,5,6,7,5,8,7,8,7,8,7,8,7,7,7,7,7,6,5],
     [4,1,2,5,6,7,8,9,7,8,7,8,7,7,7,7,7,7,6,6,4,4,4]
     ]

from sklearn import linear_model
clf = linear_model.LinearRegression()
arr_x = np.array(x)
clf.fit(arr_x.transpose(),
         y)
coeff = clf.coef_
const = clf.intercept_
y_model = np.zeros([len(y)])
for i in range(len(y)):
    y_model[i] = np.sum(arr_x.transpose()[i]*coeff) + const
print coeff
print const
#print y - y_model
print 'score', clf.score(arr_x.transpose(),
         y)
print 'predict', clf.predict(arr_x.transpose())
print 'ymodel', y_model
