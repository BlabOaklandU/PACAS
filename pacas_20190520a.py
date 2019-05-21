#c:\Pacas\test_20190511a.py
#Kenneth

#import matplotlib



import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd


#from mpl_toolkits.mplot3d import Axes3D
#from scipy.optimize import curve_fit


DataFrame = pd.read_csv('pacas.csv', skiprows=3, header=0)

print("\n\nFirst five rows\n\n")

print(DataFrame.head())

print("\n\nLast five rows\n\n")

print(DataFrame.tail())

#x = DataFrame.sum(axis = 'columns')
#print("Sums =\n\n", x)






#y = DataFrame.describe()
#print("Describe =\n\n", y)



DataFrame["Accession_Combined"] = DataFrame["Accession #1"] + ' : ' + DataFrame["Accession #2"]

for col in DataFrame.columns:
    print(col)

#hist = DataFrame.hist(bins=5)

#ax = plt.gca()
#z = DataFrame.plot(kind='line',x='Accession #1', y='Comparison D Mismatches')

labels = list(DataFrame.Accession_Combined)

print("labels = \n", labels)

#fig = plt.figure()

#ax1 = fig.add_subplot(221)
#ax2 = fig.add_subplot(222)
#ax3 = fig.add_subplot(223)
#ax4 = fig.add_subplot(224)

print(DataFrame.head( ))

DataFrame.sort_values('Comparison A Mismatches', inplace=True)

print("\n\nSorting: \n\n")


print(DataFrame.head())



ax1 = DataFrame.plot(x='Accession_Combined',   y='Comparison A Mismatches',   kind = 'barh')
#plt.xlabel("labels")
#c='DarkBlue',
#xticks = 'labels'

#ax1.set_yticklabels(labels, rotation=0)

#ax2 = DataFrame.plot(x='Accession #1',  y='Comparison B Mismatches',  c='DarkGreen')

#ax3 = DataFrame.plot(x='Accession #1',  y='Comparison C Mismatches',  c='DarkRed')

#ax4 = DataFrame.plot(x='Accession #1',  y='Comparison D Mismatches',  c='Yellow')

#Accession_Combined

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#n = 200

#xs = np.random.rand(n)
#ys = np.random.rand(n)
#zs = np.random.rand(n)
#rs = np.sqrt(xs*xs + ys*ys + zs*zs)

#ax.scatter(xs, ys, zs, c=rs, marker='o')

#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')

#plt.show()





#func = lambda x, a, b, c: a + b*x * np.exp(-c*x)

#x = np.linspace(0, 6, 20)
#noise = 0.05 * np.random.normal(size=len(x))
#y = func(x,1,2,1) + noise

#sigma = np.ones(len(x))
#p0 = np.ones(3)          # initial guesses for a, b, c

#p, pcov = curve_fit(func,x,y,p0,sigma)

#perr = np.sqrt(np.diag(pcov))

#print ('Best-fit:')
#print u'a = {:g} \u00B1 {:g}'.format(p[0],perr[0])
#print u'b = {:g} \u00B1 {:g}'.format(p[1],perr[1])
#print u'c = {:g} \u00B1 {:g}'.format(p[2],perr[2])

#xfit = np.linspace(0,6,100)

#plt.plot(x,y,'ko',label='Data')
#plt.plot(xfit, func(xfit,p[0],p[1],p[2]),'b-',label='Fit')
#plt.legend()
#plt.text(4,1.5,"$y = a + bx e^{-cx}$",fontsize=20)
#plt.axis([-0.5,6.5,0.92,1.85])
plt.show()
