'''Plot grafic 4a'''

import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns

path = str(pathlib.Path(__file__).parent.resolve())
df = pd.read_excel(path + '/fig4a.xlsx')
# df = df[df['density']<10]
g = sns.regplot(data=df,x="count", y="length",scatter=True, color = '#2D708EFF')
plt.ylim(0,50)
sns.regplot(data=df,x="count", y="density",scatter=True, ax=g.axes.twinx(), color='#404788FF')
g.set_xlabel('Seminal root count')
g.set_ylabel('Primary root length (cm)', color = '#2D708EFF')
plt.ylabel('Lateral root density (no./cm)', color = '#404788FF')
g.text(2, 45, '$\it{r} = -0.2, \it{p} = 0.003$', fontsize=11)
g.text(2, 2, '$\it{r} = -0.16, \it{p} = 0.017$', fontsize=11)
plt.xlim(0,10)

plt.show()