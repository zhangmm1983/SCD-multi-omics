import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull

df = pd.read_excel("模型参数.xlsx",index_col=0)
df_D = df[df.group.str.contains("D")]
df_S = df[df.group.str.contains("S")]

print(df_D)
print(df_S)

##################################
fig = plt.figure(figsize=(10,8))
ax = fig.add_axes([0.2,0.2,0.6,0.6])
# plt.xlim(min(df_D.PC1)*2,max(df_D.PC1)*2)
# plt.ylim(min(df_D.PC2)*2,max(df_D.PC2)*2)


def draw(df,c,shape):
    points = np.asarray([[df.PC1.tolist()[i],df.PC2.tolist()[i]] for i in range(len(df.PC1))])

    hull = ConvexHull(points)  ###计算外接凸图案的边界点
    hull1=hull.vertices.tolist()
    hull1.append(hull1[0])
    ax.plot(points[hull1,0], points[hull1,1],lw=2,color=c,alpha=0.3)
    ax.fill_between(points[hull1,0],points[hull1,1],alpha=0.1,facecolor=c)
    ax.scatter(df.PC1,df.PC2,color=c,marker=shape)

draw(df_D,"#941225",'o')
draw(df_S,"#226E4B","^")


ax.spines['bottom'].set_linewidth(3)
ax.spines['left'].set_linewidth(3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xticks(size=15)
plt.yticks(size=15)
plt.xlabel("PC1(0.0225%)",size=25)
plt.ylabel("PC01(0.305%)",size=25)
plt.title("Metabolites-LC",size=25)

plt.xlim([-130,130])

plt.savefig("Metabolites-LC_opls.pdf")
plt.savefig("Metabolites-LC_opls.png")