import matplotlib.pyplot as plt
import seaborn as sn
import math
import numpy as np
def pred_vs_labels(preds, labels, text):

        file = 'graphs/predictions_{}.txt'.format(text)  #('_pdb' if self.test_set == 'pdb' else ''))
        correlation = np.corrcoef(preds,labels)[0,1]


        
        sn.set(style="darkgrid")

        rmse = math.sqrt(((preds - labels)**2).mean())
        plt.figure()
        h = sn.jointplot(
            x=labels,
            y=preds,
            xlim=(1,13),
            ylim=(1,13),
            kind='reg',
            joint_kws={"scatter_kws":{'alpha':0.2, 'edgecolors':'white'}, "label":f"correlation: {correlation}",} )
        h.set_axis_labels('Experimental pK', 'Predicted pK')
        h.fig.suptitle(f"Test Predictions - R: {correlation:.2f}, RMSE: {rmse:.2f}")
       
        plt.tight_layout()

        # h.fig.plot([0,1],[0,1],  'r', transform=h.fig.transAxes)
