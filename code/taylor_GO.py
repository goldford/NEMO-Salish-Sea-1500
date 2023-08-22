import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import skill_metrics as sm

# Generate some example data
observed = np.random.rand(50)
model1 = np.round(np.random.rand(50)+100,2)
model2 = np.round(np.random.rand(50)+100,2)
model3 = np.round(np.random.rand(50)+100,2)
model4 = np.round(np.random.rand(50)+100,2)

taylor_stats1 = sm.taylor_statistics(model1, observed, 'data')
taylor_stats2 = sm.taylor_statistics(model2, observed, 'data')
taylor_stats3 = sm.taylor_statistics(model3, observed, 'data')
taylor_stats4 = sm.taylor_statistics(model4, observed, 'data')

sdev = np.array([taylor_stats1['sdev'][0], taylor_stats1['sdev'][1],
                 taylor_stats2['sdev'][1], taylor_stats3['sdev'][1]])
crmsd = np.array([taylor_stats1['crmsd'][0], taylor_stats1['crmsd'][1],
                  taylor_stats2['crmsd'][1], taylor_stats3['crmsd'][1]])
ccoef = np.array([taylor_stats1['ccoef'][0], taylor_stats1['ccoef'][1],
                  taylor_stats2['ccoef'][1], taylor_stats3['ccoef'][1]])

# Create a figure and a 2x2 grid of subplots using GridSpec
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 2, figure=fig)

# Create Taylor diagrams for each subplot
ax0 = fig.add_subplot(gs[0, 0])
ax0 = sm.taylor_diagram(sdev, crmsd, ccoef, numberPanels=1,
                        styleCOR='-', colCOR='#585858',
                        colRMS='#585858', titleRMS='off', labelRMS='NCRMSD', rmsLabelFormat='0:.1f')

ax1 = fig.add_subplot(gs[0, 1])
ax1 = sm.taylor_diagram(sdev, crmsd, ccoef, numberPanels=1, styleCOR='-', colCOR='#585858',
                        colRMS='#585858', labelRMS='NCRMSD', rmsLabelFormat='0:.1f')
ax2 = fig.add_subplot(gs[1, 0])
ax2 = sm.taylor_diagram(sdev, crmsd, ccoef, numberPanels=1, styleCOR='-', colCOR='#585858',
                        colRMS='#585858', labelRMS='NCRMSD', rmsLabelFormat='0:.1f')
ax3 = fig.add_subplot(gs[1, 1])
ax3 = sm.taylor_diagram(sdev, crmsd, ccoef, numberPanels=1, styleCOR='-', colCOR='#585858',
                        colRMS='#585858', labelRMS='NCRMSD', rmsLabelFormat='0:.1f')

# Set titles and adjust layout
# axes[0].set_title('Model 1')
# axes[1].set_title('Model 2')
# axes[2].set_title('Model 3')
# axes[3].set_title('Model 4')

plt.tight_layout()
plt.show()
