## detect_puffs ##
detect\_puffs is a plugin for the image processing program [Flika](http://brettjsettle.pythonanywhere.com/).  It contains a set of algorithms for detecting local calcium signals from movies.

### Tutorials ###
[How to detect subcellular Ca2+ signals using Flika](
http://htmlpreview.github.io/?https://github.com/kyleellefsen/detect_puffs/blob/master/docs/How%20to%20detect%20subcellular%20Ca2%2B%20signals%20using%20Flika.html)

### Sample Code ###
```python
from plugins.detect_puffs.puff_simulator.puff_simulator import simulate_puffs
from plugins.detect_puffs.threshold_cluster import threshold_cluster
data_window=simulate_puffs(nFrames=1000,nPuffs=20)
data_window=g.m.currentWindow
norm_window=data_window
norm_window.setName('Normalized Window')
binary_window=threshold(1.1, keepSourceWindow=True)
binary_window.setName('Binary Window')
threshold_cluster(binary_window,data_window,norm_window,density_threshold=5.8)
```
